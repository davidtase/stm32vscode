// Decoder of a saved ECG + PCG recording for viewer.html, with no DOM, so node can test it. The file format is the one in
// recorder.js and Tools/load_recording.py. Raw bytes are not kept: only the decoded samples, which are 28 times smaller.
import { createParser } from './pcg_parser.js';
import { createParser as createEcgParser } from './ecg_dsp.js';
import { createUpsampler } from './upsample.js';

const PCG = 1, ECG = 2, EVENT = 3, ECG_1K = 4; // same numbers as recorder.js
const HEADER = 13;
const SHOWN = new Set(['start', 'stop', 'lost', 'read-error']);

// names: the file names of the chosen folder. Returns the indices of the numbered segments in number order, and what is wrong
export function segmentList(names) {
  const errors = [], found = [];
  names.forEach((name, i) => {
    const m = /^(.*)\.ecgpcg$/.exec(name);
    if (!m) return;
    if (/^\d+$/.test(m[1])) found.push({ i, k: Number(m[1]) });
    else errors.push(`${name}: not a numbered segment`);
  });
  found.sort((a, b) => a.k - b.k);
  if (!found.length) errors.push('no segment files');
  const have = new Set(found.map((f) => f.k));
  const missing = [];
  for (let k = 1; k <= (found.at(-1)?.k ?? 0); k++) if (!have.has(k)) missing.push(k);
  if (missing.length) errors.push(`missing segment numbers: ${missing.join(', ')}`);
  return { order: found.map((f) => f.i), errors };
}

// A growing typed array, so the sample count need not be known in advance
function sink(Type) {
  let a = new Type(1 << 16), n = 0;
  return {
    push(s) {
      if (n + s.length > a.length) {
        const b = new Type(Math.max(2 * a.length, n + s.length));
        b.set(a.subarray(0, n));
        a = b;
      }
      a.set(s, n);
      n += s.length;
    },
    get length() { return n; },
    done: () => a.slice(0, n),
  };
}

// The recorder's ECG filters (protocol section 12.2, [b, a] at 500 Hz), copied because ecg_dsp.js keeps them private. The tests check them
const HP = [[0.997781024102941, -1.99556204820588, 0.997781024102941], [1, -1.99555712434579, 0.995566972065975]];
const NF = [
  [0.951558090821252, -2.1284207536999, 2.85527558178394, -2.1284207536999, 0.951558090821252],
  [1, -2.18192506761969, 2.85407296314357, -2.07491643978012, 0.904318800282871],
];
const HF = [
  [0.987761389276826, -3.9510455571073, 5.92656833566095, -3.9510455571073, 0.987761389276826],
  [1, -3.97537191256092, 5.92641855596542, -3.92671919775679, 0.975672562146085],
];
const H40 = [
  -0.00812031626499363, 0.00279430037155961, 0.0363062486230826, 0.0960723789024121, 0.163761187856556, 0.209186200511384,
  0.209186200511384, 0.163761187856556, 0.0960723789024121, 0.0363062486230826, 0.00279430037155961, -0.00812031626499363,
];
const PAD = 15; // filtfilt pads 3 * the longest coefficient list (5) on both ends: a run of PAD frames or fewer cannot be filtered

// scipy.signal.lfilter, direct form II transposed, from state z (changed)
function lfilter([b, a], x, z) {
  const n = z.length, y = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) {
    const o = b[0] * x[i] + z[0];
    for (let k = 0; k < n; k++) z[k] = (k + 1 < n ? z[k + 1] : 0) + b[k + 1] * x[i] - a[k + 1] * o;
    y[i] = o;
  }
  return y;
}

// scipy.signal.filtfilt with its defaults (odd extension of 3 * len(b), start from the steady state), on a run of more than 3 * len(b) samples
function filtfilt(ba, x) {
  const [b, a] = ba, n = x.length, pad = 3 * b.length, g = b.reduce((s, v) => s + v) / a.reduce((s, v) => s + v);
  const zi = new Float64Array(b.length - 1); // lfilter_zi: the state of a constant input of 1
  for (let k = b.length - 2, s = 0; k >= 0; k--) zi[k] = s += b[k + 1] - a[k + 1] * g;
  const e = new Float64Array(n + 2 * pad);
  for (let i = 0; i < pad; i++) {
    e[i] = 2 * x[0] - x[pad - i];
    e[pad + n + i] = 2 * x[n - 1] - x[n - 2 - i];
  }
  e.set(x, pad);
  const pass = (s) => lfilter(ba, s, zi.map((v) => v * s[0])).reverse();
  return pass(pass(e)).subarray(pad, pad + n);
}

// One run of raw lead II counts to 2n - 1 samples at 1000 Hz, sample 2m being frame m. The three IIR stages go forward and backward;
// the 12-tap FIR is symmetric, so one causal pass and a shift of its 5.5 samples = 11 samples at 1000 Hz centre it. Its ends repeat the edge values
function zeroPhase(x) {
  let y = Float64Array.from(x);
  for (const stage of [HP, NF, HF]) y = filtfilt(stage, y);
  const M = y.length, up = createUpsampler(2), out = new Float32Array(2 * (M + 6) - 1);
  let o = 0;
  for (let m = 0; m < M + 6; m++) {
    let c = 0;
    for (let k = 0; k < 12; k++) c += H40[k] * y[Math.min(Math.max(m - k, 0), M - 1)];
    for (const v of up(c)) out[o++] = v;
  }
  return out.subarray(11, 11 + 2 * M - 1);
}

// Least squares line t = t0 + a + b * (index - i0) through (index, t) pairs, as running sums so no pair is kept.
// b is ms per sample. Not a correction of the recording: it only places the trace on the host clock for display
function fitter() {
  let n = 0, i0 = 0, t0 = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
  return {
    add(i, t) {
      if (!n) [i0, t0] = [i, t];
      const x = i - i0, y = t - t0;
      n++; sx += x; sy += y; sxx += x * x; sxy += x * y;
    },
    done() {
      const d = n * sxx - sx * sx;
      if (n < 2 || !d) return null;
      const b = (n * sxy - sx * sy) / d;
      return { i0, t0, a: (sy - b * sx) / n, b, hz: 1000 / b };
    },
  };
}

// chunks: async iterable of ArrayBuffer, one per segment in number order. Never throws on a damaged recording
export async function loadRecording(chunks) {
  const errors = [], events = [], pcgSink = sink(Int16Array), ecgSink = sink(Float32Array), rawSink = sink(Int16Array);
  const pcgFit = fitter(), ecgFit = fitter(), rawFit = fitter(), parser = createParser(), ecgParser = createEcgParser();
  const ecgResyncs = []; // frame numbers where a run of raw ECG frames starts after a resync
  let segments = 0, first = null, startT = null, stale = true, resync = false;

  for await (const buf of chunks) {
    const dv = new DataView(buf), name = `file ${++segments}`;
    for (let i = 0; i < buf.byteLength; ) {
      const n = i + HEADER <= buf.byteLength ? dv.getUint32(i + 8, true) : null;
      if (n === null || i + HEADER + n > buf.byteLength) {
        errors.push(`${name}: the record at byte ${i} overruns the file`);
        break;
      }
      const t = dv.getFloat64(i, true), kind = dv.getUint8(i + 12), at = i + HEADER;
      i = at + n;
      if (!Number.isFinite(t)) {
        errors.push(`${name}: the record before byte ${i} has no valid time`);
        continue;
      }
      first ??= t;
      if (kind === PCG) {
        let kept = 0;
        for (const f of parser.push(new Uint8Array(buf, at, n))) {
          if (stale) stale = false; // the device queues one stale frame while the port is closed
          else { pcgSink.push(f); kept++; }
        }
        if (kept) pcgFit.add(pcgSink.length - 1, t);
      } else if (kind === ECG) {
        const { stats } = ecgParser, [locked, skipped] = [stats.frames > 0, stats.skipped];
        const frames = ecgParser.push(new Uint8Array(buf, at, n));
        // shortcut: the parser does not say where in a chunk it skipped, so a cut can be one chunk (a few frames) early, upgrade only if a recording with resyncs matters
        if (locked && stats.skipped > skipped) resync = true;
        if (frames.length) {
          if (resync) ecgResyncs.push(rawSink.length);
          resync = false;
          rawSink.push(frames.map((f) => f[1])); // lead II
          rawFit.add(2 * (rawSink.length - 1), t); // in the 1000 Hz index space, where frame m is sample 2m
        }
      } else if (kind === ECG_1K) {
        if (n % 4) errors.push(`${name}: a filtered ECG record of ${n} bytes is not whole samples`);
        if (n >= 4) {
          ecgSink.push(new Float32Array(buf.slice(at, at + (n & ~3)))); // a copy: the payload is not 4-byte aligned
          ecgFit.add(ecgSink.length - 1, t);
        }
      } else if (kind === EVENT) {
        try {
          const e = JSON.parse(new TextDecoder().decode(new Uint8Array(buf, at, n)));
          if (e === null || typeof e !== 'object' || Array.isArray(e)) throw new Error('not an object');
          if (SHOWN.has(e.e)) events.push({ t, ...e });
          if (e.e === 'start') startT ??= t;
        } catch (err) {
          errors.push(`${name}: the event before byte ${i}: ${err.message}`);
        }
      } else {
        errors.push(`${name}: the record before byte ${i} has unknown kind ${kind}`);
      }
    }
    await null; // lets the page paint its progress line between files
  }
  // Raw frames give the ECG filtered offline, each run between resyncs on its own and a NaN between runs. The kind 4 trace is the
  // fallback for a recording with no raw frames
  const raw = rawSink.done();
  let ecgShort = 0; // runs too short to filter, left NaN
  let ecg = ecgSink.done();
  if (raw.length) {
    ecg = new Float32Array(2 * raw.length - 1).fill(NaN);
    [0, ...ecgResyncs].forEach((s, r) => {
      const x = raw.subarray(s, ecgResyncs[r] ?? raw.length);
      if (x.length > PAD) ecg.set(zeroPhase(x), 2 * s);
      else ecgShort++;
    });
  }
  return {
    pcg: pcgSink.done(), ecg, pcgFit: pcgFit.done(),
    ecgFit: raw.length ? rawFit.done() : ecgFit.done(), ecgSource: raw.length ? 'raw' : 'stored', ecgResyncs, ecgShort,
    t0: startT ?? first, events, segments, drops: parser.stats.drops, errors,
  };
}

// Host time in ms of sample i, and the other way round
export const timeAt = (fit, i) => fit.t0 + fit.a + fit.b * (i - fit.i0);
export const indexAt = (fit, t) => fit.i0 + (t - fit.t0 - fit.a) / fit.b;

// For each of width pixel columns between host times t0 and t1 (ms): the min and max of the samples that fall in it,
// NaN where none does. A single-sample spike is therefore always the max of its column, at any zoom
export function columns(samples, fit, t0, t1, width) {
  const min = new Float32Array(width).fill(NaN), max = new Float32Array(width).fill(NaN);
  let lo = Math.max(0, Math.ceil(indexAt(fit, t0)));
  for (let x = 0; x < width; x++) {
    const hi = Math.min(samples.length, Math.ceil(indexAt(fit, t0 + ((x + 1) * (t1 - t0)) / width)));
    if (hi > lo) {
      let a = samples[lo], b = a;
      for (let i = lo + 1; i < hi; i++) {
        const v = samples[i];
        if (v < a) a = v;
        else if (v > b) b = v;
      }
      min[x] = a;
      max[x] = b;
    }
    lo = Math.max(lo, hi);
  }
  return { min, max };
}
