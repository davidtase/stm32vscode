// Decoder of a saved ECG + PCG recording for viewer.html, with no DOM, so node can test it. The file format is the one in
// recorder.js and Tools/load_recording.py. Raw bytes are not kept, except the PCG's (1.5 MB per 10 minutes), which are walked as the loader does.
// Time is PCG time: milliseconds = PCG sample index, index 0 = the first live frame's first sample, as in the loader's pcg.npy.
import { buildMap, filtfilt, findBeats, alignCheck, HP05 } from './align.js';
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
export const HP = [[0.997781024102941, -1.99556204820588, 0.997781024102941], [1, -1.99555712434579, 0.995566972065975]];
export const NF = [
  [0.951558090821252, -2.1284207536999, 2.85527558178394, -2.1284207536999, 0.951558090821252],
  [1, -2.18192506761969, 2.85407296314357, -2.07491643978012, 0.904318800282871],
];
export const HF = [
  [0.987761389276826, -3.9510455571073, 5.92656833566095, -3.9510455571073, 0.987761389276826],
  [1, -3.97537191256092, 5.92641855596542, -3.92671919775679, 0.975672562146085],
];
const H40 = [
  -0.00812031626499363, 0.00279430037155961, 0.0363062486230826, 0.0960723789024121, 0.163761187856556, 0.209186200511384,
  0.209186200511384, 0.163761187856556, 0.0960723789024121, 0.0363062486230826, 0.00279430037155961, -0.00812031626499363,
];
export const BLOCK_S = 30; // seconds per block of the ECG-to-PCG map
const PAD = 15; // filtfilt pads 3 * the longest coefficient list (5) on both ends: a run of PAD frames or fewer cannot be filtered

// One run of raw lead II counts to 2n - 1 samples at 1000 Hz, sample 2m being frame m. The IIR stages go forward and backward;
// the 12-tap FIR is symmetric, so one causal pass and a shift of its 5.5 samples = 11 samples at 1000 Hz centre it. Its ends repeat the edge values.
// stages: the viewer's trace by default, the labeling trace (one 0.5 Hz high-pass in place of the 0.25 and 0.75 Hz ones) for the detectors
export function zeroPhase(x, stages = [HP, NF, HF]) {
  let y = Float64Array.from(x);
  for (const stage of stages) y = filtfilt(stage, y);
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
// b is ms per sample. Not a correction of the recording: it only tells the host clock from the PCG one
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

// The PCG frames as Tools/load_recording.py finds them: 24 bytes [0x8000][seq][10 int16], a frame counts when the header of the next
// one follows 24 bytes later, the first (stale) one is left out. Each frame is stamped with the chunk that held its last byte, and a
// dropped frame stays in the array as 10 NaN, so an index is a time. ends: cumulative byte count at the end of each chunk, ts: its stamp
function pcgFrames(b, ends, ts) {
  const ok = (c) => b[c] === 0 && b[c + 1] === 0x80 && b[c + 3] < 0x80, seq = (c) => b[c + 2] | (b[c + 3] << 8), starts = [];
  for (let c = 0; c + 28 <= b.length; c++) {
    if (ok(c) && ok(c + 24) && seq(c) !== seq(c + 24)) starts.push(c), c += 23;
  }
  const live = starts.slice(1), index = new Float64Array(live.length), t = new Float64Array(live.length);
  let n = 0, drops = 0, k = 0;
  live.forEach((c, j) => {
    if (j) {
      const gap = (seq(c) - seq(live[j - 1]) - 1) & 0x7fff;
      drops += gap;
      n += gap + 1;
    }
    index[j] = 10 * n;
    while (ends[k] < c + 24) k++;
    t[j] = ts[k];
  });
  const dv = new DataView(b.buffer, b.byteOffset, b.byteLength), pcg = new Float32Array(live.length ? 10 * (n + 1) : 0).fill(NaN);
  live.forEach((c, j) => { for (let i = 0; i < 10; i++) pcg[index[j] + i] = dv.getInt16(c + 4 + 2 * i, true); });
  return { pcg, index, t, drops };
}

// chunks: async iterable of ArrayBuffer, one per segment in number order. Never throws on a damaged recording
export async function loadRecording(chunks) {
  const errors = [], events = [], pcgBytes = sink(Uint8Array), ecgSink = sink(Float32Array), rawSink = sink(Int16Array), stampSink = sink(Float64Array);
  const ecgFit = fitter(), ecgParser = createEcgParser(), pcgEnds = [], pcgTs = [];
  const ecgResyncs = []; // frame numbers where a run of raw ECG frames starts after a resync
  let segments = 0, first = null, startT = null, resync = false;

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
        pcgBytes.push(new Uint8Array(buf, at, n));
        pcgEnds.push(pcgBytes.length);
        pcgTs.push(t);
      } else if (kind === ECG) {
        const { stats } = ecgParser, [locked, skipped] = [stats.frames > 0, stats.skipped];
        const frames = ecgParser.push(new Uint8Array(buf, at, n));
        // shortcut: the parser does not say where in a chunk it skipped, so a cut can be one chunk (a few frames) early, upgrade only if a recording with resyncs matters
        if (locked && stats.skipped > skipped) resync = true;
        if (frames.length) {
          if (resync) ecgResyncs.push(rawSink.length);
          resync = false;
          rawSink.push(frames.map((f) => f[1])); // lead II
          stampSink.push(new Float64Array(frames.length).fill(t)); // each frame: the chunk that held its last byte
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
  let ecg = ecgSink.done(), ecgLabel = null;
  if (raw.length) {
    ecg = new Float32Array(2 * raw.length - 1).fill(NaN);
    ecgLabel = ecg.slice();
    [0, ...ecgResyncs].forEach((s, r) => {
      const x = raw.subarray(s, ecgResyncs[r] ?? raw.length);
      if (x.length > PAD) [ecg.set(zeroPhase(x), 2 * s), ecgLabel.set(zeroPhase(x, [HP05, NF]), 2 * s)];
      else ecgShort++;
    });
  }
  const { pcg, index, t: pcgT, drops } = pcgFrames(pcgBytes.done(), pcgEnds, pcgTs), pcgLine = fitter();
  index.forEach((i, j) => pcgLine.add(i, pcgT[j]));
  const pcgFit = pcgLine.done(), line = pcgFit ?? { i0: 0, t0: first ?? 0, a: 0, b: 1, hz: 1000 }; // no PCG: time is ms since the first record
  const stamps = stampSink.done(), map = raw.length ? buildMap(stamps, ecgResyncs, line, BLOCK_S) : null, stored = raw.length ? null : ecgFit.done();
  // R and T-end proposals on the labeling trace, and the check of the map against the PCG's S1
  const beats = map ? findBeats(raw, ecgLabel, [0, ...ecgResyncs], map) : null, check = beats?.length ? alignCheck(pcg, beats.map((b) => 1000 * b.t)) : null;
  // Both traces on PCG time, in ms (the PCG's is its own index). The stored forward-only trace goes through the PCG line instead of the map
  const toPcgMs = (host) => (host - line.t0 - line.a) / line.b + line.i0, toHost = (t) => line.t0 + line.a + line.b * (t - line.i0);
  const ecgAxis = map ? { timeAt: (i) => 1000 * map.timeOf(i / 2), indexAt: (t) => 2 * map.frameAt(t / 1000) }
    : stored && { timeAt: (i) => toPcgMs(stored.t0 + stored.a + stored.b * (i - stored.i0)), indexAt: (t) => stored.i0 + (toHost(t) - stored.t0 - stored.a) / stored.b };
  return {
    pcg, ecg, ecgLabel, beats, check, pcgFit, pcgFrameT: pcgT, pcgFrameIndex: index, ecgStamps: stamps, map, ecgHz: map ? 2 * map.frameHz : stored?.hz,
    axes: { pcg: { timeAt: (i) => i, indexAt: (t) => t }, ecg: ecgAxis }, toPcgMs,
    ecgSource: raw.length ? 'raw' : 'stored', ecgResyncs, ecgShort, t0: startT ?? first, events, segments, drops, errors,
  };
}

// For each of width pixel columns between PCG times t0 and t1 (ms): the min and max of the samples that fall in it, NaN where none does.
// axis.indexAt(t) is the fractional index of the samples at time t. A single-sample spike is therefore always the max of its column, at any zoom
export function columns(samples, axis, t0, t1, width) {
  const min = new Float32Array(width).fill(NaN), max = new Float32Array(width).fill(NaN);
  let lo = Math.max(0, Math.ceil(axis.indexAt(t0)));
  for (let x = 0; x < width; x++) {
    const hi = Math.min(samples.length, Math.ceil(axis.indexAt(t0 + ((x + 1) * (t1 - t0)) / width)));
    if (hi > lo) {
      let a = Infinity, b = -Infinity; // NaN samples (a dropped frame, the join of two ECG runs) are skipped
      for (let i = lo; i < hi; i++) {
        const v = samples[i];
        if (v < a) a = v;
        if (v > b) b = v;
      }
      if (a <= b) [min[x], max[x]] = [a, b];
    }
    lo = Math.max(lo, hi);
  }
  return { min, max };
}
