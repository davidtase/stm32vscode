// Puts the ECG frames on the PCG clock and proposes R peaks and T-end points, with no DOM, so node can test it. The PCG counter is exact
// (0.6 ppm from the host clock) and the ECG module's clock is not, so each ECG frame's host stamp gets a straight line per run plus a slow
// wander read off the stamps. Display and labels only: nothing in the recording is corrected, and every proposal is the user's to review.
const FRAME_HZ = 500; // nominal ECG frame rate, used to size blocks and windows
const SHORT = 10 * FRAME_HZ; // a run under 10 s is too short for blocks: its line only

const median = (x) => {
  const s = Float64Array.from(x).sort(), h = s.length >> 1; // typed, so the sort is numeric
  return s.length % 2 ? s[h] : (s[h - 1] + s[h]) / 2;
};

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
export function filtfilt(ba, x) {
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

// stamps: host ms of each ECG frame (the chunk that held its last byte). runStarts: frame numbers where a run starts after a resync.
// pcgLine: {i0, t0, a, b}, host ms = t0 + a + b * (PCG index - i0). blockS: seconds per block of residual medians.
// Per run: the least squares line of stamp against frame number, plus the median of its residual in blocks, joined linearly between
// block centres and held at the ends. The median makes a stalled stamp harmless.
export function buildMap(stamps, runStarts, pcgLine, blockS) {
  const M = stamps.length, starts = [0, ...runStarts], B = Math.round(blockS * FRAME_HZ), runs = [];
  starts.forEach((s, r) => {
    const n = (starts[r + 1] ?? M) - s, y0 = stamps[s], xm = (n - 1) / 2;
    let ym = 0, sxy = 0, sxx = 0;
    for (let k = 0; k < n; k++) ym += stamps[s + k] - y0;
    ym /= n;
    for (let k = 0; k < n; k++) [sxy, sxx] = [sxy + (k - xm) * (stamps[s + k] - y0 - ym), sxx + (k - xm) ** 2];
    const b = n > 1 ? sxy / sxx : 1000 / FRAME_HZ, a = y0 + ym - b * xm, run = { s, n, a, b, xs: [s], cs: [0], short: n < SHORT };
    if (!run.short) {
      const res = new Float64Array(n).map((_, k) => stamps[s + k] - a - b * k), edges = [];
      for (let e = 0; e < n; e += B) edges.push(e);
      if (edges.length > 1 && n - edges.at(-1) < B / 2) edges.pop(); // a short last block joins the one before
      run.xs = edges.map((e, k) => s + (e + (edges[k + 1] ?? n) - 1) / 2);
      run.cs = edges.map((e, k) => median(res.subarray(e, edges[k + 1] ?? n)));
    }
    runs.push(run);
  });
  const interp = ({ xs, cs }, m) => {
    let lo = 0, hi = xs.length - 1;
    if (m <= xs[0]) return cs[0];
    if (m >= xs[hi]) return cs[hi];
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xs[mid] <= m) lo = mid; else hi = mid; }
    return cs[lo] + ((cs[hi] - cs[lo]) * (m - xs[lo])) / (xs[hi] - xs[lo]);
  };
  // PCG seconds of a fractional frame number
  const timeOf = (m) => {
    let r = runs.length - 1;
    while (r && m < runs[r].s) r--;
    const run = runs[r], host = run.a + run.b * (m - run.s) + interp(run, m);
    return (pcgLine.i0 + (host - pcgLine.t0 - pcgLine.a) / pcgLine.b) / 1000;
  };
  // The inverse, by bisection (timeOf only rises). Outside the recording it returns -1 or M, which is just "outside"
  const frameAt = (t) => {
    let lo = -1, hi = M;
    for (let k = 0; k < 36; k++) { const mid = (lo + hi) / 2; if (timeOf(mid) < t) lo = mid; else hi = mid; }
    return (lo + hi) / 2;
  };
  const long = runs.reduce((p, q) => (q.n > p.n ? q : p));
  return {
    timeOf, frameAt, runs, short: runs.filter((r) => r.short).length, frameHz: 1000 / long.b,
    maxC: Math.max(...runs.flatMap((r) => r.cs.map(Math.abs))), // largest move of the map from the straight line, ms
  };
}

// ---- Labeling trace, R peaks, T-end points and their flags. ECG frames are 500 Hz, so a position is a frame number and 1 frame = 2 ms.
// scipy.signal.butter(2, 0.5 / 250, 'high'): the labeling trace's one high-pass (run forward and backward: -6 dB at 0.5 Hz) in place of the viewer's two
export const HP05 = [[0.9955669720176472, -1.9911339440352944, 0.9955669720176472], [1, -1.9911142922016536, 0.9911535958689355]];
// scipy.signal.butter(2, [5 / 250, 15 / 250], 'band'), the R detector's QRS band
const BAND = [
  [0.003621681514928638, 0, -0.007243363029857276, 0, 0.003621681514928638],
  [1, -3.8000503652844575, 5.439339787793407, -3.4763426471814802, 0.8371816512560225],
];
const MS = 1000 / FRAME_HZ, sec = (s) => Math.round(s * FRAME_HZ); // sec(0.3): 300 ms as a whole number of frames
const MIN_RUN = 5 * FRAME_HZ; // shortcut: a run under 5 s gets no beats, upgrade if a recording with many short runs matters

// R peaks of a run of raw counts, as fractional frame numbers. x: the labeling trace of the same frames (the peaks are refined on it).
// Pan-Tompkins: band-pass, derivative, square, 150 ms mean, peaks 200 ms apart, adaptive signal and noise thresholds
export function detectR(raw, x) {
  const n = x.length, f = filtfilt(BAND, Float64Array.from(raw)), sq = new Float64Array(n + 1), half = 37;
  for (let i = 0; i < n; i++) {
    const g = (f[Math.min(i + 1, n - 1)] - f[Math.max(i - 1, 0)]) / (i && i < n - 1 ? 2 : 1);
    sq[i + 1] = sq[i] + g * g;
  }
  const mwi = Float64Array.from({ length: n }, (_, i) => (sq[Math.min(n, i + half + 1)] - sq[Math.max(0, i - half)]) / (2 * half + 1));
  const pk = [], keep = [];
  for (let i = 1; i < n - 1; i++) if (mwi[i - 1] < mwi[i] && mwi[i] > mwi[i + 1]) pk.push(i); // a flat-topped peak is skipped
  for (const j of pk.map((_, j) => j).sort((a, b) => mwi[pk[a]] - mwi[pk[b]]).reverse()) { // scipy find_peaks(distance=100): tallest first
    if (keep[j] === 0) continue;
    for (let k = j - 1; k >= 0 && pk[j] - pk[k] < sec(0.2); k--) keep[k] = 0;
    for (let k = j + 1; k < pk.length && pk[k] - pk[j] < sec(0.2); k++) keep[k] = 0;
  }
  const peaks = pk.filter((_, j) => keep[j] !== 0);
  if (!peaks.length) return [];
  const first = Float64Array.from(peaks.slice(0, 200), (p) => mwi[p]).sort(), pct = (q) => {
    const h = q * (first.length - 1), i = Math.floor(h);
    return first[i] + (first[Math.min(i + 1, first.length - 1)] - first[i]) * (h - i);
  };
  let spki = 0.5 * pct(0.9), npki = pct(0.1); // signal and noise levels, 0.125 smoothing
  const cand = [];
  for (const p of peaks) {
    if (mwi[p] > npki + 0.25 * (spki - npki)) {
      cand.push(p);
      spki = 0.125 * mwi[p] + 0.875 * spki;
    } else npki = 0.125 * mwi[p] + 0.875 * npki;
  }
  const pol = median(cand.map((c) => x[c])) > median(x) ? 1 : -1, R = [];
  for (const c of cand) { // the extreme of the trace within +-60 ms, then a parabola through three samples
    const lo = Math.max(0, c - sec(0.06)), hi = Math.min(n, c + sec(0.06) + 1);
    let k = lo;
    for (let i = lo + 1; i < hi; i++) if (pol * x[i] > pol * x[k]) k = i;
    if (k > 0 && k < n - 1) {
      const [y0, y1, y2] = [pol * x[k - 1], pol * x[k], pol * x[k + 1]], den = y0 - 2 * y1 + y2;
      R.push(k + (den ? (0.5 * (y0 - y2)) / den : 0));
    }
  }
  return R;
}

// Martinez wavelet: the scale 2^5 detail of the quadratic spline filters h = [1 3 3 1]/8 and g = [2 -2] (a trous), 62 taps
const KER = (() => {
  const conv = (a, b) => a.reduce((o, p, i) => (b.forEach((q, j) => (o[i + j] += p * q)), o), new Array(a.length + b.length - 1).fill(0));
  const up = (f, k) => f.flatMap((v, i) => (i ? [...Array(2 ** k - 1).fill(0), v] : [v]));
  let s = [1], ker;
  for (let k = 0; k < 5; k++) [ker, s] = [conv(s, up([2, -2], k)), conv(s, up([1, 3, 3, 1].map((v) => v / 8), k))];
  return ker;
})();

// T-end of every beat of a run: wavelet (the default), tangent (the second opinion), and the T amplitude. x: trace, R: positions from detectR.
// Per beat j, between R[j] and R[j + 1]: baseline = median of the TP segment 300 to 220 ms before the next R (no proposal if it would
// start under 450 ms after R), T peak = the largest deviation 120 ms after R to 0.5 of the RR. Positions are frames, NaN = none
export function tEnds(x, R) {
  const n = x.length, off = (KER.length - 1) >> 1, W = new Float64Array(n), dx = new Float64Array(n); // W: numpy convolve(x, KER, 'same')
  for (let i = 0; i < n; i++) {
    let w = 0;
    for (let k = Math.max(0, i + off - n + 1); k < KER.length && i + off - k >= 0; k++) w += KER[k] * x[i + off - k];
    W[i] = w;
    dx[i] = (x[Math.min(i + 1, n - 1)] - x[Math.max(i - 1, 0)]) / (i && i < n - 1 ? 2 : 1);
  }
  return R.map((rj, j) => {
    const none = { e: NaN, tan: NaN, amp: NaN };
    if (j + 1 >= R.length) return none;
    const r = Math.round(rj), rn = Math.round(R[j + 1]), rr = rn - r, a = rn - sec(0.3), b = rn - sec(0.22);
    if (a < r + sec(0.45)) return none;
    const base = median(x.subarray(a, b)), lo = r + sec(0.12), hi = Math.trunc(r + Math.min(sec(0.5), 0.55 * rr));
    if (hi <= lo) return none;
    let p = lo;
    for (let i = lo + 1; i < hi; i++) if (Math.abs(x[i] - base) > Math.abs(x[p] - base)) p = i;
    const s = Math.sign(x[p] - base), out = { e: NaN, tan: NaN, amp: Math.abs(x[p] - base) };
    // tangent: where the steepest fall after the peak, run on down to the baseline
    let k = p;
    for (let i = p + 1; i < Math.min(n - 2, p + sec(0.16)); i++) if (-s * dx[i] > -s * dx[k]) k = i;
    if (dx[k] * s < 0) out.tan = k + (base - x[k]) / dx[k];
    // wavelet: after the zero crossing at the T peak, the extreme of the opposite sign, then the first point back above 25 % of it
    const wl = r + sec(0.08), wh = Math.min(n - 2, r + Math.trunc(0.65 * rr));
    let z = wl;
    while (z < wh && !(s * W[z] > 0 && s * W[z + 1] <= 0 && Math.abs(z - p) < sec(0.08))) z++;
    if (z >= wh) return out;
    let m = z;
    for (let i = z + 1; i <= wh; i++) if (s * W[i] < s * W[m]) m = i;
    for (let i = m; i <= wh && s * W[m] < 0; i++) if (s * W[i] > 0.25 * s * W[m]) { out.e = i; break; }
    return out;
  });
}

// The proposed beats of a recording. raw: all raw counts, label: the 1000 Hz labeling trace (frame m is sample 2m, NaN where a run is
// too short to filter), starts: where each run starts (0 and the resyncs), map: buildMap's. Each beat: R and its T-end as a frame number
// (m, e, tan) and in PCG seconds (t, te), the T amplitude in counts, the RR to the next beat in ms, and doubt: why it should be looked at first
export function findBeats(raw, label, starts, map) {
  const beats = [];
  starts.forEach((s, r) => {
    const end = starts[r + 1] ?? raw.length;
    if (end - s < MIN_RUN || Number.isNaN(label[2 * s])) return;
    const x = Float32Array.from({ length: end - s }, (_, k) => label[2 * (s + k)]), R = detectR(raw.subarray(s, end), x), T = tEnds(x, R);
    R.forEach((m, j) => beats.push({
      m: s + m, t: map.timeOf(s + m), e: s + T[j].e, te: map.timeOf(s + T[j].e), tan: s + T[j].tan, amp: T[j].amp,
      rr: j + 1 < R.length ? (R[j + 1] - m) * MS : NaN, doubt: [],
    }));
  });
  const med = (f) => median(beats.map(f).filter(Number.isFinite)), dMed = med((b) => b.e - b.tan), rrMed = med((b) => b.rr), ampMed = med((b) => b.amp);
  for (const b of beats) {
    if (b.amp < 0.5 * ampMed) [b.e, b.te] = [NaN, NaN], b.doubt.push('small T'); // no T wave: the wavelet would only follow noise, so no proposal
    if (!Number.isFinite(b.e)) b.doubt.push('no T-end'); // also the last beat of a run, which has no next R
    if (Math.abs(b.e - b.tan - dMed) * MS > 20) b.doubt.push('methods disagree');
    if (b.rr < 0.8 * rrMed || b.rr > 1.2 * rrMed) b.doubt.push('RR outlier');
  }
  return beats;
}

// Where a hand-placed R belongs: the sample of x with the largest deviation from the baseline within 20 samples of c (the baseline is the mean
// of 600 samples around c, the QRS being narrow), refined by a parabola when it is a peak. A fractional index, or null when x has no data there
export function snapPeak(x, c) {
  let base = 0, n = 0, k = -1, best = -1;
  for (let i = Math.max(0, c - 300); i < Math.min(x.length, c + 300); i++) if (x[i] === x[i]) base += x[i], n++;
  base /= n;
  for (let i = Math.max(1, c - 20); i <= Math.min(x.length - 2, c + 20); i++) if (Math.abs(x[i] - base) > best) [k, best] = [i, Math.abs(x[i] - base)];
  if (k < 0) return null;
  const [y0, y1, y2] = [-1, 0, 1].map((d) => Math.sign(x[k] - base) * (x[k + d] - base)), den = y0 - 2 * y1 + y2;
  return k + (y1 >= y0 && y1 >= y2 && den ? (0.5 * (y0 - y2)) / den : 0); // no parabola unless k is a peak: on a slope it would fly off
}

// Independent check of the map: the PCG envelope (sqrt of the 20 ms mean of the squared samples) averaged on the R train per 60 s block,
// the block's lag read as the best correlation over the first 150 ms after R against the whole recording's, +-15 ms in 0.25 ms steps and a
// parabola. rMs: the R positions in PCG ms. Returns the lag of each block with at least 30 beats against the median of the blocks, a
// report and no adjustment: worst is the largest of them, and warn says it is over limitMs
export function alignCheck(pcg, rMs, limitMs = 2) {
  const n = pcg.length, mid = median(pcg.filter(Number.isFinite)), cum = new Float64Array(n + 1), LEN = 311; // the window is -80 to +230 ms
  for (let i = 0; i < n; i++) cum[i + 1] = cum[i] + ((pcg[i] - mid) ** 2 || 0); // a dropped frame (NaN) counts as nothing
  const env = Float64Array.from({ length: n }, (_, i) => Math.sqrt((cum[Math.min(n, i + 10)] - cum[Math.max(0, i - 10)]) / 20));
  const ensemble = (rs) => {
    const e = new Float64Array(LEN);
    let c = 0;
    for (const r of rs) {
      const i0 = Math.floor(r) - 80, f = r - Math.floor(r);
      if (i0 < 0 || i0 + LEN >= n) continue;
      for (let k = 0; k < LEN; k++) e[k] += env[i0 + k] * (1 - f) + env[i0 + k + 1] * f;
      c++;
    }
    return e.map((v) => v / c);
  };
  const at = (E, p) => { const x = p + 80, i = Math.floor(x); return E[i] + (E[i + 1] - E[i]) * (x - i); }; // p: ms after R
  const whole = ensemble(rMs), w = Array.from({ length: 151 }, (_, i) => i), tw = w.map((p) => at(whole, p)), blocks = new Map();
  for (const r of rMs) {
    const b = Math.floor(r / 60000);
    (blocks.get(b) ?? blocks.set(b, []).get(b)).push(r);
  }
  const lags = [];
  for (const [, rs] of [...blocks].sort((a, b) => a[0] - b[0])) {
    if (rs.length < 30) continue;
    const E = ensemble(rs), c = [];
    for (let L = -15; L <= 15; L += 0.25) {
      const y = w.map((p) => at(E, p + L)), my = y.reduce((s, v) => s + v) / y.length, mt = tw.reduce((s, v) => s + v) / tw.length;
      let sxy = 0, sxx = 0, syy = 0;
      y.forEach((v, i) => [sxy, sxx, syy] = [sxy + (v - my) * (tw[i] - mt), sxx + (v - my) ** 2, syy + (tw[i] - mt) ** 2]);
      c.push(sxy / Math.sqrt(sxx * syy));
    }
    const k = c.indexOf(Math.max(...c)), den = c[k - 1] - 2 * c[k] + c[k + 1];
    lags.push(-(-15 + 0.25 * k + (k > 0 && k < c.length - 1 ? (0.5 * (c[k - 1] - c[k + 1])) / den * 0.25 : 0)));
  }
  const m = median(lags), dev = lags.map((v) => v - m), worst = Math.max(...dev.map(Math.abs));
  return { lags: dev, worst, warn: lags.length >= 2 && worst > limitMs, limitMs };
}
