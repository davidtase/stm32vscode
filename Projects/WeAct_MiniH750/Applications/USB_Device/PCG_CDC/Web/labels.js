// The R and T-end labels of a recording, their review states, undo and the label file (a CSV with # header lines), with no DOM, so node can test it.
// A label is { kind: 'R' | 'T', t: PCG seconds, state, doubt }. The list is kept sorted by time. A beat is a live R and the first live T after it.
// Nothing here knows the recording: proposals come from align.js findBeats, and the page passes in the traces' answers (a snapped R, a T-end).
export const STATES = ['proposed', 'confirmed', 'edited', 'added', 'deleted'];
const REVIEWED = new Set(['confirmed', 'edited', 'added']), COLUMNS = 'beat,kind,t_pcg_s,state,doubt', HISTORY = 100, GAP_S = 0.15;
const TIME_NOTE = 'seconds = PCG sample index / 1000, index 0 = first live PCG frame, loader convention (frame\'s first sample stamped at the arrival of its last byte)';
const us = (t) => Math.round(t * 1e6) / 1e6; // times live at the file's six decimals, so a saved and reloaded label is equal to the bit
const copy = (l) => ({ ...l });
const byTime = (a, b) => a.t - b.t || (a.kind === b.kind ? 0 : a.kind === 'R' ? -1 : 1);

// What the detectors did, written in the header so a label set says how it was made. Keep in step with align.js (findBeats, detectR, tEnds)
export const DETECT = {
  R: 'Pan-Tompkins on raw counts: 5-15 Hz zero-phase band-pass, derivative, square, 150 ms mean, peaks 200 ms apart, adaptive thresholds (noise + 0.25 (signal - noise)), refined on the labeling trace within 60 ms and a parabola',
  T_end: 'Martinez wavelet, quadratic spline, scale 2^5 at 500 Hz, end at 25 % of the last extreme; second opinion: tangent at the steepest down slope; baseline: median of 300 to 220 ms before the next R',
  trace: 'labeling trace: 2nd-order Butterworth high-pass at 0.5 Hz forward and backward, 50 Hz notch, 40 Hz FIR',
  doubt: 'T-end methods differ by over 20 ms from their recording median; RR under 0.8 or over 1.2 times the median; T amplitude under half the median; no T-end',
};

// The map's numbers for the header: block length, the PCG line and, per run, its line and block medians. map: buildMap's, pcgLine: loadRecording's pcgFit
export const describeMap = (map, blockS, pcgLine) => ({
  block_s: blockS, pcg_line: pcgLine && { i0: pcgLine.i0, t0: pcgLine.t0, a: pcgLine.a, b: pcgLine.b },
  runs: map.runs.map((r) => ({ start_frame: r.s, frames: r.n, line_a_ms: r.a, line_b_ms: r.b, short: r.short,
    block_centre_frames: r.xs.map((x) => Math.round(x * 10) / 10), block_median_ms: r.cs.map((c) => Math.round(c * 1e4) / 1e4) })),
});

// findBeats' output as labels, every one 'proposed'. A beat with no T-end gets no T row. The beat's doubt goes on its R
export const fromBeats = (found) => found.flatMap((b) => [
  { kind: 'R', t: us(b.t), state: 'proposed', doubt: b.doubt.join(';') },
  ...(Number.isFinite(b.te) ? [{ kind: 'T', t: us(b.te), state: 'proposed', doubt: '' }] : []),
]).sort(byTime);

// Beats of a list: each live R with the first live T before the next live R (deleted labels count nowhere, a T with no R before it is left out)
export function groupBeats(list) {
  const beats = [];
  for (const l of list) {
    if (l.state === 'deleted') continue;
    if (l.kind === 'R') beats.push({ r: l, t: null });
    else if (beats.length && !beats.at(-1).t) beats.at(-1).t = l;
  }
  return beats;
}

// A beat is reviewed when its R and its T (if it has one) are confirmed, edited or added
export const reviewed = (b) => REVIEWED.has(b.r.state) && (!b.t || REVIEWED.has(b.t.state));

// The counter: beats, reviewed beats and flagged beats not yet reviewed
export function counter(list) {
  const beats = groupBeats(list);
  return { beats: beats.length, reviewed: beats.filter(reviewed).length, flagged: beats.filter((b) => !reviewed(b) && b.r.doubt).length };
}

// A label list with its operations and history (100 snapshots of the whole list). Operations take an index into `list`, return the label's
// new index (the list is re-sorted) or -1 when nothing changed, and every change is one undo step.
export function createLabels(initial = []) {
  let cur = initial.map(copy).sort(byTime), past = [], future = [];
  const change = (fn) => {
    const next = cur.map(copy), hit = fn(next);
    if (!hit) return -1;
    past.push(cur);
    if (past.length > HISTORY) past.shift();
    cur = next.sort(byTime);
    future = [];
    return cur.indexOf(hit);
  };
  const live = (next, i) => next[i] && next[i].state !== 'deleted' && next[i];
  const step = (from, to) => { if (!from.length) return false; to.push(cur); cur = from.pop(); return true; };
  return {
    get list() { return cur; }, // read only: change it through the operations
    // proposed to confirmed; a confirmed, edited or added label is already reviewed. Several indices are one step (Enter confirms R and T)
    confirm: (...ix) => change((next) => ix.map((i) => live(next, i)).filter((l) => l?.state === 'proposed').map((l) => (l.state = 'confirmed', l)).at(0)),
    move: (i, t) => change((next) => {
      const l = live(next, i);
      if (!l || us(t) === l.t) return null;
      l.t = us(t);
      if (l.state === 'proposed' || l.state === 'confirmed') l.state = 'edited'; // an added label stays added
      return l;
    }),
    // a deleted label stays in the list, so a re-run does not bring it back. Deleting an R deletes the T of its beat too, unless alone
    delete: (i, alone = false) => change((next) => {
      const l = live(next, i), beat = l && groupBeats(next).find((b) => b.r === l);
      if (!l) return null;
      l.state = 'deleted';
      if (beat?.t && !alone) beat.t.state = 'deleted';
      return l;
    }),
    // a hand-made R (the page snaps it to a peak first) and its T-end proposal, te: PCG seconds or NaN for none; no proposal when a live T
    // already lies between it and the next live R (the T of a beat whose R was deleted alone)
    add: (t, te = NaN) => change((next) => {
      const r = { kind: 'R', t: us(t), state: 'added', doubt: '' }, nr = next.find((l) => l.kind === 'R' && l.state !== 'deleted' && l.t > r.t)?.t ?? Infinity;
      const own = next.some((l) => l.kind === 'T' && l.state !== 'deleted' && l.t > r.t && l.t < nr);
      next.push(r);
      if (Number.isFinite(te) && !own) next.push({ kind: 'T', t: us(te), state: 'proposed', doubt: '' });
      return r;
    }),
    // the T-end of the beat to the left of t: moves its T, or adds one (state added) when it has none or it was deleted. -1 with no beat there
    setT: (t) => change((next) => {
      const beat = groupBeats(next).filter((b) => b.r.t < t).at(-1);
      if (!beat) return null;
      if (!beat.t) return next[next.push({ kind: 'T', t: us(t), state: 'added', doubt: '' }) - 1];
      beat.t.t = us(t);
      if (beat.t.state === 'proposed' || beat.t.state === 'confirmed') beat.t.state = 'edited';
      return beat.t;
    }),
    // Replaces the proposed labels with `found` (findBeats' output) and nothing else. No proposal lands within 150 ms of a label that stays,
    // a deleted one included. A proposed T of a beat whose R is not proposed stays too, or a hand-added R would lose its T-end. dry: only the counts
    proposeAgain(found, dry = false) {
      const own = new Set(groupBeats(cur).filter((b) => b.t && b.r.state !== 'proposed').map((b) => b.t));
      const keep = cur.filter((l) => l.state !== 'proposed' || own.has(l)), fresh = [];
      const free = (kind, t) => !keep.some((l) => l.kind === kind && Math.abs(l.t - t) < GAP_S) && !fresh.some((l) => l.kind === kind && Math.abs(l.t - t) < GAP_S);
      for (const l of fromBeats(found)) if (free(l.kind, l.t)) fresh.push(l);
      const out = { kept: keep.length, replaced: cur.length - keep.length, proposed: fresh.length };
      if (!dry && (out.replaced || out.proposed)) change((next) => (next.splice(0, next.length, ...keep.map(copy), ...fresh), next[0] ?? {}));
      return out;
    },
    undo: () => step(past, future),
    redo: () => step(future, past),
    get canUndo() { return past.length > 0; },
    get canRedo() { return future.length > 0; },
  };
}

// meta: { folder, segments, pcg_samples, ecg_frames, pcg_drops, map, detect } (map and detect are JSON objects). One header line each, then the rows
export function writeLabels(list, meta) {
  const one = (s) => String(s).replace(/[\r\n]+/g, ' ');
  let beat = -1;
  const rows = [...list].sort(byTime).map((l) => `${l.kind === 'R' ? ++beat : Math.max(beat, 0)},${l.kind},${l.t.toFixed(6)},${l.state},${l.doubt}`);
  return [
    '# pcg-ecg-labels 1',
    `# recording: folder=${one(meta.folder)} segments=${meta.segments} pcg_samples=${meta.pcg_samples} ecg_frames=${meta.ecg_frames} pcg_drops=${meta.pcg_drops}`,
    `# time: ${TIME_NOTE}`,
    `# map: ${JSON.stringify(meta.map ?? null)}`,
    `# detect: ${JSON.stringify(meta.detect ?? DETECT)}`,
    COLUMNS, ...rows,
  ].join('\n') + '\n';
}

// Reads a label file. Throws an Error with a plain message for a file that is not one, is damaged, or belongs to another recording
// (want: { segments, pcg_samples, ecg_frames } of the open one). Returns { list, meta } and meta goes straight back into writeLabels
export function readLabels(text, want) {
  const lines = text.split(/\r?\n/);
  if (lines[0] !== '# pcg-ecg-labels 1') throw new Error('This is not a label file (it should start with "# pcg-ecg-labels 1").');
  const meta = {}, list = [];
  let table = false;
  const json = (s, name, at) => { try { return JSON.parse(s); } catch { throw new Error(`Label file line ${at}: the ${name} header is not valid JSON.`); } };
  lines.forEach((line, k) => {
    const at = k + 1;
    if (!line.trim() || k === 0) return;
    if (!table) {
      if (line === COLUMNS) table = true;
      else if (line.startsWith('# recording: ')) {
        const m = /^# recording: folder=(.*) segments=(\d+) pcg_samples=(\d+) ecg_frames=(\d+) pcg_drops=(\d+)$/.exec(line);
        if (!m) throw new Error(`Label file line ${at}: the recording header is damaged.`);
        Object.assign(meta, { folder: m[1], segments: +m[2], pcg_samples: +m[3], ecg_frames: +m[4], pcg_drops: +m[5] });
      } else if (line.startsWith('# map: ')) meta.map = json(line.slice(7), 'map', at);
      else if (line.startsWith('# detect: ')) meta.detect = json(line.slice(10), 'detect', at);
      else if (!line.startsWith('#')) throw new Error(`Label file line ${at}: expected the column line "${COLUMNS}".`);
      return;
    }
    const [, kind, t, state, ...rest] = line.split(',');
    if ((kind !== 'R' && kind !== 'T') || !STATES.includes(state) || !/^-?\d+(\.\d+)?$/.test(t)) throw new Error(`Label file line ${at}: not a label row.`);
    list.push({ kind, t: us(+t), state, doubt: rest.join(',') });
  });
  if (!table || meta.segments === undefined) throw new Error('The label file has no recording header or no label table.');
  if (want) for (const k of ['segments', 'pcg_samples', 'ecg_frames']) {
    if (meta[k] !== want[k]) throw new Error(`This label file is for another recording (${k}: file ${meta[k]}, this recording ${want[k]}).`);
  }
  return { list: list.sort(byTime), meta };
}
