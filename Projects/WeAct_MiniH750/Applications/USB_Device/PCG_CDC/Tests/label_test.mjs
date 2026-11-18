// Self-check of Web/align.js (the ECG-to-PCG map, R and T-end proposals, the S1 check), on synthetic data only: a recording of a person never goes into a test or a fixture.
// From the project folder: node --test Tests/label_test.mjs
import test from 'node:test';
import assert from 'node:assert/strict';
import { buildMap, detectR, tEnds, findBeats, alignCheck, snapPeak, HP05 } from '../Web/align.js';
import { zeroPhase, NF } from '../Web/viewer_data.js';
import { createLabels, fromBeats, groupBeats, counter, readLabels, writeLabels, describeMap, DETECT, STATES } from '../Web/labels.js';

const PCG = { i0: 0, t0: 0, a: 0, b: 1 }; // host ms = PCG ms, so the map's seconds are the host's
let seed = 12345;
const noise = () => ((seed = (seed * 1103515245 + 12345) % 2 ** 31) / 2 ** 31 - 0.5) * 2; // uniform in -1..1 ms, deterministic
const median = (x) => { const s = [...x].sort((p, q) => p - q), h = s.length >> 1; return s.length % 2 ? s[h] : (s[h - 1] + s[h]) / 2; };

// 600 s at a clock 1000 ppm fast, with a 3 ms sine wander (one and a half periods) and 1 ms of noise on every stamp
const N = 300000, truth = Float64Array.from({ length: N }, (_, m) => 1000 + (m * 2) / 1.001 + 3 * Math.sin((2 * Math.PI * m) / 200000));
const stamps = truth.map((v) => v + noise());

test('30 s blocks follow a wander that one line cannot', () => {
  const map = buildMap(stamps, [], PCG, 30), run = map.runs[0];
  let worstMap = 0, worstLine = 0;
  for (let e = 0; e + 15000 <= N; e += 15000) {
    const at = (f) => Array.from({ length: 15000 }, (_, k) => f(e + k));
    worstMap = Math.max(worstMap, Math.abs(median(at((m) => map.timeOf(m) * 1000 - truth[m]))));
    worstLine = Math.max(worstLine, Math.abs(median(at((m) => run.a + run.b * m - truth[m]))));
  }
  assert.ok(worstMap < 0.5, `map ${worstMap}`);
  assert.ok(worstLine > 2, `line ${worstLine}`);
  assert.ok(Math.abs(map.frameHz / (1000 * 1.001 / 2) - 1) < 1e-5, `rate ${map.frameHz}`);
  assert.ok(map.maxC > 1 && map.maxC < 4, `max |c| ${map.maxC}`); // how far the map moves from the line, ms
  assert.equal(map.short, 0);
});

test('time and frame are inverse, also between frames', () => {
  const map = buildMap(stamps, [], PCG, 30);
  for (const m of [0, 0.25, 1, 7777.5, 150000, 299998.75, 299999]) assert.ok(Math.abs(map.frameAt(map.timeOf(m)) - m) < 1e-3, `frame ${m}`);
  assert.ok(map.frameAt(-5) < 0 && map.frameAt(1e6) > N - 1); // outside the recording
});

test('a resync starts a run with its own line, and a run under 10 s keeps its line and is flagged', () => {
  // frames 0..99999 as above, then 500 ms of lost bytes and 100000 more, then a 5 s run (2500 frames) after another 100 ms
  const s2 = Float64Array.from({ length: 202500 }, (_, m) => (m < 100000 ? stamps[m] : m < 200000 ? stamps[m] + 500 : stamps[m] + 600));
  const map = buildMap(s2, [100000, 200000], PCG, 30);
  assert.equal(map.runs.length, 3);
  assert.equal(map.short, 1);
  assert.equal(map.runs[2].short, true);
  assert.deepEqual(map.runs[2].cs, [0]); // a short run has the line only
  for (const m of [99999, 100000, 150000, 199999, 200000, 202499]) assert.ok(Math.abs(map.timeOf(m) * 1000 - s2[m]) < 2, `frame ${m}: ${map.timeOf(m) * 1000} against ${s2[m]}`);
  assert.ok(map.timeOf(100000) - map.timeOf(99999) > 0.5); // the lost time is in the map
  assert.ok(map.timeOf(202499) > map.timeOf(200000));
});

test('a stalled stamp moves nothing: the block median ignores it', () => {
  const s2 = stamps.slice();
  for (let m = 70000; m < 70033; m++) s2[m] += 66; // a 66 ms stall of the host
  const a = buildMap(stamps, [], PCG, 30), b = buildMap(s2, [], PCG, 30);
  for (const m of [10000, 69000, 71000, 200000]) assert.ok(Math.abs(a.timeOf(m) - b.timeOf(m)) * 1000 < 0.05, `frame ${m}`);
});

// ---- R peaks, T-end points and the S1 check. 120 s of ECG at 500 Hz: a Gaussian QRS (sigma 8 ms, 60 counts) at known fractional times, a half-sine T wave
// 180 to 340 ms after each R (20 counts), 0.3 Hz wander, 50 Hz hum, 1 count of noise; beat 20 has no T wave; the RR wanders between 0.8 and 1.0 s
const R0 = [];
for (let k = 0, t = 0.7; t < 119; k++) R0.push(t), t += 0.9 + 0.1 * Math.sin(k * 0.7) + 0.0003 * k; // the 0.3 ms steps make the times fall between samples
const ecg = (wander) => Int16Array.from({ length: 60000 }, (_, i) => {
  const t = i / 500;
  let v = -1300 + wander * Math.sin(2 * Math.PI * 0.3 * t) + 5 * Math.sin(2 * Math.PI * 50 * t) + noise();
  R0.forEach((r, k) => {
    v += 60 * Math.exp(-0.5 * ((t - r) / 0.008) ** 2);
    const u = t - r - 0.18;
    if (u > 0 && u < 0.16 && k !== 20) v += 20 * Math.sin((Math.PI * u) / 0.16);
  });
  return Math.round(v);
});
const labelTrace = (raw) => zeroPhase(raw, [HP05, NF]), at500 = (tr) => Float32Array.from({ length: (tr.length + 1) / 2 }, (_, m) => tr[2 * m]);

test('every R is found, to the millisecond, with hum and wander', () => {
  const raw = ecg(40), R = detectR(raw, at500(labelTrace(raw)));
  assert.equal(R.length, R0.length);
  R.forEach((m, k) => assert.ok(Math.abs(m / 500 - R0[k]) < 0.002, `beat ${k}: ${m / 500} against ${R0[k]}`));
});

test('T-end: within 30 ms of the model, steady under a wander twice the T wave, none for a missing T wave or the last beat', () => {
  const raw = ecg(0), wavy = ecg(40), R = detectR(raw, at500(labelTrace(raw)));
  const T = tEnds(at500(labelTrace(raw)), R), Tw = tEnds(at500(labelTrace(wavy)), R);
  T.forEach((o, k) => {
    if (k === 20 || k === R.length - 1) return; // no T wave, no next R: below
    assert.ok(Math.abs(o.e / 500 - (R0[k] + 0.34)) < 0.03, `beat ${k}: T-end ${o.e / 500} against ${R0[k] + 0.34}`);
    assert.ok(Math.abs(Tw[k].e - o.e) * 2 < 10, `beat ${k} moved ${Math.abs(Tw[k].e - o.e) * 2} ms by the wander`);
  });
  assert.ok(Number.isNaN(T.at(-1).e) && Number.isNaN(T.at(-1).tan));
  const stamps = Float64Array.from({ length: 60000 }, (_, m) => 2 * m), beats = findBeats(raw, labelTrace(raw), [0], buildMap(stamps, [], PCG, 30));
  assert.ok(beats.at(-1).doubt.includes('no T-end'));
  assert.ok(Number.isNaN(beats[20].e) && beats[20].doubt.includes('no T-end') && beats[20].doubt.includes('small T'), `beat 20 has no T wave: ${beats[20].e} ${beats[20].doubt}`);
  assert.ok(Math.abs(beats[5].t - R0[5]) < 0.002); // the time is the map's: here the stamps are the frame numbers times 2 ms
});

// 600 s of PCG at 1 kHz: a 60 Hz damped burst (S1) 20 ms after each R, 5 counts of noise
const rP = [];
for (let k = 0, t = 0.5; t < 595; k++) rP.push(1000 * t), t += 0.9 + 0.1 * Math.sin(k * 0.7);
const pcg = Float32Array.from({ length: 600000 }, (_, i) => 5 * noise());
rP.forEach((r) => { for (let j = 20; j < 200; j++) pcg[Math.round(r) + j] += 100 * Math.sin((2 * Math.PI * 60 * (j - 20)) / 1000) * Math.exp(-(j - 20) / 30); });

test('the S1 check: a flat lag is quiet, a ramp of the R train warns', () => {
  const flat = alignCheck(pcg, rP);
  assert.equal(flat.lags.length, 10);
  assert.ok(!flat.warn && flat.worst < 0.5, `flat: worst ${flat.worst}`);
  const ramp = alignCheck(pcg, rP.map((r) => r + (6 * r) / 600000)); // the R train drifts 6 ms over the recording (a block is within 2.7 ms of the median)
  assert.ok(ramp.warn && ramp.worst > 2, `ramp: worst ${ramp.worst}`);
  assert.ok(!alignCheck(pcg, rP.slice(0, 40)).warn); // under 2 blocks of 30 beats: no verdict
});

// ---- Labels: states, operations, undo, the file. Synthetic: 6 beats, R every 0.9 s and its T-end 0.34 s later; a placeholder folder name
const found = (n = 6) => Array.from({ length: n }, (_, k) => ({ t: 1 + 0.9 * k, te: k === 3 ? NaN : 1.34 + 0.9 * k, doubt: k === 2 ? ['RR outlier'] : k === 3 ? ['no T-end'] : [] }));
const META = { folder: 'recording-folder', segments: 4, pcg_samples: 6000, ecg_frames: 3000, pcg_drops: 0 };
const at = (L, kind, t) => L.list.findIndex((l) => l.kind === kind && Math.abs(l.t - t) < 1e-6);
const snap = (L) => JSON.stringify(L.list);

test('move: proposed and confirmed become edited, added stays added, a deleted label cannot move', () => {
  const L = createLabels(fromBeats(found()));
  L.confirm(at(L, 'R', 1.9));
  const i = L.move(at(L, 'R', 1), 1.005), j = L.move(at(L, 'R', 1.9), 1.95);
  assert.deepEqual([L.list[i].t, L.list[i].state, L.list[j].t, L.list[j].state], [1.005, 'edited', 1.95, 'edited']);
  const k = L.add(20, 20.34), a = L.move(k, 20.01);
  assert.deepEqual([L.list[a].state, L.list[a].t], ['added', 20.01]);
  const d = at(L, 'T', 1.34);
  L.delete(d);
  assert.equal(L.move(at(L, 'T', 1.34), 1.4), -1);
  assert.equal(L.move(a, 20.01), -1); // not moved: no change, no history step
});

test('delete keeps the row as deleted, takes the T-end of the beat along and leaves the counter', () => {
  const L = createLabels(fromBeats(found())), n = L.list.length;
  L.delete(at(L, 'R', 1.9));
  assert.equal(L.list.length, n);
  assert.deepEqual(L.list.filter((l) => l.state === 'deleted').map((l) => [l.kind, l.t]), [['R', 1.9], ['T', 2.24]]);
  assert.deepEqual(counter(L.list), { beats: 5, reviewed: 0, flagged: 2 }); // beats 2 and 3 are still there, the deleted one was beat 1
  assert.equal(L.delete(at(L, 'R', 1.9)), -1); // already deleted: nothing
});

test('add makes an added R with a proposed T-end; setT moves a T or adds one where the beat has none', () => {
  const L = createLabels(fromBeats(found()));
  const i = L.add(10.5, 10.84);
  assert.deepEqual([L.list[i].kind, L.list[i].state, L.list[i + 1].kind, L.list[i + 1].state], ['R', 'added', 'T', 'proposed']);
  const t3 = L.setT(4.2); // beat 3 (R at 3.7) has no T-end
  assert.deepEqual([L.list[t3].kind, L.list[t3].state, L.list[t3].t], ['T', 'added', 4.2]);
  const t0 = L.setT(1.4);
  assert.deepEqual([L.list[t0].state, L.list[t0].t], ['edited', 1.4]);
  assert.equal(L.setT(0.5), -1); // no beat to its left
  assert.equal(groupBeats(L.list).filter((b) => b.t).length, 7);
});

test('delete alone takes only that R and leaves its T-end; one argument still takes both; undo restores the list', () => {
  const L = createLabels(fromBeats(found())), before = snap(L);
  const k = L.delete(at(L, 'R', 1.9), true);
  assert.deepEqual(L.list.filter((l) => l.state === 'deleted').map((l) => [l.kind, l.t]), [['R', 1.9]]);
  assert.equal(L.list[at(L, 'T', 2.24)].state, 'proposed');
  assert.ok(k >= 0 && L.undo() && snap(L) === before);
  L.delete(at(L, 'R', 1.9));
  assert.equal(L.list.filter((l) => l.state === 'deleted').length, 2);
});

test('add proposes no T-end when a live T already lies before the next live R, and both when none does', () => {
  const L = createLabels(fromBeats(found()));
  L.delete(at(L, 'R', 1.9), true); // its T at 2.24 stays
  const i = L.add(1.9, 2.241);
  assert.equal(L.list[i].kind, 'R');
  assert.equal(L.list.filter((l) => l.kind === 'T' && l.state !== 'deleted').length, 5); // no second T (beat 3 never had one)
  assert.equal(L.list.filter((l) => l.state === 'added').length, 1);
  const j = L.add(10.5, 10.84);
  assert.deepEqual([L.list[j].kind, L.list[j + 1].kind, L.list[j + 1].state], ['R', 'T', 'proposed']);
});

test('undo and redo restore the exact list after any 10 operations, the history is 100 deep', () => {
  const L = createLabels(fromBeats(found())), states = [snap(L)];
  const ops = [() => L.confirm(0, 1), () => L.move(0, 1.01), () => L.delete(at(L, 'R', 1.9)), () => L.add(7, 7.3), () => L.setT(4.2),
    () => L.move(at(L, 'T', 3.14), 3.2), () => L.proposeAgain(found(8)), () => L.delete(at(L, 'T', 1.34)), () => L.add(30), () => L.confirm(at(L, 'R', 2.8))];
  for (const op of ops) { op(); states.push(snap(L)); }
  assert.equal(new Set(states).size, 11); // every operation changed something
  for (let k = 10; k > 0; k--) { assert.ok(L.undo()); assert.equal(snap(L), states[k - 1], `undo to ${k - 1}`); }
  assert.ok(!L.undo() && !L.canUndo);
  for (let k = 1; k <= 10; k++) { assert.ok(L.redo()); assert.equal(snap(L), states[k], `redo to ${k}`); }
  assert.ok(!L.redo());
  L.undo();
  L.move(0, 1.02); // a new change drops the redo steps
  assert.ok(!L.canRedo);
  const M = createLabels(fromBeats(found()));
  for (let k = 0; k < 105; k++) M.move(0, 1 + 0.001 * (k + 1));
  let u = 0;
  while (M.undo()) u++;
  assert.equal(u, 100);
});

test('save then load returns the same labels and states, and the same file again', () => {
  const L = createLabels(fromBeats(found()));
  L.confirm(at(L, 'R', 1)); L.move(at(L, 'R', 1.9), 1.93); L.add(8, 8.3); L.delete(at(L, 'T', 2.24)); L.delete(at(L, 'R', 4.6));
  assert.deepEqual([...new Set(L.list.map((l) => l.state))].sort(), [...STATES].sort()); // all five states are in the file
  const map = describeMap(buildMap(stamps, [], PCG, 30), 30, PCG), text = writeLabels(L.list, { ...META, map, detect: DETECT });
  const back = readLabels(text, { segments: 4, pcg_samples: 6000, ecg_frames: 3000 });
  assert.deepEqual(back.list, L.list);
  assert.deepEqual(back.meta, { ...META, map, detect: DETECT });
  assert.equal(writeLabels(back.list, back.meta), text);
  const head = text.split('\n');
  assert.equal(head[0], '# pcg-ecg-labels 1');
  assert.match(head[1], /^# recording: folder=recording-folder segments=4 pcg_samples=6000 ecg_frames=3000 pcg_drops=0$/);
  assert.match(head[2], /^# time: .*PCG sample index \/ 1000.*loader convention/);
  assert.equal(head[5], 'beat,kind,t_pcg_s,state,doubt');
  assert.equal(head.length - 7, L.list.length); // one row per label, deleted ones included, and a closing newline
  assert.match(head[6], /^0,R,1\.000000,confirmed,$/);
  assert.ok(back.meta.map.runs[0].block_median_ms.length > 10 && back.meta.map.block_s === 30);
  assert.equal(readLabels(text).list.length, L.list.length); // no identity given: no check
});

test('a label file of another recording, or a damaged one, is refused with a message', () => {
  const text = writeLabels(fromBeats(found()), META), want = { segments: 4, pcg_samples: 6000, ecg_frames: 3000 };
  assert.throws(() => readLabels(text, { ...want, ecg_frames: 3001 }), /another recording \(ecg_frames: file 3000, this recording 3001\)/);
  assert.throws(() => readLabels(text, { ...want, segments: 5 }), /another recording/);
  assert.throws(() => readLabels('hello\n'), /not a label file/);
  assert.throws(() => readLabels(text.replace('1.000000,proposed', '1.000000,maybe')), /line 7: not a label row/);
  assert.throws(() => readLabels(text.replace('0,R,1.000000', '0,X,1.000000')), /not a label row/);
  assert.throws(() => readLabels(text.replace('0,R,1.000000', '0,R,abc')), /not a label row/);
  assert.throws(() => readLabels(text.replace('# map: null', '# map: {')), /map header is not valid JSON/);
  assert.throws(() => readLabels(text.replace(/# recording:.*\n/, '')), /no recording header/);
  assert.throws(() => readLabels(text.split('\n').slice(0, 5).join('\n')), /no recording header or no label table|no label table/);
});

test('propose again keeps every label that is not proposed and every label near one, and replaces the rest', () => {
  const L = createLabels(fromBeats(found()));
  L.move(at(L, 'R', 1), 1.02); L.confirm(at(L, 'T', 1.34)); L.delete(at(L, 'R', 2.8)); L.add(8, 8.3); // a hand-made R where the detector found none
  const kept = L.list.filter((l) => l.state !== 'proposed');
  const again = found(7).map((b, k) => { const d = k === 5 ? 0.3 : 0.002; return { ...b, t: b.t + d, te: b.te + d }; }); // shifted a little, a 7th beat, beat 5 far
  const dry = L.proposeAgain(again, true), listBefore = snap(L);
  assert.equal(snap(L), listBefore); // a dry run changes nothing
  const out = L.proposeAgain(again);
  assert.deepEqual(out, dry);
  assert.deepEqual(out, { kept: 6, replaced: 7, proposed: 9 }); // hand count: 5 labels not proposed and the hand-made R's proposed T stay
  for (const k of kept) assert.ok(L.list.some((l) => l.kind === k.kind && l.t === k.t && l.state === k.state), `kept ${k.kind} ${k.t}`);
  const times = (kind) => L.list.filter((l) => l.kind === kind && l.state === 'proposed').map((l) => l.t);
  // R at 1.002 and 2.802 are within 150 ms of the edited and the deleted R, the T at 1.342 and 3.142 of the confirmed and the deleted T
  assert.deepEqual(times('R'), [1.902, 3.702, 4.602, 5.8, 6.402]);
  assert.deepEqual(times('T'), [2.242, 4.942, 6.14, 6.742, 8.3]);
  assert.equal(out.kept + out.replaced, 13); // the list had 13 labels
  assert.ok(L.undo() && snap(L) === listBefore); // one undo step
});

test('the counter equals a hand count', () => {
  const L = createLabels(fromBeats(found())); // 6 beats, 5 T-ends, doubt on beats 2 and 3
  assert.deepEqual(counter(L.list), { beats: 6, reviewed: 0, flagged: 2 });
  L.confirm(at(L, 'R', 1), at(L, 'T', 1.34)); // beat 0 reviewed
  L.confirm(at(L, 'R', 1.9)); // beat 1: only the R, so not reviewed
  L.confirm(at(L, 'R', 3.7)); // beat 3 has no T-end: its R is enough
  L.move(at(L, 'R', 2.8), 2.81); L.move(at(L, 'T', 3.14), 3.15); // beat 2 edited, flagged but now reviewed
  L.add(9, 9.3); L.delete(at(L, 'R', 5.5)); // a new beat (T proposed, not reviewed), one beat gone
  assert.deepEqual(counter(L.list), { beats: 6, reviewed: 3, flagged: 0 });
});

test('a hand-placed R snaps to the peak within 20 samples, to a fraction of one, and never flies off on a slope', () => {
  const x = Float32Array.from({ length: 2000 }, (_, i) => 5 + 50 * Math.exp(-(((i - 1000.3) / 4) ** 2))); // a spike at 1000.3 on a baseline of 5
  assert.ok(Math.abs(snapPeak(x, 1012) - 1000.3) < 0.15, `${snapPeak(x, 1012)}`);
  assert.ok(Math.abs(snapPeak(x.map((v) => -v), 988) - 1000.3) < 0.15); // negative polarity
  const slope = Float32Array.from({ length: 2000 }, (_, i) => Math.exp((i - 1000) / 50)); // no peak: the largest deviation is at the window's edge (980), where a parabola would fly off
  assert.ok(Math.abs(snapPeak(slope, 1000) - 980) < 1, `${snapPeak(slope, 1000)}`);
  assert.equal(snapPeak(new Float32Array(2000).fill(NaN), 1000), null);
  assert.equal(snapPeak(x, -5000), null); // outside the trace
});
