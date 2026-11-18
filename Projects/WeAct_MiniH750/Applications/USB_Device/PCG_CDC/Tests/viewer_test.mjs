// Self-check of Web/viewer_data.js, on synthetic data only: a recording of a person never goes into a test or a fixture.
// From the project folder: node --test Tests/viewer_test.mjs
import test from 'node:test';
import assert from 'node:assert/strict';
import { record, createLog } from '../Web/recorder.js';
import { loadRecording, segmentList, columns } from '../Web/viewer_data.js';

const FRAMES = 20000; // 200 s of PCG, 480 KB: the log rolls into 3 segments with the ECG records
const SAMPLE = (k, j) => (k * 10 + j) % 20000;

function pcgFrame(seq, k) {
  const f = new Uint8Array(24), dv = new DataView(f.buffer);
  dv.setUint16(0, 0x8000, true);
  dv.setUint16(2, seq, true);
  for (let j = 0; j < 10; j++) dv.setInt16(4 + 2 * j, SAMPLE(k, j), true);
  return f;
}
const ecgValue = (i) => Math.fround(i / 8 - 500);

// A stale frame, then FRAMES live ones, cut into chunks of odd sizes. The PCG stamps run 1000 ppm fast, the ECG ones 500 ppm slow
async function synthetic() {
  const files = [], log = createLog(async (_, blob) => files.push(await blob.arrayBuffer()));
  log.event(500, { e: 'start' });
  log.event(501, { e: 'keep-alive' });
  const pcg = new Uint8Array(24 * (FRAMES + 1));
  pcg.set(pcgFrame(5, 0), 0);
  for (let k = 0; k < FRAMES; k++) pcg.set(pcgFrame(100 + k, k), 24 * (k + 1));
  const sizes = [1, 5, 17, 37, 23];
  let p = 0, c = 0, e = 0;
  const pcgMs = 1 / 1.001, ecgMs = 1 / 0.9995; // host ms per sample
  while (p < pcg.length) {
    const end = Math.min(pcg.length, p + sizes[c++ % sizes.length]);
    log.add(1000 + (end / 24) * 10 * pcgMs, [1, pcg.subarray(p, end)]);
    p = end;
    if (c % 4 === 0) { // one ECG chunk of 1 to 9 samples, a raw record that is ignored, and an empty kind 4
      const m = 1 + (c % 9), s = new Float32Array(m).map((_, i) => ecgValue(e + i));
      e += m;
      const t = 1000 + e * ecgMs;
      log.add(t, [2, new Uint8Array(40)], [4, new Uint8Array(s.buffer)]);
      log.add(t, [4, new Uint8Array(0)]);
    }
  }
  log.event(1000 + 200000, { e: 'stop' });
  assert.ok(await log.flush());
  return { files, ecgSamples: e };
}

async function* feed(files) {
  yield* files;
}

test('a recording decodes to its samples and its two rates', async () => {
  const { files, ecgSamples } = await synthetic();
  assert.ok(files.length >= 3);
  const r = await loadRecording(feed(files));
  assert.deepEqual(r.errors, []);
  assert.equal(r.segments, files.length);
  assert.equal(r.drops, 0);
  assert.equal(r.pcg.length, 10 * (FRAMES - 1)); // the stale frame is dropped, and the last frame waits for a next header
  assert.ok(r.pcg.every((v, i) => v === SAMPLE(Math.floor(i / 10), i % 10)));
  assert.equal(r.ecg.length, ecgSamples);
  assert.ok(r.ecg.every((v, i) => v === ecgValue(i)));
  assert.deepEqual(r.events.map((e) => e.e), ['start', 'stop']);
  assert.equal(r.t0, 500);
  assert.ok(Math.abs(r.pcgFit.hz / 1001 - 1) < 5e-6, `pcg ${r.pcgFit.hz}`);
  assert.ok(Math.abs(r.ecgHz / 999.5 - 1) < 5e-6, `ecg ${r.ecgHz}`); // no raw ECG here: the stored trace, through a line
  assert.equal(r.map, null);
  assert.equal(r.axes.pcg.timeAt(12345), 12345); // PCG time in ms is the sample index
  assert.ok(Math.abs(r.axes.ecg.indexAt(r.axes.ecg.timeAt(777)) - 777) < 1e-6);
  const span = r.axes.ecg.timeAt(2000) - r.axes.ecg.timeAt(0); // 2000 ECG samples at 999.5 Hz last 2000 * 1001 / 999.5 PCG ms
  assert.ok(Math.abs(span - 2000 * (1001 / 999.5)) < 0.5, `span ${span}`);
});

// PCG frames cut in odd chunks, one frame missing by seq: every frame is stamped by the chunk that held its last byte, and a time is an index
test('PCG frames are stamped like the loader and a dropped frame stays as 10 NaN', async () => {
  const N = 200, HOLE = 60, order = [-1, ...Array.from({ length: N }, (_, k) => k).filter((k) => k !== HOLE)]; // -1: the stale frame
  const bytes = new Uint8Array(24 * order.length);
  order.forEach((k, f) => bytes.set(k < 0 ? pcgFrame(5, 0) : pcgFrame(100 + k, k), 24 * f));
  const files = [], log = createLog(async (_, blob) => files.push(await blob.arrayBuffer())), byteChunk = [];
  const sizes = [1, 5, 17, 37, 23];
  for (let p = 0, c = 0; p < bytes.length; c++) {
    const end = Math.min(bytes.length, p + sizes[c % sizes.length]);
    log.add(1000 + 3 * c, [1, bytes.subarray(p, end)]);
    byteChunk.push(...Array(end - p).fill(c));
    p = end;
  }
  assert.ok(await log.flush());
  const r = await loadRecording(feed(files));
  assert.deepEqual(r.errors, []);
  assert.equal(r.drops, 1);
  assert.equal(r.pcg.length, 10 * (N - 1)); // the last frame waits for a next header; the hole is in the middle
  assert.equal(r.pcgFrameT.length, N - 2);
  r.pcgFrameT.forEach((t, j) => {
    const k = j + (j + 1 > HOLE ? 1 : 0), f = order.indexOf(k); // frame k: the j-th live one, counting past the hole
    assert.equal(t, 1000 + 3 * byteChunk[24 * f + 23], `frame ${k}`);
    assert.equal(r.pcgFrameIndex[j], 10 * k);
  });
  for (let i = 0; i < r.pcg.length; i++) {
    if (i >= 10 * HOLE && i < 10 * HOLE + 10) assert.ok(Number.isNaN(r.pcg[i]), `index ${i}`);
    else assert.equal(r.pcg[i], SAMPLE(Math.floor(i / 10), i % 10), `index ${i}`);
  }
});

test('a cut segment and a missing number are reported, and nothing throws', async () => {
  const { files } = await synthetic();
  const cut = files.map((f, i) => (i === 1 ? f.slice(0, f.byteLength - 5) : f));
  const r = await loadRecording(feed(cut));
  assert.equal(r.errors.length, 1);
  assert.match(r.errors[0], /file 2: the record at byte \d+ overruns the file/);
  assert.ok(r.events.some((e) => e.e === 'stop')); // the other files still load
  const list = segmentList(['00004.ecgpcg', 'notes.txt', '00001.ecgpcg', 'x.ecgpcg', '00002.ecgpcg']);
  assert.deepEqual(list.order, [2, 4, 0]);
  assert.deepEqual(list.errors, ['x.ecgpcg: not a numbered segment', 'missing segment numbers: 3']);
  assert.deepEqual(segmentList(['a.txt']).errors, ['no segment files']);
});

test('columns keeps a single spike in a million samples', () => {
  const s = new Int16Array(1e6);
  s[777777] = 5000;
  const axis = { indexAt: (t) => t }; // 1 sample per ms of PCG time
  const c = columns(s, axis, 0, 1e6, 1000);
  assert.equal(c.max[777], 5000);
  assert.equal(c.max.reduce((m, v, x) => m + (x === 777 || v === 0 ? 0 : 1), 0), 0);
  assert.ok(c.min.every((v) => v === 0));
  assert.ok(columns(s, axis, 2e6, 3e6, 10).max.every(Number.isNaN)); // outside the data
});

// 40-byte ECG frames (protocol section 8: 01 C0, 02 C0, 03 C0 markers, lead II big-endian at byte 8), each part cut in chunks of odd sizes.
// The stamp of a chunk is 2.002 ms per frame completed so far (junk bytes aside), so the fitted rate is known
const MS_PER_FRAME = 2.002;
function ecgFrames(values) {
  const out = new Uint8Array(40 * values.length);
  values.forEach((v, m) => {
    out.set([1, 0xc0], 40 * m);
    out.set([2, 0xc0], 40 * m + 12);
    out.set([3, 0xc0], 40 * m + 24);
    new DataView(out.buffer).setInt16(40 * m + 8, v);
  });
  return out;
}
async function ecgRecording(...parts) {
  const files = [], log = createLog(async (_, blob) => files.push(await blob.arrayBuffer()));
  const sizes = [1, 5, 17, 37, 23];
  let c = 0, bytes = 0;
  for (const part of parts) {
    for (let p = 0; p < part.length; ) {
      const end = Math.min(part.length, p + sizes[c++ % sizes.length]);
      bytes += end - p;
      log.add(1000 + Math.floor(bytes / 40) * MS_PER_FRAME, [2, part.subarray(p, end)]);
      p = end;
    }
  }
  assert.ok(await log.flush());
  return loadRecording(feed(files));
}
const wave = (m) => Math.round(-3000 + 0.5 * m + 300 * Math.sin((2 * Math.PI * 1.3 * m) / 500) + 100 * Math.sin((2 * Math.PI * 12 * m) / 500) + 40 * Math.sin((2 * Math.PI * 50 * m) / 500));
const frames400 = Array.from({ length: 400 }, (_, m) => wave(m));

test('the ECG is filtered forward and backward like scipy, then upsampled to 1000 Hz', async () => {
  const r = await ecgRecording(ecgFrames(frames400));
  assert.deepEqual(r.errors, []);
  assert.equal(r.ecgSource, 'raw');
  assert.equal(r.ecg.length, 799);
  // scipy 1.11: y = filtfilt(HP), filtfilt(NF), filtfilt(HF) of the 400 counts; c = lfilter(H40, [1], y extended by 11 copies of its first and 6 of its
  // last value)[11:]; linear x2 of c; drop the first 11 samples. The coefficients are the ones in ecg_dsp.js
  const expected = { 0: -30.963467, 1: -24.992767, 2: -17.349284, 3: -9.705801, 400: -103.054814, 401: -102.642974, 555: -263.35002, 556: -266.767549,
    795: 144.242767, 796: 138.623187, 797: 133.003607, 798: 129.122544 };
  for (const [i, v] of Object.entries(expected)) assert.ok(Math.abs(r.ecg[i] - v) < 1e-3, `sample ${i}: ${r.ecg[i]} against ${v}`);
  assert.ok(Math.abs(r.ecgHz / (1000 / (MS_PER_FRAME / 2)) - 1) < 5e-6, `ecg ${r.ecgHz}`);
  assert.equal(r.ecgStamps.length, 400); // one host stamp per frame, rising
});

test('a symmetric pulse comes out symmetric about its own time', async () => {
  // 2000 frames: filtfilt's start-up differs a little between the two ends, and less the longer the record (0.002 counts here, 0.2 at 600 frames)
  const r = await ecgRecording(ecgFrames(Array.from({ length: 2000 }, (_, m) => Math.max(0, 200 - 20 * Math.abs(m - 1000)))));
  // the top of the pulse is flat to 0.003 counts, so the peak sample itself is not pinned down: the symmetry is what a delay would break
  assert.ok(Math.max(...r.ecg.map(Math.abs)) - r.ecg[2000] < 0.01); // frame 1000 is sample 2000
  for (let d = 1; d <= 40; d++) assert.ok(Math.abs(r.ecg[2000 - d] - r.ecg[2000 + d]) < 0.01, `d ${d}`);
});

test('a resync cuts the ECG in two runs joined by one empty sample, each filtered on its own', async () => {
  const bytes = ecgFrames(frames400), junk = new Uint8Array([9, 8, 7, 6, 5, 4, 3]);
  const r = await ecgRecording(bytes.subarray(0, 8000), Uint8Array.of(...junk, ...bytes.subarray(8000)));
  const [a, b] = await Promise.all([ecgRecording(bytes.subarray(0, 8000)), ecgRecording(bytes.subarray(8000))]);
  assert.deepEqual(r.ecgResyncs, [200]);
  assert.equal(r.ecg.length, 799);
  assert.ok(Number.isNaN(r.ecg[399]));
  assert.ok(r.ecgStamps.every((t, m) => m === 0 || t >= r.ecgStamps[m - 1])); // the stamps are per frame and run past the cut
  assert.equal(r.ecgStamps.length, 400);
  assert.equal(r.map.runs.length, 2);
  assert.deepEqual(r.ecg.slice(0, 399), a.ecg);
  assert.deepEqual(r.ecg.slice(400), b.ecg);
});
