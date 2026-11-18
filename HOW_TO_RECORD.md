# How to record ECG + PCG

The recorder page records ECG lead II and PCG together, with one Start button. It runs in Chrome or Edge, because it needs Web Serial.

## 1. Connect

1. Plug the WeAct STM32H750 PCG board into USB. It shows up as `0483:5740`.
2. Plug the Biomedic ECG module into USB. It shows up as `04d8:00df`.
3. Put the electrodes on: RA, LL and RL for lead II. Connect the PCG PCB output to PA6, with the ground wire connected first.

To check that both devices are there: `lsusb | grep -E "0483:5740|04d8:00df"`. You should see two lines.

## 2. Start the page server

From the repo root:

```
python3 -m http.server -b 127.0.0.1 -d Projects/WeAct_MiniH750/Applications/USB_Device/PCG_CDC/Web 8001
```

Leave it running. Port 8000 is taken by `strategy_app.py`, so use 8001. Stop the server with Ctrl+C when you've finished.

## 3. Record

1. Open http://localhost:8001/recorder.html in Chrome. After a code change, hard-reload with Ctrl+Shift+R.
2. The first time only, click **Allow device** and choose each device.
3. Pick **Record** to save, because the page starts in **Measure**, which shows the chosen signals (ECG, PCG or Both) live and saves nothing. Both choices are remembered, and locked while running. Click **Start** and choose a folder **outside the repo**, for example `~/Documents/recordings`. Each recording becomes a new `ecgpcg-YYYYMMDD-HHMMSS` folder in it.
4. Optional controls:
   - **Window** sets how many seconds are on screen, from 1 to 10.
   - **Start audio** plays the PCG on headphones. **Volume** goes from 0 to 300 %. Start with the computer volume low, and use headphones only, because a speaker can reach the piezo.
5. Click **Stop** and wait for "Saved …".

Keep the tab in view and don't let the PC sleep during a recording. Don't have `ecg.html` or `index.html` connected at the same time, because each port can be open in only one page.

## 4. Get the four files

From `Projects/WeAct_MiniH750/Applications/USB_Device/PCG_CDC`:

```
python3 -I -B Tools/load_recording.py "<recording folder>"            # summary only
python3 -I -B Tools/load_recording.py "<recording folder>" --export   # writes <folder>-export
```

The export holds `pcg.npy`, `ecg_raw.npy`, `ecg_filtered.npy` (float32, 1000 Hz except `ecg_raw` at 500 Hz) and `timestamps.npz`. The recording folder stays the master copy.

ECG and PCG run on separate clocks, so correct the drift offline with a piecewise fit on `pcg_frame_t` and `ecg_frame_t`. A single straight line is not enough.

## If something goes wrong

- **A device is not found:** replug it and click Start again.
- **"Not started: … already open":** close the other page that has the port open.
- **The page crashed:** the files saved so far are kept, about 1.6 s is lost, and the next Start works normally.

More detail is in `Projects/WeAct_MiniH750/Applications/USB_Device/PCG_CDC/README.md` and `Acquisition_systems_summary.md`.

## 5. View a recording

From `Projects/WeAct_MiniH750/Applications/USB_Device/PCG_CDC/Electron`, run `pnpm start`. Choose **Recording > View a recording...** (Ctrl+O) and pick your recording folder. A second window opens, so a running Record is not disturbed. No board is needed, and the folder is only read.

- ECG (lead II) is on top and PCG below, on one time axis. The ECG is the raw 500 Hz lead II filtered forward and backward on the computer (QRS on time, ST segments flattened), shown at 1000 Hz. Opening takes about 8-10 s per 10 minutes of recording.
- The mouse wheel zooms about the cursor, dragging pans, and a double-click shows the whole recording. The readout under the plots gives the time and the value of each trace at the cursor.
- The time axis is PCG seconds, as in the loader's `pcg.npy` (PCG sample number divided by 1000, 0 s = the first live PCG sample). The ECG is placed on it with its own host time stamps in 30 s blocks, which removes the slow drift between the two clocks. The check line (S1 against R, 60 s blocks) is a warning, not a correction. The fixed ECG-to-PCG offset (about +-3 ms) is not measured and not corrected, and the loader's convention (a PCG frame's first sample is stamped when its last byte arrives) still applies.

**Labeling R and T-end.** The page proposes an R peak and a T-end for every beat. They are only proposals: you must look at each one and correct the wrong ones. Flagged beats (marked "?") are the doubtful ones.

- Colours: orange dashed = proposed, green = confirmed, blue = edited (you moved it), purple = added by you. A deleted label is not drawn but stays in the file as `deleted`.
- Keys: `n` / `p` next / previous beat, `u` next unreviewed, `f` next flagged, Enter confirms the beat on screen and goes to the next unreviewed one. Click a label to select it, then drag it, or use the arrow keys (1 ms, Shift 10 ms). Tab swaps between R and T, Esc lets go, Delete removes it, Ctrl+Z / Ctrl+Y undo and redo. Shift+click adds an R, Ctrl+click sets the T-end of the beat on its left (these two and the click-to-select work with the mark mode off).
- Mark mode: the **Mark mode** button (key `m`, off when the page opens) lets the mouse correct marks of the kind chosen with `r` or `t`: click an R (or T) to remove it, drag it to move it, click anywhere else to add one; drag on empty space still pans, and it works with under 15 s on screen. Off, a click selects a label and a drag moves it, Shift+click adds an R, Ctrl+click sets a T-end. Ctrl+Z undoes any of it.
- **Save labels** downloads `<recording folder name>.labels.csv`. Put it beside the recording folder, never inside it. **Load labels** reads it back, and **Propose again** replaces only labels you have not touched. Changes are also kept in the browser automatically, and the page offers to restore them next time.
- The file has five `#` header lines (map, method, PCG time convention) and then the columns `beat,kind,t_pcg_s,state,doubt`: `kind` is R or T, the time is in PCG seconds. In Python, plain `numpy.loadtxt` stops at the column line; this works: `numpy.genfromtxt(path, skip_header=5, names=True, dtype=None, delimiter=',', encoding='utf-8')`. pandas was not tried.
