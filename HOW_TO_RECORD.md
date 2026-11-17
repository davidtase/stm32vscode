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

- ECG (lead II) is on top and PCG below, on one time axis in seconds from the start. The ECG is the raw 500 Hz lead II filtered forward and backward on the computer (QRS on time, ST segments flattened), shown at 1000 Hz. Opening takes about 10 s per 10 minutes of recording.
- The mouse wheel zooms about the cursor, dragging pans, and a double-click shows the whole recording. The readout under the plots gives the time and the value of each trace at the cursor.
- The PCG time can be up to about 10 ms late, and each trace is placed with one straight line, so this is for looking, not for aligning. Use the loader and the export in section 4 for that.
