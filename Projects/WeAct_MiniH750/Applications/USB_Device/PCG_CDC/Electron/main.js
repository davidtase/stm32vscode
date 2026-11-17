const { app, BrowserWindow, session } = require('electron');
const http = require('http');
const fs = require('fs');
const path = require('path');

const WEB = path.resolve(__dirname, '..', 'Web');
const PORT = 8002;
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript' };

// Serve ../Web read-only over http: the page uses ES modules and audioWorklet.addModule.
const server = http.createServer((req, res) => {
  let file;
  try {
    const name = decodeURIComponent(new URL(req.url, 'http://x').pathname);
    file = path.join(WEB, name);
  } catch { res.writeHead(400).end(); return; }
  if (!file.startsWith(WEB + path.sep)) { res.writeHead(403).end(); return; }
  fs.readFile(file, (err, data) => {
    if (err) { res.writeHead(404).end(); return; }
    res.writeHead(200, { 'Content-Type': TYPES[path.extname(file)] || 'application/octet-stream' }).end(data);
  });
});

// Only the PCG board (0483:5740) and the ECG module (04d8:00df). Electron gives vendorId and productId as decimal strings.
const KNOWN = [[0x0483, 0x5740], [0x04d8, 0x00df]];

app.whenReady().then(() => {
  const ses = session.defaultSession;
  ses.setPermissionCheckHandler(() => true); // same as Electron's default; serial and the folder picker need it
  ses.setDevicePermissionHandler(({ deviceType, device }) =>
    deviceType === 'serial' && KNOWN.some(([v, p]) => Number(device.vendorId) === v && Number(device.productId) === p));
  ses.on('select-serial-port', (event, ports, webContents, callback) => {
    event.preventDefault();
    callback(ports[0]?.portId ?? ''); // '' cancels, so the page's Allow device gets its usual "no device" answer
  });
  server.listen(PORT, '127.0.0.1', () => {
    const win = new BrowserWindow({ width: 1200, height: 800, webPreferences: { backgroundThrottling: false } });
    win.loadURL(`http://127.0.0.1:${PORT}/recorder.html`);
  });
});
app.on('window-all-closed', () => app.quit());

// Recording > View a recording...: the viewer opens in a second window, the recorder window is not touched.
const { Menu, MenuItem } = require('electron');
const openViewer = () => new BrowserWindow({ width: 1400, height: 900 }).loadURL(`http://127.0.0.1:${PORT}/viewer.html`);
app.whenReady().then(() => {
  const menu = Menu.getApplicationMenu(); // the default menu, so Reload and Developer Tools stay
  menu.append(new MenuItem({ label: 'Recording', submenu: [{ label: 'View a recording...', accelerator: 'CmdOrCtrl+O', click: openViewer }] }));
  Menu.setApplicationMenu(menu);
});
