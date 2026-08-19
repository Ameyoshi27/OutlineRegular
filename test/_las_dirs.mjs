// 临时分析:用点对方向直方图(类Hough, mod 90°)估计支撑点云的墙面主方向
import { readFileSync } from 'fs';

const path = process.argv[2];
const buf = readFileSync(path);
const offsetToPts = buf.readUInt32LE(96);
const ptFmtId = buf.readUInt8(104);
const ptLen = buf.readUInt16LE(105);
let n = buf.readUInt32LE(107);
if (n === 0) n = Number(buf.readBigUInt64LE(247));
const sx = buf.readDoubleLE(131), sy = buf.readDoubleLE(139);
const ox = buf.readDoubleLE(155), oy = buf.readDoubleLE(163);
console.log(`ptFmt=${ptFmtId} ptLen=${ptLen} count=${n}`);

const xs = new Float64Array(n), ys = new Float64Array(n);
let minX = 1e18, minY = 1e18, maxX = -1e18, maxY = -1e18;
for (let i = 0; i < n; i++) {
  const p = offsetToPts + i * ptLen;
  xs[i] = buf.readInt32LE(p) * sx + ox;
  ys[i] = buf.readInt32LE(p + 4) * sy + oy;
  if (xs[i] < minX) minX = xs[i]; if (xs[i] > maxX) maxX = xs[i];
  if (ys[i] < minY) minY = ys[i]; if (ys[i] > maxY) maxY = ys[i];
}
console.log(`range x[${minX.toFixed(1)},${maxX.toFixed(1)}] y[${minY.toFixed(1)},${maxY.toFixed(1)}]`);

// 网格哈希,半径R内的点对方向直方图
const cell = 1.0, R = 2.5, R2 = R * R;
const grid = new Map();
const key = (cx, cy) => cx * 100000 + cy;
for (let i = 0; i < n; i++) {
  const k = key(Math.floor(xs[i] / cell), Math.floor(ys[i] / cell));
  let a = grid.get(k); if (!a) { a = []; grid.set(k, a); }
  a.push(i);
}
const bins = new Float64Array(90);
const kc = Math.ceil(R / cell);
let pairs = 0;
for (let i = 0; i < n; i++) {
  const cx = Math.floor(xs[i] / cell), cy = Math.floor(ys[i] / cell);
  for (let dx = -kc; dx <= kc; dx++) for (let dy = -kc; dy <= kc; dy++) {
    const a = grid.get(key(cx + dx, cy + dy));
    if (!a) continue;
    for (const j of a) {
      if (j <= i) continue;
      const ddx = xs[j] - xs[i], ddy = ys[j] - ys[i];
      const d2 = ddx * ddx + ddy * ddy;
      if (d2 > R2 || d2 < 0.01) continue;
      let ang = Math.atan2(ddy, ddx) * 180 / Math.PI;
      ang = ((ang % 90) + 90) % 90;
      bins[Math.min(89, Math.floor(ang))] += 1;
      pairs++;
    }
  }
}
console.log(`pairs=${pairs}`);
// 平滑+找峰值
const smooth = new Float64Array(90);
for (let i = 0; i < 90; i++)
  for (let d = -2; d <= 2; d++) smooth[i] += bins[((i + d) % 90 + 90) % 90];
let total = 0; for (const v of smooth) total += v;
const sorted = Array.from(smooth).map((v, i) => [v, i]).sort((a, b) => b[0] - a[0]);
console.log('方向直方图峰(mod 90°, 前8个, 含占比):');
for (const [v, i] of sorted.slice(0, 8)) console.log(`  ${i}°-${i + 1}°: ${(v / total * 100).toFixed(2)}%`);
// 检查特定角
for (const t of [0, 9, 13, 23, 32]) {
  let s = 0; for (let d = -2; d <= 2; d++) s += bins[((t + d) % 90 + 90) % 90];
  console.log(`  ±2°@${t}°: ${(s / total * 100).toFixed(2)}%`);
}
