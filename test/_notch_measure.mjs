// 临时分析: 简化楼梯轮廓并测量缺口宽深
import { readFileSync } from 'fs';

function parseShp(path) {
  const buf = readFileSync(path);
  const recs = [];
  let off = 100;
  while (off + 8 <= buf.length) {
    const cl = buf.readInt32BE(off + 4) * 2;
    const b = buf.subarray(off + 8, off + 8 + cl);
    if (cl <= 0 || off + 8 + cl > buf.length) break;
    const t = b.readInt32LE(0);
    if (t === 5 || t === 15 || t === 25) {
      let p = 36;
      const np = b.readInt32LE(p); p += 4;
      const npts = b.readInt32LE(p); p += 4;
      const parts = [];
      for (let i = 0; i < np; i++) { parts.push(b.readInt32LE(p)); p += 4; }
      parts.push(npts);
      const ring = [];
      for (let j = parts[0]; j < parts[1]; j++) { ring.push([b.readDoubleLE(p), b.readDoubleLE(p + 8)]); p += 16; }
      if (ring.length > 2) recs.push({ ring });
    }
    off += 8 + cl;
  }
  return recs;
}

const dir = 'D:/outlineRegular/outlineRegular/build_deps_release/Release/';
const ini = parseShp(dir + 'initial_building_outline.shp').filter(r => r.ring)[1254];
const ring = [...ini.ring];
if (ring.length > 1 && Math.hypot(ring[0][0]-ring[ring.length-1][0], ring[0][1]-ring[ring.length-1][1]) < 1e-9) ring.pop();
ini.ring = ring;
console.log('initial[1254] 顶点:', ini.ring.length);

function dp(points, tol) {
  if (points.length < 3) return points;
  const keep = new Array(points.length).fill(false);
  keep[0] = keep[points.length - 1] = true;
  const stack = [[0, points.length - 1]];
  while (stack.length) {
    const [a, b] = stack.pop();
    let maxD = -1, idx = -1;
    const ax = points[a][0], ay = points[a][1], bx = points[b][0], by = points[b][1];
    const dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy) || 1e-12;
    for (let i = a + 1; i < b; i++) {
      const d = Math.abs((points[i][0] - ax) * dy - (points[i][1] - ay) * dx) / L;
      if (d > maxD) { maxD = d; idx = i; }
    }
    if (maxD > tol) { keep[idx] = true; stack.push([a, idx], [idx, b]); }
  }
  return points.filter((_, i) => keep[i]);
}

const simp = dp(ini.ring, 0.45);
console.log('DP(0.45m) 简化后顶点:', simp.length);

const n = simp.length;
let turnSum = 0;
const turns = [];
for (let i = 0; i < n; i++) {
  const a = simp[(i + n - 1) % n], b = simp[i], c = simp[(i + 1) % n];
  const cr = (b[0]-a[0])*(c[1]-b[1]) - (b[1]-a[1])*(c[0]-b[0]);
  turns.push(cr); turnSum += cr;
}
const sgn = turnSum >= 0 ? 1 : -1;
const reflex = [];
for (let i = 0; i < n; i++) if (turns[i] * sgn < 0) reflex.push(i);
console.log('凹顶点:', reflex.length, '位置:', reflex.join(','));

for (let k = 0; k < reflex.length; k++) {
  const i = reflex[k];
  const j = reflex[(k + 1) % reflex.length];
  const gap = (j - i + n) % n;
  if (gap >= 1 && gap <= 2) {
    let bottom = 0;
    for (let s = i; s !== j; s = (s + 1) % n)
      bottom += Math.hypot(simp[(s+1)%n][0]-simp[s][0], simp[(s+1)%n][1]-simp[s][1]);
    const inLen = Math.hypot(simp[i][0]-simp[(i+n-1)%n][0], simp[i][1]-simp[(i+n-1)%n][1]);
    const outLen = Math.hypot(simp[(j+1)%n][0]-simp[j][0], simp[(j+1)%n][1]-simp[j][1]);
    console.log(`缺口@v${i}: 宽=${bottom.toFixed(1)}m 入深=${inLen.toFixed(1)}m 出深=${outLen.toFixed(1)}m`);
  }
}
// 边列表
const edges = simp.map((p, i) => {
  const q = simp[(i+1)%n];
  return Math.hypot(q[0]-p[0], q[1]-p[1]).toFixed(1);
});
console.log('简化边长序列:', edges.join(' '));
