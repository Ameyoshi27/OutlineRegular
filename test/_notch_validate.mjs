// 临时验证: 在 hypothesis id=98 上模拟 DetectEvidenceBackedNotches + 支撑点证据
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
const hyp = parseShp(dir + 'debug_best_hypothesis.shp').filter(r => r.ring)[98];
const poly = hyp.ring;
const n = poly.length;
console.log('hypothesis[98] 顶点:', n);

// 与 C++ 相同的缺口检测 (minDepth = max(1.0, 0.6*2.47) = 1.48, maxWidth=15)
const minDepth = 1.48, maxWidth = 15.0;
let turnSum = 0;
const turns = [];
for (let i = 0; i < n; i++) {
  const a = poly[(i + n - 1) % n], b = poly[i], c = poly[(i + 1) % n];
  const cr = (b[0]-a[0])*(c[1]-b[1]) - (b[1]-a[1])*(c[0]-b[0]);
  turns.push(cr); turnSum += cr;
}
const majority = turnSum >= 0 ? 1 : -1;
const reflex = [];
for (let i = 0; i < n; i++) if (turns[i] * majority < 0) reflex.push(i);
console.log('凹顶点:', reflex.join(','));

const distPtSeg = (p, a, b) => {
  const dx = b[0]-a[0], dy = b[1]-a[1];
  const L2 = dx*dx+dy*dy || 1e-12;
  let t = ((p[0]-a[0])*dx+(p[1]-a[1])*dy)/L2;
  t = Math.max(0, Math.min(1, t));
  return Math.hypot(p[0]-a[0]-t*dx, p[1]-a[1]-t*dy);
};

const notches = [];
for (let k = 0; k < reflex.length; k++) {
  const entry = reflex[k], exit = reflex[(k+1)%reflex.length];
  const gap = (exit - entry + n) % n;
  if (gap === 0 || gap > 2) continue;
  let width = 0;
  for (let s = entry; s !== exit; s = (s+1)%n)
    width += Math.hypot(poly[(s+1)%n][0]-poly[s][0], poly[(s+1)%n][1]-poly[s][1]);
  const wallIn = Math.hypot(poly[entry][0]-poly[(entry+n-1)%n][0], poly[entry][1]-poly[(entry+n-1)%n][1]);
  const wallOut = Math.hypot(poly[(exit+1)%n][0]-poly[exit][0], poly[(exit+1)%n][1]-poly[exit][1]);
  const depth = Math.min(wallIn, wallOut);
  if (depth < minDepth || width < 0.5 || width > maxWidth) continue;
  notches.push({ entry, exit, width, depth });
  console.log(`候选缺口@v${entry}-v${exit}: 宽=${width.toFixed(1)}m 深=${depth.toFixed(1)}m`);
}
if (!notches.length) console.log('未检测到缺口!');

// 加载支撑点LAS, 统计每个缺口3条边 0.8m 内的点数
const lasPath = dir + 'debug_support_points.las';
const buf = readFileSync(lasPath);
const offsetToPts = buf.readUInt32LE(96);
const ptFmt = buf.readUInt8(104);
const ptLen = buf.readUInt16LE(105);
let cnt = buf.readUInt32LE(107);
if (cnt === 0) cnt = Number(buf.readBigUInt64LE(247));
const sx = buf.readDoubleLE(131), sy = buf.readDoubleLE(139);
const ox = buf.readDoubleLE(155), oy = buf.readDoubleLE(163);
console.log(`LAS: fmt=${ptFmt} len=${ptLen} count=${cnt}`);

// 建筑bbox预筛
let bx1=1e18,by1=1e18,bx2=-1e18,by2=-1e18;
for (const p of poly) { bx1=Math.min(bx1,p[0]);by1=Math.min(by1,p[1]);bx2=Math.max(bx2,p[0]);by2=Math.max(by2,p[1]); }
const local = [];
for (let i = 0; i < cnt; i++) {
  const p = offsetToPts + i * ptLen;
  const x = buf.readInt32LE(p) * sx + ox, y = buf.readInt32LE(p + 4) * sy + oy;
  if (x < bx1 - 5 || x > bx2 + 5 || y < by1 - 5 || y > by2 + 5) continue;
  local.push([x, y]);
}
console.log('建筑±5m内支撑点:', local.length);

for (const notch of notches) {
  let sc = 0;
  for (const p of local) {
    let best = 1e18;
    for (let s = (notch.entry + n - 1) % n; ; s = (s + 1) % n) {
      best = Math.min(best, distPtSeg(p, poly[s], poly[(s + 1) % n]));
      if (s === (notch.exit + 1) % n) break;
    }
    if (best <= 0.8) sc++;
  }
  console.log(`缺口@v${notch.entry}: 宽=${notch.width.toFixed(1)} 深=${notch.depth.toFixed(1)} 支撑点=${sc} → ${sc >= 8 ? '受保护' : '证据不足'}`);
}
