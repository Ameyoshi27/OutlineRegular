// 诊断: 为何指定 id 未在 GeometryOnly 合并中合并
// 输出每个 id 的 parent/面积/边界, 以及指定对之间的缝隙距离
import { readFileSync } from 'fs';

const dir = process.argv[2] || 'test';
const wantIds = (process.argv[3] || '133,148,134,138,240,263,334,260').split(',').map(Number);
const pairs = process.argv[4] ? process.argv[4].split(';').map(p => p.split(',').map(Number)) : [[133,148],[134,138]];

function parseShp(path) {
  const buf = readFileSync(path);
  const feats = [];
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
      const rings = [];
      for (let i = 0; i < np; i++) {
        const r = [];
        for (let j = parts[i]; j < parts[i + 1]; j++) {
          r.push([b.readDoubleLE(p), b.readDoubleLE(p + 8)]);
          p += 16;
        }
        rings.push(r);
      }
      feats.push({ rings });
    }
    off += 8 + cl;
  }
  return feats;
}

function parseDbf(path) {
  const buf = readFileSync(path);
  const recCount = buf.readInt32LE(4);
  const headerLen = buf.readInt16LE(8);
  const recLen = buf.readInt16LE(10);
  let p = 32;
  const fields = [];
  while (buf[p] !== 0x0d) {
    const name = buf.toString('ascii', p, p + 11).replace(/[\0 ]+$/, '').trim();
    const type = String.fromCharCode(buf[p + 11]);
    const len = buf[p + 16];
    fields.push({ name, type, len, off: 0 });
    p += 32;
  }
  let off = 1;
  for (const f of fields) { f.off = off; off += f.len; }
  const recs = [];
  for (let r = 0; r < recCount; r++) {
    const base = headerLen + r * recLen;
    if (base + recLen > buf.length) break;
    const rec = {};
    for (const f of fields) {
      rec[f.name] = buf.toString('ascii', base + f.off, base + f.off + f.len).trim();
    }
    recs.push(rec);
  }
  return recs;
}

const feats = parseShp(`${dir}/initial_building_outline.shp`);
const dbf = parseDbf(`${dir}/initial_building_outline.dbf`);
console.log(`features=${feats.length} dbf=${dbf.length}`);

// id 字段为矢量分配时的稳定顺序号
function realId(rec) {
  return parseInt(rec.id, 10);
}

const byId = new Map();
for (let k = 0; k < Math.min(feats.length, dbf.length); k++) {
  const id = realId(dbf[k]);
  byId.set(id, { idx: k, ...dbf[k], ring: feats[k].rings[0], rings: feats[k].rings });
}

function area(r) {
  let a = 0;
  for (let i = 0; i < r.length; i++) {
    const j = (i + 1) % r.length;
    a += r[i][0] * r[j][1] - r[j][0] * r[i][1];
  }
  return Math.abs(a / 2);
}
function bbox(r) {
  let x0 = 1e18, y0 = 1e18, x1 = -1e18, y1 = -1e18;
  for (const p of r) { x0 = Math.min(x0, p[0]); y0 = Math.min(y0, p[1]); x1 = Math.max(x1, p[0]); y1 = Math.max(y1, p[1]); }
  return [x0, y0, x1, y1];
}
// 两环最小距离(边-边)
function segSegDist(p1, p2, p3, p4) {
  const d = (a, b, c) => (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1]];
  const norm = a => Math.hypot(a[0], a[1]);
  const proj = (p, a, b) => {
    const ab = sub(b, a), ap = sub(p, a);
    const t = Math.max(0, Math.min(1, dot(ap, ab) / Math.max(dot(ab, ab), 1e-12)));
    return [a[0] + t * ab[0], a[1] + t * ab[1]];
  };
  const ptSeg = (p, a, b) => Math.hypot(...sub(p, proj(p, a, b)));
  if (d(p1, p2, p3) * d(p1, p2, p4) < 0 && d(p3, p4, p1) * d(p3, p4, p2) < 0) return 0;
  return Math.min(ptSeg(p1, p3, p4), ptSeg(p2, p3, p4), ptSeg(p3, p1, p2), ptSeg(p4, p1, p2));
}
function ringDist(r1, r2) {
  let best = 1e18;
  for (let i = 0; i < r1.length; i++)
    for (let j = 0; j < r2.length; j++) {
      best = Math.min(best, segSegDist(r1[i], r1[(i + 1) % r1.length], r2[j], r2[(j + 1) % r2.length]));
      if (best === 0) return 0;
    }
  return best;
}
// 共享边长度(近似): 距离<0.05m 的边对, 投影重叠长度
function sharedEdgeLen(r1, r2, tol = 0.08) {
  let len = 0;
  for (let i = 0; i < r1.length; i++) {
    const a = r1[i], b = r1[(i + 1) % r1.length];
    const e = Math.hypot(b[0] - a[0], b[1] - a[1]);
    if (e < 1e-9) continue;
    for (let j = 0; j < r2.length; j++) {
      const c = r2[j], d = r2[(j + 1) % r2.length];
      // 中点到对面环的最小距离
      const m = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];
      let dm = 1e18;
      for (let q = 0; q < r2.length; q++) {
        const p1 = r2[q], p2 = r2[(q + 1) % r2.length];
        dm = Math.min(dm, segSegDist(m, m, p1, p2));
      }
      if (dm < tol) { len += e; break; }
    }
  }
  return len;
}

for (const id of wantIds) {
  const f = byId.get(id);
  if (!f) { console.log(`id=${id}: NOT FOUND`); continue; }
  const [x0, y0, x1, y1] = bbox(f.ring);
  console.log(`id=${id} parent=${f.parent} mask=${f.mask} verts=${f.ring.length} area=${area(f.ring).toFixed(1)} bbox=[${x0.toFixed(1)},${y0.toFixed(1)} ~ ${x1.toFixed(1)},${y1.toFixed(1)}] w=${(x1-x0).toFixed(1)} h=${(y1-y0).toFixed(1)}`);
}

console.log('\n--- 指定对分析 ---');
for (const [ida, idb] of pairs) {
  const fa = byId.get(ida), fb = byId.get(idb);
  if (!fa || !fb) { console.log(`${ida}-${idb}: missing`); continue; }
  const d = ringDist(fa.ring, fb.ring);
  const sl = sharedEdgeLen(fa.ring, fb.ring);
  console.log(`${ida}(parent=${fa.parent}) - ${idb}(parent=${fb.parent}): 最小距离=${d.toFixed(3)}m 近似共享边长=${sl.toFixed(2)}m 同parent=${fa.parent === fb.parent && +fa.parent > 0}`);
}

// 对每个单独提到的 id, 找全库最近的其它要素
console.log('\n--- 各 id 最近邻 ---');
const allIds = [...byId.keys()];
for (const id of wantIds) {
  const fa = byId.get(id);
  if (!fa) continue;
  const [x0, y0, x1, y1] = bbox(fa.ring);
  let bestId = -1, bestD = 1e18, bestParent = '';
  for (const oid of allIds) {
    if (oid === id) continue;
    const fb = byId.get(oid);
    const [ox0, oy0, ox1, oy1] = bbox(fb.ring);
    if (ox0 > x1 + 50 || x0 > ox1 + 50 || oy0 > y1 + 50 || y0 > oy1 + 50) continue;
    const d = ringDist(fa.ring, fb.ring);
    if (d < bestD) { bestD = d; bestId = oid; bestParent = fb.parent; }
  }
  console.log(`id=${id}(parent=${fa.parent}) 最近: id=${bestId}(parent=${bestParent}) 距离=${bestD.toFixed(3)}m`);
}
