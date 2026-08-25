// 方向接入验收:
// 1) 单方向建筑: 所有长边(>1.5m)必须贴近主系统族(theta/theta+90, 容差8°)
// 2) 多方向建筑: 至少保留两个方向族(按边长加权聚类, 间隔>20°)
// 3) 未解释长斜边: 不属于任何系统族且 >3m 的边
// 4) 汇总自交/面积比在 run console 已有, 此处补边级审计
import { readFileSync } from 'fs';

function parseShp(path) {
  const buf = readFileSync(path);
  const feats = []; let off = 100;
  while (off + 8 <= buf.length) {
    const cl = buf.readInt32BE(off + 4) * 2;
    const b = buf.subarray(off + 8, off + 8 + cl);
    if (cl <= 0 || off + 8 + cl > buf.length) break;
    const t = b.readInt32LE(0);
    if (t === 5 || t === 15 || t === 25) {
      let p = 36; const np = b.readInt32LE(p); p += 4; const npts = b.readInt32LE(p); p += 4;
      const parts = []; for (let i = 0; i < np; i++) { parts.push(b.readInt32LE(p)); p += 4; }
      parts.push(npts);
      const rings = []; for (let i = 0; i < np; i++) {
        const r = []; for (let j = parts[i]; j < parts[i + 1]; j++) {
          r.push([b.readDoubleLE(p), b.readDoubleLE(p + 8)]); p += 16;
        } rings.push(r);
      }
      feats.push({ rings });
    }
    off += 8 + cl;
  }
  return feats;
}
function parseDbfIds(path) {
  const buf = readFileSync(path);
  const recCount = buf.readInt32LE(4), headerLen = buf.readInt16LE(8), recLen = buf.readInt16LE(10);
  let p = 32; const fields = [];
  while (buf[p] !== 0x0d) {
    const name = buf.toString('ascii', p, p + 11).replace(/[\0 ]+$/, '').trim();
    const type = String.fromCharCode(buf[p + 11]); const len = buf[p + 16];
    fields.push({ name, type, len, off: 0 }); p += 32;
  }
  let off = 1; for (const f of fields) { f.off = off; off += f.len; }
  const idf = fields.find(f => f.name === 'id'); const ids = [];
  for (let r = 0; r < recCount; r++) {
    const base = headerLen + r * recLen;
    ids.push(parseInt(buf.toString('ascii', base + idf.off, base + idf.off + idf.len).trim(), 10));
  }
  return ids;
}

const dir = process.argv[2];
const res = parseShp(`${dir}/regularized_building.shp`);
const ids = parseDbfIds(`${dir}/regularized_building.dbf`);
const byId = new Map(ids.map((id, i) => [id, res[i]]));
const log = readFileSync(`${dir}/run_console.txt`, 'utf8').split('\n');
const modeByFid = {}, sysByFid = {};
for (const l of log) {
  let m = l.match(/\[DirectionApply\] fid=(\d+) mode=(\w+)/);
  if (m) modeByFid[m[1]] = m[2];
  if (l.includes('systems=') && l.includes('angle_deg=')) {
    const fm = l.match(/fid=(\d+)/);
    const angs = [...l.matchAll(/angle_deg=([0-9.]+)/g)].map(x => +x[1]);
    if (fm && angs.length) sysByFid[fm[1]] = angs;
  }
}
// 边到方向系统族的距离(展开角, mod 90)
const distToSystem = (ang, sysDeg) => {
  let d = Math.abs(((ang - sysDeg) % 90 + 90) % 90);
  if (d > 45) d = 90 - d;
  return d;
};
function allEdges(f) {
  const out = [];
  for (const r of f.rings) for (let i = 0; i < r.length; i++) {
    const a = r[i], b = r[(i + 1) % r.length];
    const len = Math.hypot(b[0] - a[0], b[1] - a[1]);
    out.push({ len, ang: Math.atan2(b[1] - a[1], b[0] - a[0]) * 180 / Math.PI });
  }
  return out;
}
const fold = a => { a = ((a % 90) + 90) % 90; return a > 45 ? 90 - a : a; };

let singleViolations = [], unexplained = [], multiOK = 0, multiFail = [];
let singleCount = 0, multiCount = 0;
for (const [fidS, mode] of Object.entries(modeByFid)) {
  const fid = +fidS;
  const f = byId.get(fid); if (!f) continue;
  const sys = sysByFid[fidS]; if (!sys) continue;
  const edges = allEdges(f);
  if (mode === 'single') {
    singleCount++;
    for (const e of edges) {
      if (e.len < 1.5) continue;
      if (distToSystem(e.ang, sys[0]) > 8) {
        singleViolations.push({ fid, len: e.len, dev: distToSystem(e.ang, sys[0]) });
        break;
      }
    }
  } else if (mode === 'multi') {
    multiCount++;
    // 边长加权的方向族: 长边折叠角聚类(间隔>20°)
    const long = edges.filter(e => e.len >= 1.5).map(e => fold(e.ang)).sort((a, b) => a - b);
    const cl = [];
    for (const a of long) {
      if (!cl.length || a - cl[cl.length - 1] > 20) cl.push(a);
    }
    if (cl.length >= 2) multiOK++;
    else multiFail.push({ fid, clusters: cl });
  }
  // 未解释长斜边(所有模式): 不属任何系统族且 >3m
  for (const e of edges) {
    if (e.len < 3.0) continue;
    const best = Math.min(...sys.map(s => distToSystem(e.ang, s)));
    if (best > 10) unexplained.push({ fid, len: e.len, dev: best, mode });
  }
}
console.log(`单方向建筑: ${singleCount}, 违规(>1.5m边偏离主系统>8°): ${singleViolations.length}`);
singleViolations.slice(0, 10).forEach(v =>
  console.log(`  fid=${v.fid} len=${v.len.toFixed(1)}m dev=${v.dev.toFixed(1)}°`));
console.log(`多方向建筑: ${multiCount}, 输出保留>=2方向族: ${multiOK}, 不达标: ${multiFail.length}`);
multiFail.forEach(v => console.log(`  fid=${v.fid} 簇=${v.clusters.map(c => c.toFixed(0)).join('/')}`));
console.log(`未解释长斜边(>3m, 偏离所有系统>10°): ${unexplained.length} 条`);
unexplained.slice(0, 10).forEach(u =>
  console.log(`  fid=${u.fid} len=${u.len.toFixed(1)}m dev=${u.dev.toFixed(1)}° (${u.mode})`));
