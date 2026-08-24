// 对比新旧结果的凹凸保持: 按 fid 对齐
import { readFileSync } from 'fs';

function parseShp(path) {
  const buf = readFileSync(path);
  const feats = [];
  let off = 100;
  while (off + 8 <= buf.length) {
    const cl = buf.readInt32BE(off + 4) * 2;
    const b = buf.subarray(off + 8, off + 8 + cl);
    if (cl <= 0 || off + 8 + cl > buf.length) break;
    if (b.readInt32LE(0) === 5) {
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

function parseDbfFids(path) {
  const buf = readFileSync(path);
  const recCount = buf.readInt32LE(4);
  const headerLen = buf.readInt16LE(8);
  const recLen = buf.readInt16LE(10);
  let p = 32;
  const fields = [];
  while (buf[p] !== 0x0d) {
    const name = buf.toString('ascii', p, p + 11).replace(/\0+$/, '').trim();
    const type = String.fromCharCode(buf[p + 11]);
    const len = buf[p + 16];
    fields.push({ name, type, len });
    p += 32;
  }
  const fidField = fields.find(f => /fid/i.test(f.name) && f.type === 'N') || fields.find(f => f.type === 'N');
  const fidOff = 32 + fields.slice(0, fields.indexOf(fidField)).reduce((s, f) => s + f.len, 0);
  const fids = [];
  for (let r = 0; r < recCount; r++) {
    const base = headerLen + r * recLen;
    if (base + recLen > buf.length) break;
    fids.push(parseInt(buf.toString('ascii', base + fidOff, base + fidOff + fidField.len).trim(), 10));
  }
  return fids;
}

function area(r) {
  let a = 0;
  for (let i = 0; i < r.length; i++) {
    const j = (i + 1) % r.length;
    a += r[i][0] * r[j][1] - r[j][0] * r[i][1];
  }
  return a / 2;
}

const dir = process.argv[2] || 'test/maskonly_topo';
const oldDir = process.argv[3] || 'test/maskonly_out';
const init = parseShp(`${dir}/initial_building_outline.shp`);
const initFids = parseDbfFids(`${dir}/initial_building_outline.dbf`);
const newR = parseShp(`${dir}/regularized_building.shp`);
const newFids = parseDbfFids(`${dir}/regularized_building.dbf`);
const oldR = parseShp(`${oldDir}/regularized_building.shp`);
const oldFids = parseDbfFids(`${oldDir}/regularized_building.dbf`);

const newByFid = new Map(newFids.map((f, i) => [f, newR[i]]));
const oldByFid = new Map(oldFids.map((f, i) => [f, oldR[i]]));

let preserved = 0, collapsed = 0, bothKeep = 0;
const samples = [];
for (let k = 0; k < init.length; k++) {
  const fid = initFids[k];
  const iv = init[k].rings[0].length;
  const a = Math.abs(area(init[k].rings[0]));
  if (iv < 7 || a < 30 || a > 600) continue;
  const nf = newByFid.get(fid);
  const of_ = oldByFid.get(fid);
  if (!nf || !of_) continue;
  const nv = nf.rings[0].length;
  const ov = of_.rings[0].length;
  if (nv >= 6 && ov <= 4) {
    preserved++;
    if (samples.length < 6) samples.push({ fid, iv, nv, ov, area: a.toFixed(0) });
  } else if (nv <= 4 && ov <= 4) collapsed++;
  else if (nv >= 6 && ov >= 6) bothKeep++;
}
console.log(`notch-like (init>=7 verts, 30-600m2): new-preserves-old-collapsed=${preserved}, both-collapse=${collapsed}, both-keep=${bothKeep}`);
console.log('samples (fid, initV, newV, oldV, area):', JSON.stringify(samples));
