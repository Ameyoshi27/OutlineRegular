// 临时分析: 三阶段几何对比 - 缺口(凹槽)的宽度/深度测量
import { readFileSync } from 'fs';

function parseShpRecords(path) {
  const buf = readFileSync(path);
  const recs = [];
  let off = 100;
  while (off + 8 <= buf.length) {
    const contentLen = buf.readInt32BE(off + 4) * 2;
    const body = buf.subarray(off + 8, off + 8 + contentLen);
    if (contentLen <= 0 || off + 8 + contentLen > buf.length) break;
    const t = body.readInt32LE(0);
    if (t === 5 || t === 15 || t === 25) {
      let p = 36;
      const np = body.readInt32LE(p); p += 4;
      const npts = body.readInt32LE(p); p += 4;
      const parts = [];
      for (let i = 0; i < np; i++) { parts.push(body.readInt32LE(p)); p += 4; }
      parts.push(npts);
      const rings = [];
      for (let i = 0; i < np; i++) {
        const ring = [];
        for (let j = parts[i]; j < parts[i + 1]; j++) { ring.push([body.readDoubleLE(p), body.readDoubleLE(p + 8)]); p += 16; }
        if (ring.length > 1) {
          const a = ring[0], b = ring[ring.length - 1];
          if (Math.abs(a[0]-b[0]) < 1e-9 && Math.abs(a[1]-b[1]) < 1e-9) ring.pop();
        }
        if (ring.length >= 3) rings.push(ring);
      }
      recs.push({ rings, deleted: t === 0 });
    } else recs.push({ rings: null, deleted: true });
    off += 8 + contentLen;
  }
  return recs;
}
function parseDbfRecords(path) {
  const buf = readFileSync(path);
  const n = buf.readUInt16LE(4), hs = buf.readUInt16LE(8), rs = buf.readUInt16LE(10);
  const fields = [];
  let p = 32;
  while (buf[p] !== 0x0D) {
    fields.push({ name: buf.toString('ascii', p, p + 11).replace(/\0.*$/,'').trim(), type: String.fromCharCode(buf[p+11]), len: buf[p+16] });
    p += 32;
  }
  const recs = [];
  for (let r = 0; r < n; r++) {
    const base = hs + r * rs;
    if (base + rs > buf.length) break;
    const rec = { deleted: buf[base] === 0x2A };
    let q = base + 1;
    for (const f of fields) {
      const raw = buf.toString('latin1', q, q + f.len).trim();
      rec[f.name] = f.type === 'N' ? parseFloat(raw) : raw;
      q += f.len;
    }
    recs.push(rec);
  }
  return recs;
}
const dir = 'D:/outlineRegular/outlineRegular/build_deps_release/Release/';
const initRecs = parseShpRecords(dir + 'initial_building_outline.shp');
const initDbf = parseDbfRecords(dir + 'initial_building_outline.dbf');
const active = initRecs.filter(r => r.rings);
const ini = active[1254];
console.log('initial[1254] 存在:', !!ini, ini ? '顶点=' + ini.rings[0].length : '');

const hypRecs = parseShpRecords(dir + 'debug_best_hypothesis.shp').filter(r => r.rings);
const hyp = hypRecs[98];
console.log('hypothesis[98] 存在:', !!hyp, hyp ? '顶点=' + hyp.rings[0].length : '');

const resDbf = parseDbfRecords('D:/outlineRegular/outlineRegular/test/regularized_building.dbf');
const resRecs = parseShpRecords('D:/outlineRegular/outlineRegular/test/regularized_building.shp').filter(r => r.rings);
let res = null;
for (let i = 0; i < resDbf.length; i++) {
  if (resDbf[i].mask === 3283) { res = resRecs[i]; console.log('result mask=3283: rid=' + i + ' 顶点=' + res.rings[0].length); }
}

const area = r => Math.abs(r.reduce((s, p, i) => { const q = r[(i + 1) % r.length]; return s + p[0] * q[1] - q[0] * p[1]; }, 0)) / 2;
const cen = r => { let x=0,y=0; for(const p of r){x+=p[0];y+=p[1];} return [x/r.length,y/r.length]; };

// 凹槽检测: 边相对前后边"回退"的顶点序列。简化: 找所有凹顶点(叉积符号与多数相反), 
// 更直接: 输出顶点序列和边长/转角, 人工判读
function describe(r, name) {
  const c = cen(r);
  console.log(`\n=== ${name}: 面积=${area(r).toFixed(1)} 顶点=${r.length} ===`);
  // 边序列: 长度@角度, 标记凹转(与多数转向相反)
  let turnSum = 0;
  const turns = [];
  for (let i = 0; i < r.length; i++) {
    const a = r[(i + r.length - 1) % r.length], b = r[i], d = r[(i + 1) % r.length];
    const cross = (b[0]-a[0])*(d[1]-b[1]) - (b[1]-a[1])*(d[0]-b[0]);
    turns.push(cross);
    turnSum += cross;
  }
  const majoritySign = turnSum >= 0 ? 1 : -1;
  const edges = [];
  for (let i = 0; i < r.length; i++) {
    const p = r[i], q = r[(i + 1) % r.length];
    const len = Math.hypot(q[0]-p[0], q[1]-p[1]);
    let ang = Math.atan2(q[1]-p[1], q[0]-p[0]) * 180 / Math.PI;
    ang = ((ang % 90) + 90) % 90;
    edges.push(`${len.toFixed(1)}m@${ang.toFixed(0)}°${turns[i]*majoritySign < 0 ? '[凹角]' : ''}`);
  }
  console.log('  边(长度@角度mod90, [凹角]=该边起点是凹顶点):');
  console.log('  ' + edges.join(' '));
}
if (ini) describe(ini.rings[0], 'initial fid=1254 (楼梯,取前40条边)');
if (ini && ini.rings[0].length > 45) {
  // 楼梯太密, 也输出简化说明
  console.log('  (楼梯轮廓共 ' + ini.rings[0].length + ' 顶点, 边太密, 略)');
}
if (hyp) describe(hyp.rings[0], 'best_hypothesis id=98');
if (res) describe(res.rings[0], 'result mask=3283');

// 三者质心互距
if (ini && hyp && res) {
  const c1 = cen(ini.rings[0]), c2 = cen(hyp.rings[0]), c3 = cen(res.rings[0]);
  console.log(`\n质心距离: ini-hyp=${Math.hypot(c1[0]-c2[0],c1[1]-c2[1]).toFixed(1)}m ini-res=${Math.hypot(c1[0]-c3[0],c1[1]-c3[1]).toFixed(1)}m`);
}
