// 临时分析: 定位 initial fid 与 result fid, 判断旋转与路径
import { readFileSync } from 'fs';

function parseShpRecords(path) {
  const buf = readFileSync(path);
  const recs = [];
  let off = 100;
  while (off + 8 <= buf.length) {
    const contentLen = buf.readInt32BE(off + 4) * 2;
    const body = buf.subarray(off + 8, off + 8 + contentLen);
    if (contentLen <= 0 || off + 8 + contentLen > buf.length) break;
    if (body.readInt32LE(0) === 5 || body.readInt32LE(0) === 15 || body.readInt32LE(0) === 25) {
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
      recs.push({ rings, deleted: body.readInt32LE(0) === 0 });
    } else {
      recs.push({ rings: null, deleted: true });
    }
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
    const deleted = buf[base] === 0x2A; // '*'
    const rec = { deleted };
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

// OGR 语义: fid = 未删除记录中的序号(0基)
function ogrFeatures(shpRecs, dbfRecs) {
  const feats = [];
  let shpIdx = 0;
  for (let r = 0; r < dbfRecs.length && shpIdx < shpRecs.length; r++, shpIdx++) {
    const deleted = dbfRecs[r].deleted || shpRecs[shpIdx].deleted || !shpRecs[shpIdx].rings;
    if (!deleted) feats.push({ fid: feats.length, rings: shpRecs[shpIdx].rings, attrs: dbfRecs[r] });
  }
  return feats;
}

const dir = 'D:/outlineRegular/outlineRegular/build_deps_release/Release/';
const init = ogrFeatures(parseShpRecords(dir + 'initial_building_outline.shp'), parseDbfRecords(dir + 'initial_building_outline.dbf'));
const res = ogrFeatures(parseShpRecords('regularized_building.shp'), parseDbfRecords('regularized_building.dbf'));
console.log(`initial 有效要素: ${init.length}, result 有效要素: ${res.length}`);

const area = r => Math.abs(r.reduce((s, p, i) => { const q = r[(i + 1) % r.length]; return s + p[0] * q[1] - q[0] * p[1]; }, 0)) / 2;
const bboxArea = r => { let a=1e18,b=1e18,c=-1e18,d=-1e18; for(const[x,y]of r){a=Math.min(a,x);b=Math.min(b,y);c=Math.max(c,x);d=Math.max(d,y);} return (c-a)*(d-b); };
const cen = r => { let x=0,y=0; for(const p of r){x+=p[0];y+=p[1];} return [x/r.length,y/r.length]; };
const pcaAngle = r => { const c=cen(r); let sxx=0,syy=0,sxy=0; for(const p of r){const dx=p[0]-c[0],dy=p[1]-c[1];sxx+=dx*dx;syy+=dy*dy;sxy+=dx*dy;} return ((0.5*Math.atan2(2*sxy,sxx-syy)*180/Math.PI)%90+90)%90; };
function domAngle(r) {
  const bins = new Array(90).fill(0); let total = 0;
  for (let i = 0; i < r.length; i++) {
    const [x1,y1]=r[i],[x2,y2]=r[(i+1)%r.length];
    const len = Math.hypot(x2-x1,y2-y1); if (len < 0.3) continue;
    const a = ((Math.atan2(y2-y1,x2-x1)*180/Math.PI)%90+90)%90;
    bins[Math.floor(a)] += len; total += len;
  }
  let best=-1,bi=0; bins.forEach((v,i)=>{if(v>best){best=v;bi=i;}});
  return { angle: bi, conf: total>0?best/total:0 };
}

const ini = init.find(f => f.fid === 2473);
const out = res.find(f => f.fid === 189);
if (!ini || !out) { console.log('未找到要素!', !!ini, !!out); process.exit(1); }
const ir = ini.rings[0], orr = out.rings[0];
console.log(`\ninitial fid=2473: mask=${ini.attrs.mask} 面积=${area(ir).toFixed(1)} bbox面积=${bboxArea(ir).toFixed(1)} 顶点=${ir.length} PCA主轴=${pcaAngle(ir).toFixed(1)}° 边直方图=${domAngle(ir).angle}°(conf ${(domAngle(ir).conf*100).toFixed(0)}%)`);
console.log(`result  fid=189 : 面积=${area(orr).toFixed(1)} bbox面积=${bboxArea(orr).toFixed(1)} 顶点=${orr.length} PCA主轴=${pcaAngle(orr).toFixed(1)}° 边直方图=${domAngle(orr).angle}°(conf ${(domAngle(orr).conf*100).toFixed(0)}%)`);
const c1 = cen(ir), c2 = cen(orr);
console.log(`质心距离: ${Math.hypot(c1[0]-c2[0],c1[1]-c2[1]).toFixed(2)} m`);
let dd = Math.abs(domAngle(ir).angle - domAngle(orr).angle); if (dd > 45) dd = 90 - dd;
console.log(`主方向差(边直方图): ${dd.toFixed(1)}°`);
console.log(`处理路径判定(旧判据 bbox<60): ${bboxArea(ir) < 60 ? 'SmallBuilding直通道' : '完整规则化'}; (新判据 面积<60): ${area(ir) < 60 ? 'SmallBuilding直通道' : '完整规则化'}`);

// 边序列
console.log('\nresult 顶点(相对质心):');
console.log('  ' + orr.map(p => [(p[0]-c2[0]).toFixed(1), (p[1]-c2[1]).toFixed(1)]).join(' '));
console.log('initial 前20顶点(相对质心):');
console.log('  ' + ir.slice(0, 20).map(p => [(p[0]-c1[0]).toFixed(1), (p[1]-c1[1]).toFixed(1)]).join(' '));
// 每条边方向
const edgeDirs = r => r.map((p, i) => { const q = r[(i+1)%r.length]; const l = Math.hypot(q[0]-p[0], q[1]-p[1]); return l > 0.3 ? `${(((Math.atan2(q[1]-p[1],q[0]-p[0])*180/Math.PI)%90+90)%90).toFixed(0)}°/${l.toFixed(1)}m` : null; }).filter(Boolean);
console.log('\nresult 边方向/长度: ' + edgeDirs(orr).join(' '));
