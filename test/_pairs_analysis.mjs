// 临时分析: 指定mask对 + 全量相交对分析
import { readFileSync } from 'fs';

function parseShpRecords(path) {
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
      const rings = [];
      for (let i = 0; i < np; i++) {
        const ring = [];
        for (let j = parts[i]; j < parts[i + 1]; j++) { ring.push([b.readDoubleLE(p), b.readDoubleLE(p + 8)]); p += 16; }
        if (ring.length > 1) {
          const a = ring[0], b2 = ring[ring.length - 1];
          if (Math.abs(a[0]-b2[0]) < 1e-9 && Math.abs(a[1]-b2[1]) < 1e-9) ring.pop();
        }
        if (ring.length >= 3) rings.push(ring);
      }
      recs.push({ rings, deleted: false });
    } else recs.push({ rings: null, deleted: true });
    off += 8 + cl;
  }
  return recs;
}
function parseDbf(path) {
  const buf = readFileSync(path);
  const n = buf.readUInt16LE(4), hs = buf.readUInt16LE(8), rs = buf.readUInt16LE(10);
  const fields = []; let p = 32;
  while (buf[p] !== 0x0D) { fields.push({ name: buf.toString('ascii', p, p+11).replace(/\0.*$/,'').trim(), type: String.fromCharCode(buf[p+11]), len: buf[p+16] }); p += 32; }
  const recs = [];
  for (let r = 0; r < n; r++) {
    const base = hs + r * rs; if (base + rs > buf.length) break;
    const rec = { deleted: buf[base] === 0x2A }; let q = base + 1;
    for (const f of fields) { const raw = buf.toString('latin1', q, q+f.len).trim(); rec[f.name] = f.type==='N'?parseFloat(raw):raw; q += f.len; }
    recs.push(rec);
  }
  return recs;
}

const resRecs = parseShpRecords('D:/outlineRegular/outlineRegular/test/regularized_building.shp');
const resDbf = parseDbf('D:/outlineRegular/outlineRegular/test/regularized_building.dbf');
const res = [];
for (let i = 0; i < resDbf.length; i++) if (!resDbf[i].deleted && resRecs[i]?.rings) res.push({ rid: i, mask: resDbf[i].mask, ring: resRecs[i].rings[0] });
console.log(`结果要素: ${res.length}`);

const area = r => Math.abs(r.reduce((s,p,i)=>{const q=r[(i+1)%r.length];return s+p[0]*q[1]-q[0]*p[1];},0))/2;
const bbox = r => { let a=1e18,b=1e18,c=-1e18,d=-1e18; for(const[x,y]of r){a=Math.min(a,x);b=Math.min(b,y);c=Math.max(c,x);d=Math.max(d,y);} return[a,b,c,d]; };
const pip = (p, r) => { let ins=false; for(let i=0,j=r.length-1;i<r.length;j=i++){const[xi,yi]=r[i],[xj,yj]=r[j];if((yi>p[1])!=(yj>p[1])&&p[0]<(xj-xi)*(p[1]-yi)/(yj-yi)+xi)ins=!ins;} return ins; };
const segInt=(p1,p2,p3,p4)=>{const d=(a,b,c)=>(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);const d1=d(p3,p4,p1),d2=d(p3,p4,p2),d3=d(p1,p2,p3),d4=d(p1,p2,p4);return((d1>0&&d2<0)||(d1<0&&d2>0))&&((d3>0&&d4<0)||(d3<0&&d4>0));};
function overlapArea(r1, r2, step=0.25) {
  const b1=bbox(r1), b2=bbox(r2);
  const minX=Math.max(b1[0],b2[0]),maxX=Math.min(b1[2],b2[2]),minY=Math.max(b1[1],b2[1]),maxY=Math.min(b1[3],b2[3]);
  if(minX>=maxX||minY>=maxY)return 0;
  let both=0,eith=0;
  for(let x=minX+step/2;x<maxX;x+=step)for(let y=minY+step/2;y<maxY;y+=step){const a=pip([x,y],r1),b=pip([x,y],r2);if(a&&b)both++;if(a||b)eith++;}
  return both*Math.pow(step,2);
}
function edgeAngles(r) {
  const bins = new Array(90).fill(0);
  for (let i = 0; i < r.length; i++) {
    const [x1,y1]=r[i],[x2,y2]=r[(i+1)%r.length];
    const len=Math.hypot(x2-x1,y2-y1); if(len<0.3)continue;
    bins[Math.floor(((Math.atan2(y2-y1,x2-x1)*180/Math.PI)%90+90)%90)]+=len;
  }
  let best=-1,bi=0;bins.forEach((v,i)=>{if(v>best){best=v;bi=i;}});
  return bi;
}

// 全量相交检测
const pairs = [];
for (let i = 0; i < res.length; i++) {
  const bi = bbox(res[i].ring);
  for (let j = i+1; j < res.length; j++) {
    const bj = bbox(res[j].ring);
    if (bi[2]<bj[0]||bj[2]<bi[0]||bi[3]<bj[1]||bj[3]<bi[1]) continue;
    let cross = false;
    outer: for (let a = 0; a < res[i].ring.length; a++) for (let b = 0; b < res[j].ring.length; b++)
      if (segInt(res[i].ring[a],res[i].ring[(a+1)%res[i].ring.length],res[j].ring[b],res[j].ring[(b+1)%res[j].ring.length])){cross=true;break outer;}
    if (!cross && (pip(res[i].ring[0],res[j].ring)||pip(res[j].ring[0],res[i].ring))) cross = true;
    if (cross) pairs.push([i,j]);
  }
}
console.log(`\n全量边界相交/包含对: ${pairs.length}`);
for (const [i,j] of pairs) {
  const a=res[i], b=res[j];
  const oa=overlapArea(a.ring,b.ring);
  const rect = a.ring.length===4 ? '矩形' : `${a.ring.length}顶点`;
  const rect2 = b.ring.length===4 ? '矩形' : `${b.ring.length}顶点`;
  console.log(`  mask=${a.mask}(${rect},${area(a.ring).toFixed(0)}m²,${edgeAngles(a.ring)}°) x mask=${b.ring}(${rect2},${area(b.ring).toFixed(0)}m²,${edgeAngles(b.ring)}°) 重叠≈${oa.toFixed(1)}m²`);
}
// 用户点名的mask对详情
console.log(`\n点名mask对详情:`);
const targets = [[4523,4527],[2988,4426],[4312,4133],[4333,4333]];
for (const [m1,m2] of targets) {
  const as = res.filter(r=>r.mask===m1), bs = res.filter(r=>r.mask===m2);
  for (const a of as) for (const b of bs) {
    if (a===b) continue;
    const b1=bbox(a.ring),b2=bbox(b.ring);
    if (b1[2]<b2[0]||b2[2]<b1[0]||b1[3]<b2[1]||b2[3]<b1[1]) continue;
    const oa = overlapArea(a.ring,b.ring);
    if (oa > 0.1) console.log(`  mask${m1}(v=${a.ring.length},a=${area(a.ring).toFixed(0)},${edgeAngles(a.ring)}°) x mask${m2}(v=${b.ring.length},a=${area(b.ring).toFixed(0)},${edgeAngles(b.ring)}°) 重叠≈${oa.toFixed(1)}m²`);
  }
}
