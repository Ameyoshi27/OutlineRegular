// 临时分析: initial 1320 / hypothesis 116 / result 116 三阶段对比
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

const dir = 'D:/outlineRegular/outlineRegular/build_deps_release/Release/';
const initActive = parseShpRecords(dir + 'initial_building_outline.shp').filter(r => r.rings);
const hypAll = parseShpRecords(dir + 'debug_best_hypothesis.shp').filter(r => r.rings);
const resDbf = parseDbf('D:/outlineRegular/outlineRegular/test/regularized_building.dbf');
const resActive = parseShpRecords('D:/outlineRegular/outlineRegular/test/regularized_building.shp').filter(r => r.rings);
console.log(`active: initial=${initActive.length} hyp=${hypAll.length} result=${resActive.length}/${resDbf.length}`);

const ini = initActive[1320], hyp = hypAll[116], res = resActive[116];
const area = r => Math.abs(r.reduce((s,p,i)=>{const q=r[(i+1)%r.length];return s+p[0]*q[1]-q[0]*p[1];},0))/2;
const cen = r => { let x=0,y=0; for(const p of r){x+=p[0];y+=p[1];} return [x/r.length,y/r.length]; };
console.log(`initial[1320]: 顶点=${ini.rings[0].length} 面积=${area(ini.rings[0]).toFixed(1)}`);
console.log(`hyp[116]:      顶点=${hyp.rings[0].length} 面积=${area(hyp.rings[0]).toFixed(1)}`);
console.log(`result[116]:   顶点=${res.rings[0].length} 面积=${area(res.rings[0]).toFixed(1)} mask=${resDbf[116] ? resDbf[116].mask : '?'}`);
const c1=cen(ini.rings[0]), c2=cen(hyp.rings[0]), c3=cen(res.rings[0]);
console.log(`质心距: ini-hyp=${Math.hypot(c1[0]-c2[0],c1[1]-c2[1]).toFixed(1)} hyp-res=${Math.hypot(c2[0]-c3[0],c2[1]-c3[1]).toFixed(1)}`);

function describe(r, name) {
  const n = r.length;
  let turnSum = 0; const turns = [];
  for (let i = 0; i < n; i++) {
    const a = r[(i+n-1)%n], b = r[i], c = r[(i+1)%n];
    const cr = (b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0]);
    turns.push(cr); turnSum += cr;
  }
  const sgn = turnSum >= 0 ? 1 : -1;
  const edges = [];
  for (let i = 0; i < n; i++) {
    const p = r[i], q = r[(i+1)%n];
    const len = Math.hypot(q[0]-p[0], q[1]-p[1]);
    let ang = ((Math.atan2(q[1]-p[1], q[0]-p[0])*180/Math.PI)%90+90)%90;
    edges.push(`${len.toFixed(1)}@${ang.toFixed(0)}°${turns[i]*sgn<0?'[凹]':''}`);
  }
  let reflex = 0; for (const t of turns) if (t*sgn < 0) reflex++;
  console.log(`${name}: 凹顶点=${reflex}`);
  console.log('  ' + edges.join(' '));
}
describe(hyp.rings[0], 'hyp[116]');
describe(res.rings[0], 'result[116]');

// 初始轮廓的简化形状
function dp(points, tol) {
  if (points.length < 3) return points;
  const keep = new Array(points.length).fill(false);
  keep[0] = keep[points.length-1] = true;
  const stack = [[0, points.length-1]];
  while (stack.length) {
    const [a, b] = stack.pop();
    let maxD = -1, idx = -1;
    const ax=points[a][0],ay=points[a][1],bx=points[b][0],by=points[b][1];
    const dx=bx-ax,dy=by-ay,L=Math.hypot(dx,dy)||1e-12;
    for (let i=a+1;i<b;i++){const d=Math.abs((points[i][0]-ax)*dy-(points[i][1]-ay)*dx)/L;if(d>maxD){maxD=d;idx=i;}}
    if(maxD>tol){keep[idx]=true;stack.push([a,idx],[idx,b]);}
  }
  return points.filter((_,i)=>keep[i]);
}
const simp = dp(ini.rings[0], 0.45);
console.log(`initial DP(0.45): ${simp.length} 顶点`);
describe(simp, 'initial[1320]简化');

// 叠加ASCII
const all = [simp, hyp.rings[0], res.rings[0]];
let X1=1e18,Y1=1e18,X2=-1e18,Y2=-1e18;
for (const r of all) for (const [x,y] of r) { X1=Math.min(X1,x);Y1=Math.min(Y1,y);X2=Math.max(X2,x);Y2=Math.max(Y2,y); }
const W=76, scale=Math.max((X2-X1)/W,(Y2-Y1)/36), H=Math.ceil((Y2-Y1)/scale);
const g = Array.from({length:H+1},()=>Array(W+1).fill('.'));
const mark=(x,y,ch)=>{const cx=Math.round((x-X1)/scale),cy=Math.round((y-Y1)/scale);if(cy>=0&&cy<=H&&cx>=0&&cx<=W)g[cy][cx]=ch;};
for(const[x,y]of simp)mark(x,y,'i');
for(const[x,y]of hyp.rings[0])mark(x,y,'h');
for(const[x,y]of res.rings[0])mark(x,y,'r');
console.log(`叠加(i=初始 h=假设 r=结果, ${scale.toFixed(1)}m/格):`);
for(let row=H;row>=0;row--)console.log(g[row].join(''));
