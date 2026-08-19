// 临时审计2: 按记录号对齐日志fid,检查自交/重叠面积量级/飞点要素与修复日志的关联
import { readFileSync } from 'fs';

function parseShpRecords(path) {
  const buf = readFileSync(path);
  const feats = [];
  let off = 100;
  while (off + 8 <= buf.length) {
    const contentLen = buf.readInt32LE ? buf.readInt32BE(off + 4) * 2 : 0;
    const body = buf.subarray(off + 8, off + 8 + contentLen);
    if (contentLen <= 0 || off + 8 + contentLen > buf.length) break;
    if (body.readInt32LE(0) === 5 || body.readInt32LE(0) === 15 || body.readInt32LE(0) === 25) {
      let p = 36;
      const numParts = body.readInt32LE(p); p += 4;
      const numPoints = body.readInt32LE(p); p += 4;
      const parts = [];
      for (let i = 0; i < numParts; i++) { parts.push(body.readInt32LE(p)); p += 4; }
      parts.push(numPoints);
      const rings = [];
      for (let i = 0; i < numParts; i++) {
        const ring = [];
        for (let j = parts[i]; j < parts[i + 1]; j++) { ring.push([body.readDoubleLE(p), body.readDoubleLE(p + 8)]); p += 16; }
        if (ring.length > 1) {
          const a = ring[0], b = ring[ring.length - 1];
          if (Math.abs(a[0]-b[0]) < 1e-9 && Math.abs(a[1]-b[1]) < 1e-9) ring.pop();
        }
        if (ring.length >= 3) rings.push(ring);
      }
      if (rings.length) feats.push({ fid: feats.length, rings });
    }
    off += 8 + contentLen;
  }
  return feats;
}

const resPath = process.argv[2], initPath = process.argv[3], logPath = process.argv[4];
const res = parseShpRecords(resPath);
const init = parseShpRecords(initPath);

const ringArea = r => Math.abs(r.reduce((s, p, i) => { const q = r[(i + 1) % r.length]; return s + p[0] * q[1] - q[0] * p[1]; }, 0)) / 2;
const bbox = r => { let a=1e18,b=1e18,c=-1e18,d=-1e18; for (const [x,y] of r){a=Math.min(a,x);b=Math.min(b,y);c=Math.max(c,x);d=Math.max(d,y);} return {minX:a,minY:b,maxX:c,maxY:d}; };
const pointInRing = (p, r) => { let ins=false; for (let i=0,j=r.length-1;i<r.length;j=i++){const [xi,yi]=r[i],[xj,yj]=r[j]; if((yi>p[1])!=(yj>p[1]) && p[0]<(xj-xi)*(p[1]-yi)/(yj-yi)+xi) ins=!ins;} return ins; };
const segInt=(p1,p2,p3,p4)=>{const d=(a,b,c)=>(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);const d1=d(p3,p4,p1),d2=d(p3,p4,p2),d3=d(p1,p2,p3),d4=d(p1,p2,p4);return ((d1>0&&d2<0)||(d1<0&&d2>0))&&((d3>0&&d4<0)||(d3<0&&d4>0));};

// 自交检测(每环)
function isSimple(r) {
  const n = r.length;
  for (let i = 0; i < n; i++) for (let j = i + 2; j < n; j++) {
    if (i === 0 && j === n - 1) continue;
    if (segInt(r[i], r[(i+1)%n], r[j], r[(j+1)%n])) return false;
  }
  return true;
}
const nonSimple = res.filter(f => !isSimple(f.rings[0]));
console.log(`自交环要素: ${nonSimple.length}/${res.length}: fids=` + nonSimple.map(f => f.fid).slice(0, 30).join(','));

// 相交对 + 网格采样估算重叠面积
function overlapArea(r1, r2) {
  const b1 = bbox(r1), b2 = bbox(r2);
  const minX = Math.max(b1.minX, b2.minX), maxX = Math.min(b1.maxX, b2.maxX);
  const minY = Math.max(b1.minY, b2.minY), maxY = Math.min(b1.maxY, b2.maxY);
  if (minX >= maxX || minY >= maxY) return 0;
  const step = 0.25;
  let cnt = 0, tot = 0;
  for (let x = minX + step / 2; x < maxX; x += step) for (let y = minY + step / 2; y < maxY; y += step) {
    tot++;
    const p = [x, y];
    if (pointInRing(p, r1) && pointInRing(p, r2)) cnt++;
  }
  return cnt / Math.max(tot, 1) * (maxX - minX) * (maxY - minY);
}
let pairs = [];
for (let i = 0; i < res.length; i++) {
  const bi = bbox(res[i].rings[0]);
  for (let j = i + 1; j < res.length; j++) {
    const bj = bbox(res[j].rings[0]);
    if (bi.maxX<bj.minX||bj.maxX<bi.minX||bi.maxY<bj.minY||bj.maxY<bi.minY) continue;
    let cross = false;
    outer: for (let a = 0; a < res[i].rings[0].length; a++) for (let b = 0; b < res[j].rings[0].length; b++)
      if (segInt(res[i].rings[0][a], res[i].rings[0][(a+1)%res[i].rings[0].length], res[j].rings[0][b], res[j].rings[0][(b+1)%res[j].rings[0].length])) { cross = true; break outer; }
    if (!cross && (pointInRing(res[i].rings[0][0], res[j].rings[0]) || pointInRing(res[j].rings[0][0], res[i].rings[0]))) cross = true;
    if (cross) pairs.push([i, j]);
  }
}
console.log(`相交对=${pairs.length}, 重叠面积估算(前15):`);
for (const [i, j] of pairs.slice(0, 15)) {
  const oa = overlapArea(res[i].rings[0], res[j].rings[0]);
  console.log(`  fid${res[i].fid}(a=${ringArea(res[i].rings[0]).toFixed(0)}) x fid${res[j].fid}(a=${ringArea(res[j].rings[0]).toFixed(0)}) ~${oa.toFixed(1)}m2`);
}

// 日志关联:修复涉及的fid
const log = readFileSync(logPath, 'latin1');
const diffFids = [...log.matchAll(/Difference removed fid=(\d+)/g)].map(m => +m[1]);
const dropFids = [...log.matchAll(/drop clipped fid=(\d+)/g)].map(m => +m[1]);
console.log(`Difference removed fids: ${diffFids.join(',')}`);
console.log(`dropped clipped fids: ${dropFids.join(',')}`);
// 飞点要素的自交状态
const distToRing=(p,r)=>{let b=1e18;for(let i=0;i<r.length;i++){const[x1,y1]=r[i],[x2,y2]=r[(i+1)%r.length];const dx=x2-x1,dy=y2-y1;const t=Math.max(0,Math.min(1,((p[0]-x1)*dx+(p[1]-y1)*dy)/(dx*dx+dy*dy||1e-12)));b=Math.min(b,Math.hypot(p[0]-x1-t*dx,p[1]-y1-t*dy));}return b;};
const initBoxes = init.map(f => bbox(f.rings[0]));
let fly = [];
for (const f of res) {
  let maxOut = 0;
  for (const p of f.rings[0]) {
    let best = 1e18;
    for (let m = 0; m < init.length; m++) {
      const b = initBoxes[m];
      if (p[0]<b.minX-60||p[0]>b.maxX+60||p[1]<b.minY-60||p[1]>b.maxY+60) continue;
      if (pointInRing(p, init[m].rings[0])) { best = 0; break; }
      best = Math.min(best, distToRing(p, init[m].rings[0]));
    }
    maxOut = Math.max(maxOut, best);
  }
  if (maxOut > 2) fly.push({ fid: f.fid, maxOut, simple: isSimple(f.rings[0]), verts: f.rings[0].length, area: ringArea(f.rings[0]).toFixed(0) });
}
console.log(`飞点要素(>2m): ${fly.length}`);
for (const f of fly.slice(0, 20)) console.log(`  fid${f.fid} maxOut=${f.maxOut.toFixed(1)}m selfIntersect=${!f.simple} verts=${f.verts} area=${f.area}`);
