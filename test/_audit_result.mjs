// 临时审计: 新结果 vs 初始轮廓 —— 相交对、飞点、伪曲线统计
import { readFileSync } from 'fs';

function parseShp(path) {
  const buf = readFileSync(path);
  const feats = [];
  let off = 100;
  while (off + 8 <= buf.length) {
    const contentLen = buf.readInt32BE(off + 4) * 2;
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
      if (rings.length) feats.push({ rings });
    }
    off += 8 + contentLen;
  }
  return feats;
}

const ringArea = r => Math.abs(r.reduce((s, p, i) => {
  const q = r[(i + 1) % r.length];
  return s + p[0] * q[1] - q[0] * p[1];
}, 0)) / 2;

const bbox = r => {
  let minX=1e18,minY=1e18,maxX=-1e18,maxY=-1e18;
  for (const [x,y] of r){minX=Math.min(minX,x);minY=Math.min(minY,y);maxX=Math.max(maxX,x);maxY=Math.max(maxY,y);}
  return {minX,minY,maxX,maxY};
};

const pointInRing = (p, r) => {
  let inside = false;
  for (let i = 0, j = r.length-1; i < r.length; j = i++) {
    const [xi,yi]=r[i],[xj,yj]=r[j];
    if ((yi>p[1])!=(yj>p[1]) && p[0] < (xj-xi)*(p[1]-yi)/(yj-yi)+xi) inside=!inside;
  }
  return inside;
};

const distToRing = (p, r) => {
  let best = 1e18;
  for (let i = 0; i < r.length; i++) {
    const [x1,y1]=r[i],[x2,y2]=r[(i+1)%r.length];
    const dx=x2-x1, dy=y2-y1;
    const t=Math.max(0,Math.min(1,((p[0]-x1)*dx+(p[1]-y1)*dy)/(dx*dx+dy*dy||1e-12)));
    best=Math.min(best,Math.hypot(p[0]-x1-t*dx,p[1]-y1-t*dy));
  }
  return best;
};

const segIntersect=(p1,p2,p3,p4)=>{
  const d=(a,b,c)=>(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
  const d1=d(p3,p4,p1),d2=d(p3,p4,p2),d3=d(p1,p2,p3),d4=d(p1,p2,p4);
  return ((d1>0&&d2<0)||(d1<0&&d2>0))&&((d3>0&&d4<0)||(d3<0&&d4>0));
};

const resPath = process.argv[2], initPath = process.argv[3];
const res = parseShp(resPath).map(f => f.rings[0]);
const init = parseShp(initPath).map(f => f.rings[0]);
console.log(`result feats=${res.length} initial feats=${init.length}`);

// 1) 顶点统计: 曲线化程度(总顶点数、平均每要素)
const totalVerts = res.reduce((s, r) => s + r.length, 0);
console.log(`result total verts=${totalVerts}, avg=${(totalVerts/res.length).toFixed(1)} (旧版约 8-12/要素)`);

// 2) 伪曲线启发: 输出环中"相邻边转角 3°~35° 且边长<1.5m"的边占比(正交矩形应≈0)
const curveEdges = r => {
  let c = 0, straight = 0;
  for (let i = 0; i < r.length; i++) {
    const [x1,y1]=r[i],[x2,y2]=r[(i+1)%r.length],[x3,y3]=r[(i+2)%r.length];
    const l=Math.hypot(x2-x1,y2-y1);
    let t=Math.abs((Math.atan2(y3-y2,x3-x2)-Math.atan2(y2-y1,x2-x1))*180/Math.PI);
    if(t>180)t=360-t;
    if (l<1.5 && t>2 && t<40) c++; else straight++;
  }
  return c;
};
const featsWithCurves = res.filter(r => curveEdges(r) >= 3).length;
console.log(`疑似含伪曲线段的要素(≥3条小转角短边): ${featsWithCurves}/${res.length}`);

// 3) 相交对检测
let crossPairs = [];
for (let i = 0; i < res.length; i++) {
  const bi = bbox(res[i]);
  for (let j = i+1; j < res.length; j++) {
    const bj = bbox(res[j]);
    if (bi.maxX<bj.minX||bj.maxX<bi.minX||bi.maxY<bj.minY||bj.maxY<bi.minY) continue;
    let cross = false;
    outer: for (let a = 0; a < res[i].length; a++) for (let b = 0; b < res[j].length; b++)
      if (segIntersect(res[i][a],res[i][(a+1)%res[i].length],res[j][b],res[j][(b+1)%res[j].length])) { cross = true; break outer; }
    if (!cross && (pointInRing(res[i][0],res[j]) || pointInRing(res[j][0],res[i]))) cross = true;
    if (cross) crossPairs.push([i,j]);
  }
}
console.log(`边界相交/包含对: ${crossPairs.length}`);
for (const [i,j] of crossPairs.slice(0,10)) console.log(`  pair feat${i}(area=${ringArea(res[i]).toFixed(0)}) x feat${j}(area=${ringArea(res[j]).toFixed(0)})`);

// 4) 飞点: 结果顶点在所有初始轮廓外且距离>2m
const initBoxes = init.map(bbox);
let flyFeatures = 0;
const flySamples = [];
for (let k = 0; k < res.length; k++) {
  let maxOut = 0;
  for (const p of res[k]) {
    let bestDist = 1e18;
    for (let m = 0; m < init.length; m++) {
      const b = initBoxes[m];
      if (p[0] < b.minX - 50 || p[0] > b.maxX + 50 || p[1] < b.minY - 50 || p[1] > b.maxY + 50) continue;
      if (pointInRing(p, init[m])) { bestDist = 0; break; }
      bestDist = Math.min(bestDist, distToRing(p, init[m]));
    }
    maxOut = Math.max(maxOut, bestDist);
  }
  if (maxOut > 2) { flyFeatures++; if (flySamples.length < 12) flySamples.push(`feat${k}(area=${ringArea(res[k]).toFixed(0)},verts=${res[k].length},maxOut=${maxOut.toFixed(1)}m)`); }
}
console.log(`飞出初始轮廓(>2m)的要素: ${flyFeatures}/${res.length}`);
for (const s of flySamples) console.log(`  ${s}`);
