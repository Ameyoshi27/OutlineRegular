// 临时分析脚本：量化四个问题(小多边形/相交/方向偏差/曲线段丢失)
import { readFileSync } from 'fs';

function parseShp(path) {
  const buf = readFileSync(path);
  const records = [];
  let off = 100;
  while (off + 8 <= buf.length) {
    const contentLen = buf.readInt32BE(off + 4) * 2;
    const body = buf.subarray(off + 8, off + 8 + contentLen);
    if (contentLen <= 0 || off + 8 + contentLen > buf.length) break;
    const type = body.readInt32LE(0);
    if (type === 5 || type === 15 || type === 25) {
      let p = 4;
      const box = [body.readDoubleLE(p), body.readDoubleLE(p+8), body.readDoubleLE(p+16), body.readDoubleLE(p+24)];
      p += 32;
      const numParts = body.readInt32LE(p); p += 4;
      const numPoints = body.readInt32LE(p); p += 4;
      const parts = [];
      for (let i = 0; i < numParts; i++) { parts.push(body.readInt32LE(p)); p += 4; }
      parts.push(numPoints);
      const rings = [];
      for (let i = 0; i < numParts; i++) {
        const ring = [];
        for (let j = parts[i]; j < parts[i + 1]; j++) {
          ring.push([body.readDoubleLE(p), body.readDoubleLE(p + 8)]); p += 16;
        }
        // 去掉闭合重复点
        if (ring.length > 1) {
          const a = ring[0], b = ring[ring.length - 1];
          if (Math.abs(a[0]-b[0]) < 1e-9 && Math.abs(a[1]-b[1]) < 1e-9) ring.pop();
        }
        rings.push(ring);
      }
      records.push({ num: buf.readInt32BE(off), box, rings });
    }
    off += 8 + contentLen;
  }
  return records;
}

const ringArea = r => Math.abs(r.reduce((s, p, i) => {
  const q = r[(i + 1) % r.length];
  return s + p[0] * q[1] - q[0] * p[1];
}, 0)) / 2;

const centroid = r => {
  let x = 0, y = 0;
  for (const p of r) { x += p[0]; y += p[1]; }
  return [x / r.length, y / r.length];
};

const bbox = r => {
  let minX = 1e18, minY = 1e18, maxX = -1e18, maxY = -1e18;
  for (const [x, y] of r) { minX = Math.min(minX,x); minY = Math.min(minY,y); maxX = Math.max(maxX,x); maxY = Math.max(maxY,y); }
  return { minX, minY, maxX, maxY, area: (maxX-minX)*(maxY-minY) };
};

// 主方向:按边长加权的方向直方图(mod 90°)
const dominantAngle = r => {
  const bins = new Array(90).fill(0);
  for (let i = 0; i < r.length; i++) {
    const [x1,y1] = r[i], [x2,y2] = r[(i+1)%r.length];
    const len = Math.hypot(x2-x1, y2-y1);
    if (len < 0.5) continue;
    let ang = Math.atan2(y2-y1, x2-x1) * 180 / Math.PI;
    ang = ((ang % 90) + 90) % 90;
    bins[Math.floor(ang)] += len;
  }
  let best = 0, bestLen = -1, total = 0;
  bins.forEach((l, i) => { total += l; if (l > bestLen) { bestLen = l; best = i; } });
  return { angle: best, ratio: total > 0 ? bestLen / total : 0 };
};

// 最小外接矩形(旋转卡壳简化:遍历每条边方向)
const minRect = r => {
  let best = null;
  for (let i = 0; i < r.length; i++) {
    const [x1,y1] = r[i], [x2,y2] = r[(i+1)%r.length];
    const len = Math.hypot(x2-x1, y2-y1);
    if (len < 0.3) continue;
    const ux = (x2-x1)/len, uy = (y2-y1)/len;
    let minU=1e18,maxU=-1e18,minV=1e18,maxV=-1e18;
    for (const [x,y] of r) {
      const u = x*ux+y*uy, v = -x*uy+y*ux;
      minU=Math.min(minU,u);maxU=Math.max(maxU,u);minV=Math.min(minV,v);maxV=Math.max(maxV,v);
    }
    const a = (maxU-minU)*(maxV-minV);
    if (!best || a < best.area) best = { area: a, w: maxU-minU, h: maxV-minV, angle: Math.atan2(uy,ux)*180/Math.PI };
  }
  return best;
};

function segIntersect(p1,p2,p3,p4) {
  const d = (a,b,c) => (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
  const d1=d(p3,p4,p1),d2=d(p3,p4,p2),d3=d(p1,p2,p3),d4=d(p1,p2,p4);
  return ((d1>0&&d2<0)||(d1<0&&d2>0)) && ((d3>0&&d4<0)||(d3<0&&d4>0));
}

function pointInRing(p, r) {
  let inside = false;
  for (let i = 0, j = r.length-1; i < r.length; j = i++) {
    const [xi,yi]=r[i],[xj,yj]=r[j];
    if ((yi>p[1])!=(yj>p[1]) && p[0] < (xj-xi)*(p[1]-yi)/(yj-yi)+xi) inside=!inside;
  }
  return inside;
}

function ringsCross(r1, r2) {
  for (let i=0;i<r1.length;i++) for (let j=0;j<r2.length;j++)
    if (segIntersect(r1[i],r1[(i+1)%r1.length],r2[j],r2[(j+1)%r2.length])) return true;
  return pointInRing(r1[0],r2) || pointInRing(r2[0],r1);
}

const fmt = (recs, name) => {
  const feats = recs.map(r => ({ ring: r.rings[0], area: ringArea(r.rings[0]) })).filter(f => f.ring && f.ring.length >= 3);
  console.log(`\n=== ${name}: ${feats.length} features ===`);
  const areas = feats.map(f => f.area).sort((a,b)=>a-b);
  const small = feats.filter(f => f.area < 60);
  console.log(`area min/median/max: ${areas[0]?.toFixed(1)}/${areas[Math.floor(areas.length/2)]?.toFixed(1)}/${areas[areas.length-1]?.toFixed(1)}; <60m²: ${small.length}`);
  for (const f of small.slice(0, 12)) {
    const b = bbox(f.ring);
    const mr = minRect(f.ring);
    console.log(`  small: area=${f.area.toFixed(1)} verts=${f.ring.length} bboxArea=${b.area.toFixed(1)} minRect=${mr.area.toFixed(1)} (${mr.w.toFixed(1)}x${mr.h.toFixed(1)})`);
  }
  return feats;
};

// ---- 问题1+2: 全量结果 ----
if (process.argv[2] && process.argv[2] !== 'skip') {
  const res = fmt(parseShp(process.argv[2]), process.argv[2]);
  let crossPairs = 0, checked = 0;
  const crosses = [];
  for (let i = 0; i < res.length; i++) {
    const bi = bbox(res[i].ring);
    for (let j = i+1; j < res.length; j++) {
      const bj = bbox(res[j].ring);
      if (bi.maxX<bj.minX||bj.maxX<bi.minX||bi.maxY<bj.minY||bj.maxY<bi.minY) continue;
      checked++;
      if (ringsCross(res[i].ring, res[j].ring)) { crossPairs++; if (crosses.length<15) crosses.push([i,j,res[i].area,res[j].area]); }
    }
  }
  console.log(`\n相交检测: bbox候选对=${checked}, 边界相交/包含对=${crossPairs}`);
  for (const [i,j,a1,a2] of crosses) console.log(`  overlap pair: feat${i}(area=${a1.toFixed(0)}) x feat${j}(area=${a2.toFixed(0)})`);
}

// ---- 问题3: 1_initial vs 1 方向对比 ----
if (process.argv[3] && process.argv[4]) {
  const init = parseShp(process.argv[3]).map(r=>r.rings[0]).filter(r=>r&&r.length>=3);
  const out = parseShp(process.argv[4]).map(r=>r.rings[0]).filter(r=>r&&r.length>=3);
  console.log(`\n=== 方向对比: initial=${init.length} result=${out.length} ===`);
  for (const o of out) {
    const c = centroid(o);
    let best = null, bd = 1e18;
    for (const s of init) { const sc = centroid(s); const d = Math.hypot(sc[0]-c[0],sc[1]-c[1]); if (d<bd){bd=d;best=s;} }
    if (!best || bd > 15) continue;
    const da = dominantAngle(o), db = dominantAngle(best);
    let delta = Math.abs(da.angle-db.angle); if (delta>45) delta=90-delta;
    const oa = ringArea(o), ob = ringArea(best);
    console.log(`feat: initArea=${ob.toFixed(0)} outArea=${oa.toFixed(0)} areaRatio=${(oa/ob).toFixed(2)} initVerts=${best.length} outVerts=${o.length} initAngle=${db.angle.toFixed(0)}(conf ${(db.ratio*100).toFixed(0)}%) outAngle=${da.angle.toFixed(0)}(conf ${(da.ratio*100).toFixed(0)}%) delta=${delta.toFixed(1)}° centroidShift=${bd.toFixed(1)}m`);
  }
}

// ---- 问题4: 2_initial vs 2 曲线段 ----
if (process.argv[5] && process.argv[6]) {
  const init = parseShp(process.argv[5]).map(r=>r.rings[0]).filter(r=>r&&r.length>=3);
  const out = parseShp(process.argv[6]).map(r=>r.rings[0]).filter(r=>r&&r.length>=3);
  console.log(`\n=== 曲线对比: initial=${init.length} result=${out.length} ===`);
  const chainStats = ring => {
    // 统计"方向连续变化的短边链":边长<2m且相邻边转角5°~35°的链
    const n = ring.length;
    let chainEdges = 0, maxChain = 0, cur = 0;
    for (let i = 0; i < n; i++) {
      const [x1,y1]=ring[i],[x2,y2]=ring[(i+1)%n],[x3,y3]=ring[(i+2)%n];
      const l1=Math.hypot(x2-x1,y2-y1), l2=Math.hypot(x3-x2,y3-y2);
      let a1=Math.atan2(y2-y1,x2-x1), a2=Math.atan2(y3-y2,x3-x2);
      let turn=Math.abs((a2-a1)*180/Math.PI); if(turn>180)turn=360-turn;
      if (l1<2.5 && l2<2.5 && turn>4 && turn<40) { cur++; chainEdges++; maxChain=Math.max(maxChain,cur); }
      else cur=0;
    }
    return { chainEdges, maxChain };
  };
  for (const o of out) {
    const c = centroid(o);
    let best=null,bd=1e18;
    for (const s of init){const sc=centroid(s);const d=Math.hypot(sc[0]-c[0],sc[1]-c[1]);if(d<bd){bd=d;best=s;}}
    if (!best||bd>15) continue;
    const si=chainStats(best), so=chainStats(o);
    console.log(`feat: initArea=${ringArea(best).toFixed(0)} outArea=${ringArea(o).toFixed(0)} initVerts=${best.length} outVerts=${o.length} initCurveEdges=${si.chainEdges}(maxChain=${si.maxChain}) outCurveEdges=${so.chainEdges} angleInit=${dominantAngle(best).angle.toFixed(0)} angleOut=${dominantAngle(o).angle.toFixed(0)}`);
  }
}
