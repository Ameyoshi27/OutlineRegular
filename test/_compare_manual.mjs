// 临时分析: 算法结果 vs 人工结果 系统对比
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
      if (rings.length) feats.push({ rings, ring: rings[0] });
    }
    off += 8 + contentLen;
  }
  return feats;
}

const ringArea = r => Math.abs(r.reduce((s, p, i) => { const q = r[(i + 1) % r.length]; return s + p[0] * q[1] - q[0] * p[1]; }, 0)) / 2;
const bbox = r => { let a=1e18,b=1e18,c=-1e18,d=-1e18; for (const [x,y] of r){a=Math.min(a,x);b=Math.min(b,y);c=Math.max(c,x);d=Math.max(d,y);} return {minX:a,minY:b,maxX:c,maxY:d}; };
const centroid = r => { let x=0,y=0; for (const p of r){x+=p[0];y+=p[1];} return [x/r.length,y/r.length]; };
const pointInRing = (p, r) => { let ins=false; for (let i=0,j=r.length-1;i<r.length;j=i++){const[xi,yi]=r[i],[xj,yj]=r[j]; if((yi>p[1])!=(yj>p[1]) && p[0]<(xj-xi)*(p[1]-yi)/(yj-yi)+xi) ins=!ins;} return ins; };

// 主方向: 边长加权直方图 mod 90
function domAngle(r) {
  const bins = new Array(90).fill(0);
  let total = 0;
  for (let i = 0; i < r.length; i++) {
    const [x1,y1]=r[i],[x2,y2]=r[(i+1)%r.length];
    const len = Math.hypot(x2-x1,y2-y1);
    if (len < 0.3) continue;
    let a = ((Math.atan2(y2-y1,x2-x1)*180/Math.PI) % 90 + 90) % 90;
    bins[Math.floor(a)] += len; total += len;
  }
  let best=-1, bi=0;
  bins.forEach((v,i)=>{if(v>best){best=v;bi=i;}});
  return { angle: bi, conf: total>0?best/total:0 };
}

// 网格 IoU
function iou(r1, r2, step=0.5) {
  const b1=bbox(r1), b2=bbox(r2);
  const minX=Math.max(b1.minX,b2.minX), maxX=Math.min(b1.maxX,b2.maxX);
  const minY=Math.max(b1.minY,b2.minY), maxY=Math.min(b1.maxY,b2.maxY);
  if (minX>=maxX||minY>=maxY) return 0;
  let both=0, either=0;
  for(let x=minX+step/2;x<maxX;x+=step) for(let y=minY+step/2;y<maxY;y+=step){
    const in1=pointInRing([x,y],r1), in2=pointInRing([x,y],r2);
    if(in1&&in2)both++; if(in1||in2)either++;
  }
  return either?both/either:0;
}

// 边界距离采样(点到另一环的最短距离, 环上采样)
function boundaryStats(rFrom, rTo, maxSamples=200) {
  const pts=[];
  const n=rFrom.length;
  const step=Math.max(1,Math.floor(n/maxSamples));
  for(let i=0;i<n;i+=step){
    const a=rFrom[i], b=rFrom[(i+1)%n];
    pts.push(a, [(a[0]+b[0])/2,(a[1]+b[1])/2]);
  }
  const ds=pts.map(p=>{
    let best=1e18;
    for(let i=0;i<rTo.length;i++){
      const [x1,y1]=rTo[i],[x2,y2]=rTo[(i+1)%rTo.length];
      const dx=x2-x1,dy=y2-y1;
      const t=Math.max(0,Math.min(1,((p[0]-x1)*dx+(p[1]-y1)*dy)/(dx*dx+dy*dy||1e-12)));
      best=Math.min(best,Math.hypot(p[0]-x1-t*dx,p[1]-y1-t*dy));
    }
    return best;
  }).sort((a,b)=>a-b);
  const mean=ds.reduce((s,v)=>s+v,0)/ds.length;
  const q90=ds[Math.floor(ds.length*0.9)];
  return {mean,q90};
}

// 凹凸性: 面积/凸包面积
function convexity(r) {
  const pts=[...r].sort((a,b)=>a[0]-b[0]||a[1]-b[1]);
  const cross=(o,a,b)=>(a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0]);
  const lower=[]; for(const p of pts){while(lower.length>=2&&cross(lower[lower.length-2],lower[lower.length-1],p)<=0)lower.pop();lower.push(p);}
  const upper=[]; for(let i=pts.length-1;i>=0;i--){const p=pts[i];while(upper.length>=2&&cross(upper[upper.length-2],upper[upper.length-1],p)<=0)upper.pop();upper.push(p);}
  const hull=lower.length+upper.length-2;
  let twice=0;
  for(let i=1;i<lower.length-1;i++)twice+=Math.abs(cross(lower[0],lower[i],lower[i+1]));
  // 凸包面积(鞋带)
  const hp=[...lower.slice(0,-1),...upper.slice(0,-1)];
  let ha=0; for(let i=0;i<hp.length;i++){const a=hp[i],b=hp[(i+1)%hp.length];ha+=a[0]*b[1]-b[0]*a[1];}
  return ringArea(r)/Math.abs(ha/2);
}

const auto = parseShp(process.argv[2]).map(f=>({ring:f.ring}));
const human = parseShp(process.argv[3]).map(f=>({ring:f.ring}));
console.log(`自动结果: ${auto.length} 要素, 人工结果: ${human.length} 要素`);

for (const [name, arr] of [['自动',auto],['人工',human]]) {
  const areas=arr.map(f=>ringArea(f.ring)).sort((a,b)=>a-b);
  const verts=arr.map(f=>f.ring.length);
  const avgV=verts.reduce((s,v)=>s+v,0)/verts.length;
  const med=areas[Math.floor(areas.length/2)];
  console.log(`${name}: 面积 min=${areas[0].toFixed(0)} 中位=${med.toFixed(0)} max=${areas[areas.length-1].toFixed(0)} 总和=${areas.reduce((s,v)=>s+v,0).toFixed(0)} | 顶点均值=${avgV.toFixed(1)} 中位=${verts.sort((a,b)=>a-b)[Math.floor(verts.length/2)]} | <60m²:${areas.filter(a=>a<60).length} <20m²:${areas.filter(a=>a<20).length}`);
  const r=arr[0].ring;
  console.log(`  坐标范围: x[${r[0][0].toFixed(0)}..] (首要素)`);
}

// bbox 预筛 + IoU 矩阵
const pairs=[];
for(let i=0;i<auto.length;i++){
  const bi=bbox(auto[i].ring);
  for(let j=0;j<human.length;j++){
    const bh=bbox(human[j].ring);
    if(bi.maxX<bh.minX||bh.maxX<bi.minX||bi.maxY<bh.minY||bh.maxY<bi.minY)continue;
    const v=iou(auto[i].ring,human[j].ring);
    if(v>0.05)pairs.push({i,j,v});
  }
}
// 贪心互最优匹配
pairs.sort((a,b)=>b.v-a.v);
const usedA=new Set(),usedH=new Set(),matched=[];
for(const p of pairs){
  if(usedA.has(p.i)||usedH.has(p.j))continue;
  usedA.add(p.i);usedH.add(p.j);matched.push(p);
}
console.log(`\n匹配对(IoU>0.05贪心): ${matched.length}`);
const unmatchedA=auto.map((_,i)=>i).filter(i=>!usedA.has(i));
const unmatchedH=human.map((_,j)=>j).filter(j=>!usedH.has(j));

// 匹配质量分布
const ious=matched.map(m=>m.v).sort((a,b)=>a-b);
const q=p=>ious[Math.floor(ious.length*p)];
console.log(`IoU分布: min=${ious[0].toFixed(2)} p25=${q(0.25).toFixed(2)} 中位=${q(0.5).toFixed(2)} p75=${q(0.75).toFixed(2)} max=${ious[ious.length-1].toFixed(2)}`);
console.log(`IoU>=0.8: ${ious.filter(v=>v>=0.8).length}, 0.6~0.8: ${ious.filter(v=>v>=0.6&&v<0.8).length}, 0.4~0.6: ${ious.filter(v=>v>=0.4&&v<0.6).length}, <0.4: ${ious.filter(v=>v<0.4).length}`);

// 差匹配案例细节(IoU<0.5)
console.log(`\n--- IoU<0.5 的匹配对(问题样本) ---`);
const problems=matched.filter(m=>m.v<0.5).sort((a,b)=>a.v-b.v);
for(const m of problems.slice(0,15)){
  const a=auto[m.i].ring,h=human[m.j].ring;
  const aa=ringArea(a),ha=ringArea(h);
  const da=domAngle(a),dh=domAngle(h);
  let dd=Math.abs(da.angle-dh.angle); if(dd>45)dd=90-dd;
  const bs=boundaryStats(a,h);
  console.log(`IoU=${m.v.toFixed(2)} auto(a=${aa.toFixed(0)},v=${a.length},ang=${da.angle}°conf=${(da.conf*100).toFixed(0)}%) human(a=${ha.toFixed(0)},v=${h.length},ang=${dh.angle}°conf=${(dh.conf*100).toFixed(0)}%) 面积比=${(aa/ha).toFixed(2)} 角差=${dd.toFixed(0)}° 边距mean=${bs.mean.toFixed(1)} q90=${bs.q90.toFixed(1)}`);
}

// 全体匹配的面积比/角差/边距统计
let areaRatios=[],angleDiffs=[],bmeans=[],bq90s=[],vdiffs=[];
for(const m of matched){
  const a=auto[m.i].ring,h=human[m.j].ring;
  areaRatios.push(ringArea(a)/ringArea(h));
  const da=domAngle(a),dh=domAngle(h);
  let dd=Math.abs(da.angle-dh.angle); if(dd>45)dd=90-dd;
  angleDiffs.push({d:dd,confA:da.conf,confH:dh.conf});
  const bs=boundaryStats(a,h);
  bmeans.push(bs.mean);bq90s.push(bs.q90);
  vdiffs.push(h.length-a.length);
}
const sorted=arr=>[...arr].sort((a,b)=>a-b);
const [ar1,ar2,ar3]=[sorted(areaRatios)[Math.floor(areaRatios.length*0.25)],sorted(areaRatios)[Math.floor(areaRatios.length*0.5)],sorted(areaRatios)[Math.floor(areaRatios.length*0.75)]];
console.log(`\n全体匹配: 面积比(auto/human) p25=${ar1.toFixed(2)} 中位=${ar2.toFixed(2)} p75=${ar3.toFixed(2)}`);
const [b1,b2,b3]=[sorted(bmeans)[Math.floor(bmeans.length*0.25)],sorted(bmeans)[Math.floor(bmeans.length*0.5)],sorted(bmeans)[Math.floor(bmeans.length*0.75)]];
console.log(`边界平均距离 p25=${b1.toFixed(2)} 中位=${b2.toFixed(2)} p75=${b3.toFixed(2)} m`);
const [q1,q2,q3]=[sorted(bq90s)[Math.floor(bq90s.length*0.25)],sorted(bq90s)[Math.floor(bq90s.length*0.5)],sorted(bq90s)[Math.floor(bq90s.length*0.75)]];
console.log(`边界q90距离 p25=${q1.toFixed(2)} 中位=${q2.toFixed(2)} p75=${q3.toFixed(2)} m`);
const ang5=angleDiffs.filter(x=>x.d>5&&x.confA>0.5&&x.confH>0.5).length;
const ang15=angleDiffs.filter(x=>x.d>15&&x.confA>0.5&&x.confH>0.5).length;
console.log(`主方向差>5°(双方置信>50%): ${ang15>0?ang5+'个(其中>15°:'+ang15+'个)':'无'}`);
console.log(`人工顶点数-自动顶点数: 中位=${sorted(vdiffs)[Math.floor(vdiffs.length/2)]} 均值=${(vdiffs.reduce((s,v)=>s+v,0)/vdiffs.length).toFixed(1)}`);

// 凹凸复杂度
const convA=auto.map(f=>convexity(f.ring));
const convH=human.map(f=>convexity(f.ring));
console.log(`凹凸度(1=凸): auto中位=${sorted(convA)[Math.floor(convA.length/2)].toFixed(2)} human中位=${sorted(convH)[Math.floor(convH.length/2)].toFixed(2)}`);

// 多对一/一对多: 拆分合并模式
const autoToHuman=new Map(), humanToAuto=new Map();
for(const p of pairs){
  if(p.v<0.1)continue;
  if(!autoToHuman.has(p.i))autoToHuman.set(p.i,[]);
  autoToHuman.get(p.i).push(p.j);
  if(!humanToAuto.has(p.j))humanToAuto.set(p.j,[]);
  humanToAuto.get(p.j).push(p.i);
}
let overSplit=0,underSplit=0;
for(const [j,is] of humanToAuto) if(is.length>=2) overSplit++;
for(const [i,js] of autoToHuman) if(js.length>=2) underSplit++;
console.log(`\n一个人工要素对应多个自动(IoU>0.1): ${overSplit} 个 (自动过拆或人工合并)`);
console.log(`一个自动要素对应多个人工(IoU>0.1): ${underSplit} 个 (自动漏拆或人工过拆)`);

// 未匹配的自动要素(多余)
console.log(`\n--- 未匹配自动要素 (${unmatchedA.length}个) 面积分布 ---`);
const ua=unmatchedA.map(i=>ringArea(auto[i].ring)).sort((a,b)=>a-b);
if(ua.length)console.log(`min=${ua[0].toFixed(0)} 中位=${ua[Math.floor(ua.length/2)].toFixed(0)} max=${ua[ua.length-1].toFixed(0)}; <20m²:${ua.filter(a=>a<20).length} <60:${ua.filter(a=>a<60).length} >=60:${ua.filter(a=>a>=60).length}`);
// 检查未匹配自动是否其实靠近某人工(近距离=小位移未匹配)
let nearMiss=0;
for(const i of unmatchedA){
  const c=centroid(auto[i].ring);
  let mind=1e18;
  for(const j of unmatchedH){const ch=centroid(human[j].ring);mind=Math.min(mind,Math.hypot(c[0]-ch[0],c[1]-ch[1]));}
  if(mind<10)nearMiss++;
}
console.log(`未匹配自动中质心距某未匹配人工<10m的: ${nearMiss}`);
// 未匹配的人工要素(漏检)
console.log(`\n--- 未匹配人工要素 (${unmatchedH.length}个) 面积分布 ---`);
const uh=unmatchedH.map(j=>ringArea(human[j].ring)).sort((a,b)=>a-b);
if(uh.length)console.log(`min=${uh[0].toFixed(0)} 中位=${uh[Math.floor(uh.length/2)].toFixed(0)} max=${uh[uh.length-1].toFixed(0)}; <20m²:${uh.filter(a=>a<20).length} 20~60:${uh.filter(a=>a>=20&&a<60).length} >=60:${uh.filter(a=>a>=60).length}`);
// 漏检的人工要素是否有自动要素在旁边(部分重叠但IoU低)
let partial=0, none=0;
for(const j of unmatchedH){
  const bh=bbox(human[j].ring);
  let hasNear=false;
  for(const p of pairs){ if(p.j===j&&p.v>0.02){hasNear=true;break;} }
  // bbox 邻近检查
  if(!hasNear)for(let i=0;i<auto.length;i++){const ba=bbox(auto[i].ring);
    if(ba.maxX<bh.minX-3||bh.maxX+3<ba.minX||ba.maxY<bh.minY-3||bh.maxY+3<ba.minY)continue;
    if(iou(auto[i].ring,human[j].ring,1.0)>0.02){hasNear=true;break;}}
  hasNear?partial++:none++;
}
console.log(`未匹配人工: 有微弱重叠自动要素(部分检出) ${partial}, 完全无自动要素(整体漏检) ${none}`);
// 整体漏检的面积
const missedAreas=[];
for(const j of unmatchedH){
  const bh=bbox(human[j].ring);let hasNear=false;
  for(let i=0;i<auto.length;i++){const ba=bbox(auto[i].ring);
    if(ba.maxX<bh.minX-3||bh.maxX+3<ba.minX||ba.maxY<bh.minY-3||bh.maxY+3<ba.minY)continue;
    if(iou(auto[i].ring,human[j].ring,1.0)>0.02){hasNear=true;break;}}
  if(!hasNear)missedAreas.push(ringArea(human[j].ring));
}
if(missedAreas.length){missedAreas.sort((a,b)=>a-b);console.log(`整体漏检面积: min=${missedAreas[0].toFixed(0)} 中位=${missedAreas[Math.floor(missedAreas.length/2)].toFixed(0)} max=${missedAreas[missedAreas.length-1].toFixed(0)}`);}
