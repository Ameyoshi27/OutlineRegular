// 临时分析: ASCII 渲染未匹配自动要素 vs 附近人工要素
import { readFileSync } from 'fs';

function parseShp(path){const buf=readFileSync(path);const feats=[];let off=100;while(off+8<=buf.length){const cl=buf.readInt32BE(off+4)*2;const b=buf.subarray(off+8,off+8+cl);if(cl<=0||off+8+cl>buf.length)break;if(b.readInt32LE(0)===5||b.readInt32LE(0)===15||b.readInt32LE(0)===25){let p=36;const np=b.readInt32LE(p);p+=4;const npts=b.readInt32LE(p);p+=4;const parts=[];for(let i=0;i<np;i++){parts.push(b.readInt32LE(p));p+=4;}parts.push(npts);const rings=[];for(let i=0;i<np;i++){const ring=[];for(let j=parts[i];j<parts[i+1];j++){ring.push([b.readDoubleLE(p),b.readDoubleLE(p+8)]);p+=16;}if(ring.length>2)rings.push(ring);}feats.push({rings});}off+=8+cl;}return feats;}
const auto=parseShp('regularized_building.shp');
const human=parseShp('E:/LODData/tjh_shp_extract/shp/test.shp');
const cen=r=>{let x=0,y=0;for(const p of r){x+=p[0];y+=p[1];}return[x/r.length,y/r.length];};
const pip=(p,r)=>{let ins=false;for(let i=0,j=r.length-1;i<r.length;j=i++){const[xi,yi]=r[i],[xj,yj]=r[j];if((yi>p[1])!=(yj>p[1])&&p[0]<(xj-xi)*(p[1]-yi)/(yj-yi)+xi)ins=!ins;}return ins;};
const bboxOf=r=>{let a=1e18,b=1e18,c=-1e18,d=-1e18;for(const[x,y]of r){a=Math.min(a,x);b=Math.min(b,y);c=Math.max(c,x);d=Math.max(d,y);}return[a,b,c,d];};
const distToRing=(p,r)=>{let b=1e18;for(let i=0;i<r.length;i++){const[x1,y1]=r[i],[x2,y2]=r[(i+1)%r.length];const dx=x2-x1,dy=y2-y1;const t=Math.max(0,Math.min(1,((p[0]-x1)*dx+(p[1]-y1)*dy)/(dx*dx+dy*dy||1e-12)));b=Math.min(b,Math.hypot(p[0]-x1-t*dx,p[1]-y1-t*dy));}return b;};
const area=r=>Math.abs(r.reduce((s,p,i)=>{const q=r[(i+1)%r.length];return s+p[0]*q[1]-q[0]*p[1];},0))/2;

// 未匹配自动 = 质心不在任何人工内的(近似)
const unmatchedAuto=[];
for(const f of auto){
  const c=cen(f.rings[0]);
  if(!human.some(h=>h.rings.some(r=>pip(c,r)))) unmatchedAuto.push(f);
}
console.log('渲染前8个未匹配案例 (A=自动填充 a=自动边界, H=人工填充 h=人工边界, .=空)');
let shown=0;
for(const f of unmatchedAuto){
  if(shown>=8)break;
  const target=f.rings[0];
  const [tx1,ty1,tx2,ty2]=bboxOf(target);
  // 收集范围内所有要素
  const margin=15;
  const autos=auto.filter(g=>{const[b1,b2,b3,b4]=bboxOf(g.rings[0]);return b3>tx1-margin&&b1<tx2+margin&&b4>ty1-margin&&b2<ty2+margin;});
  const humans=human.filter(g=>g.rings.some(r=>{const[b1,b2,b3,b4]=bboxOf(r);return b3>tx1-margin&&b1<tx2+margin&&b4>ty1-margin&&b2<ty2+margin;}));
  const all=autos.flatMap(g=>g.rings).concat(humans.flatMap(g=>g.rings));
  let X1=1e18,Y1=1e18,X2=-1e18,Y2=-1e18;
  for(const r of all)for(const[x,y]of r){X1=Math.min(X1,x);Y1=Math.min(Y1,y);X2=Math.max(X2,x);Y2=Math.max(Y2,y);}
  const W=76, scale=Math.max((X2-X1)/W, (Y2-Y1)/40);
  const H=Math.ceil((Y2-Y1)/scale);
  const grid=Array.from({length:H+1},()=>Array(W+1).fill('.'));
  const mark=(x,y,ch)=>{const cx=Math.round((x-X1)/scale),cy=Math.round((y-Y1)/scale);if(cy>=0&&cy<=H&&cx>=0&&cx<=W){const cur=grid[cy][cx];if(cur==='.'||ch!=='.')grid[cy][cx]=ch;}};
  // 稀疏填充内部
  for(const g of autos)for(const r of g.rings)for(let x=X1;x<=X2;x+=scale)for(let y=Y1;y<=Y2;y+=scale/2){const p=[x+scale/4,y];const inpip=r.some?pip(p,r):false;if(pip(p,r))mark(x,y,'A');}
  for(const g of humans)for(const r of g.rings)for(let x=X1;x<=X2;x+=scale)for(let y=Y1;y<=Y2;y+=scale/2){const p=[x+scale/4,y];if(pip(p,r))mark(x,y,'H');}
  for(const g of autos)for(const r of g.rings)for(const[x,y]of r)mark(x,y,'a');
  for(const g of humans)for(const r of g.rings)for(const[x,y]of r)mark(x,y,'h');
  console.log(`\n案例${shown+1}: 自动面积=${area(target).toFixed(0)} 顶点=${target.length} 比例尺≈${scale.toFixed(1)}m/格 范围${((X2-X1)).toFixed(0)}x${((Y2-Y1)).toFixed(0)}m 自动要素数=${autos.length} 人工要素数=${humans.length}`);
  for(let row=H;row>=0;row--)console.log(grid[row].join(''));
  shown++;
}
