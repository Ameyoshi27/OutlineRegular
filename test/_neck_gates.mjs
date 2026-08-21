// 临时分析: 在 id 1948/1949 轮廓上复现窄颈判据的每道门槛
import { readFileSync } from 'fs';

function parseShpRecords(path){const buf=readFileSync(path);const recs=[];let off=100;while(off+8<=buf.length){const cl=buf.readInt32BE(off+4)*2;const b=buf.subarray(off+8,off+8+cl);if(cl<=0||off+8+cl>buf.length)break;const t=b.readInt32LE(0);if(t===5||t===15||t===25){let p=36;const np=b.readInt32LE(p);p+=4;const npts=b.readInt32LE(p);p+=4;const parts=[];for(let i=0;i<np;i++){parts.push(b.readInt32LE(p));p+=4;}parts.push(npts);const rings=[];for(let i=0;i<np;i++){const ring=[];for(let j=parts[i];j<parts[i+1];j++){ring.push([b.readDoubleLE(p),b.readDoubleLE(p+8)]);p+=16;}if(ring.length>2)rings.push(ring);}recs.push({rings,deleted:false});}else recs.push({rings:null,deleted:true});off+=8+cl;}return recs;}
function parseDbf(path){const buf=readFileSync(path);const n=buf.readUInt16LE(4),hs=buf.readUInt16LE(8),rs=buf.readUInt16LE(10);const fields=[];let p=32;while(buf[p]!==0x0D){fields.push({name:buf.toString('ascii',p,p+11).replace(/\0.*$/,'').trim(),type:String.fromCharCode(buf[p+11]),len:buf[p+16]});p+=32;}const recs=[];for(let r=0;r<n;r++){const base=hs+r*rs;if(base+rs>buf.length)break;const rec={deleted:buf[base]===0x2A};let q=base+1;for(const f of fields){const raw=buf.toString('latin1',q,q+f.len).trim();rec[f.name]=f.type==='N'?parseFloat(raw):raw;q+=f.len;}recs.push(rec);}return recs;}
const dir='D:/outlineRegular/outlineRegular/build_deps_release/Release/';
const recs=parseShpRecords(dir+'initial_building_outline.shp');
const dbf=parseDbf(dir+'initial_building_outline.dbf');

const pip=(p,r)=>{let ins=false;for(let i=0,j=r.length-1;i<r.length;j=i++){const[xi,yi]=r[i],[xj,yj]=r[j];if((yi>p[1])!=(yj>p[1])&&p[0]<(xj-xi)*(p[1]-yi)/(yj-yi)+xi)ins=!ins;}return ins;};
const distToSeg=(p,a,b)=>{const dx=b[0]-a[0],dy=b[1]-a[1];const L2=dx*dx+dy*dy||1e-12;let t=((p[0]-a[0])*dx+(p[1]-a[1])*dy)/L2;t=Math.max(0,Math.min(1,t));return Math.hypot(p[0]-a[0]-t*dx,p[1]-a[1]-t*dy);};
// 切线在多边形内: 采样点全部在内部或极贴近边界
const cutInside=(r,i,j)=>{
  const a=r[i],b=r[j];
  const steps=Math.max(2,Math.ceil(Math.hypot(b[0]-a[0],b[1]-a[1])/0.05));
  for(let s=1;s<steps;s++){
    const t=s/steps;
    const p=[a[0]+t*(b[0]-a[0]),a[1]+t*(b[1]-a[1])];
    if(pip(p,r))continue;
    // 边界容差 0.05m
    let near=false;
    for(let k=0;k<r.length;k++){if(distToSeg(p,r[k],r[(k+1)%r.length])<0.05){near=true;break;}}
    if(!near)return false;
  }
  return true;
};
const area=r=>Math.abs(r.reduce((s,p,i)=>{const q=r[(i+1)%r.length];return s+p[0]*q[1]-q[0]*p[1];},0))/2;

for(const target of [1948,1949]){
  let ring=null;
  for(let r=0;r<dbf.length;r++){
    if(dbf[r].deleted||!recs[r]?.rings?.length)continue;
    if(dbf[r].id===target){ring=recs[r].rings[0];break;}
  }
  if(!ring){console.log(target,'未找到');continue;}
  const n=ring.length;
  const edge=[],prefix=[0];
  for(let i=0;i<n;i++){edge.push(Math.hypot(ring[(i+1)%n][0]-ring[i][0],ring[(i+1)%n][1]-ring[i][1]));prefix.push(prefix[i]+edge[i]);}
  const per=prefix[n];
  console.log(`\n=== id=${target}: 顶点=${n} 周长=${per.toFixed(0)}m 面积=${area(ring).toFixed(0)}m² ===`);
  console.log('宽度≤4.5m 且 两侧弧长≥4m 的顶点对及各门槛状态:');
  let shown=0;
  const pairs=[];
  for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
    const fw=j-i,bw=n-fw;
    if(fw<2||bw<2)continue;
    const arcA=prefix[j]-prefix[i],arcB=per-arcA;
    const w=Math.hypot(ring[i][0]-ring[j][0],ring[i][1]-ring[j][1]);
    if(w>4.5||w<0.01)continue;
    if(arcA<4.0||arcB<4.0)continue;
    pairs.push({i,j,w,arcA,arcB});
  }
  pairs.sort((a,b)=>a.w-b.w);
  for(const p of pairs){
    if(shown>=12)break;
    const inside=cutInside(ring,p.i,p.j);
    // 两部分面积
    const ringA=ring.slice(p.i,p.j+1);
    const ringB=ring.slice(p.j).concat(ring.slice(0,p.i+1));
    const aA=area(ringA),aB=area(ringB);
    const okW=p.w<=3.0, okA=aA>=20&&aB>=20;
    const verdict=(okW&&inside&&okA)?'>>> 应被切!' : `未切(宽${okW?'过':'超'}|线${inside?'内':'外'}|积${okA?'过':'小'})`;
    console.log(`  v${p.i}~v${p.j}: 宽=${p.w.toFixed(2)}m 弧=${p.arcA.toFixed(0)}/${p.arcB.toFixed(0)}m 切线在${inside?'内':'外'} 积=${aA.toFixed(0)}/${aB.toFixed(0)}m² ${verdict}`);
    shown++;
  }
  if(!pairs.length)console.log('  (无满足 弧长>=4m 且 宽<=4.5m 的顶点对)');
}
