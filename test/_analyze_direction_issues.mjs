// 综合分析: 方向检测问题 id 的证据对齐
// 对每个 id 输出: [检测系统] [用户标注] [初始轮廓长边方向直方图] [输出边方向直方图]
import { readFileSync } from 'fs';

// ---------- DBF ----------
function parseDbf(base) {
  const dbf = readFileSync(base + '.dbf');
  const n = dbf.readUInt16LE(4), headerLen = dbf.readUInt16LE(8), recLen = dbf.readUInt16LE(10);
  const fields = []; let p = 32;
  while (dbf[p] !== 0x0D) {
    fields.push({ name: dbf.toString('latin1', p, p + 11).replace(/\0+$/, '').trim(), len: dbf[p + 16] });
    p += 32;
  }
  const rows = [];
  for (let r = 0; r < n; r++) {
    const off = headerLen + r * recLen; const row = {}; let q = 0;
    for (const f of fields) { row[f.name] = dbf.toString('latin1', off + 1 + q, off + 1 + q + f.len).trim(); q += f.len; }
    rows.push(row);
  }
  return rows;
}

// ---------- SHP ----------
function parseShp(base) {
  const shp = readFileSync(base + '.shp');
  const geoms = []; let off = 100;
  while (off + 8 <= shp.length) {
    const contentLen = shp.readInt32BE(off + 4) * 2;
    const body = shp.subarray(off + 8, off + 8 + contentLen);
    if (body.length < 4) break;
    const type = body.readInt32LE(0);
    if (type === 5 || type === 15 || type === 25) { // Polygon[M,Z]
      let q = 4 + 32;
      const numParts = body.readInt32LE(q); q += 4;
      const numPoints = body.readInt32LE(q); q += 4;
      const parts = [];
      for (let i = 0; i < numParts; i++) { parts.push(body.readInt32LE(q)); q += 4; }
      parts.push(numPoints);
      const rings = [];
      for (let i = 0; i < numParts; i++) {
        const ring = [];
        for (let j = parts[i]; j < parts[i + 1]; j++) { ring.push([body.readDoubleLE(q), body.readDoubleLE(q + 8)]); q += 16; }
        rings.push(ring);
      }
      geoms.push({ rings });
    } else if (type === 3 || type === 13 || type === 23) { // Polyline[M,Z]
      let q = 4 + 32;
      const numParts = body.readInt32LE(q); q += 4;
      const numPoints = body.readInt32LE(q); q += 4;
      const parts = [];
      for (let i = 0; i < numParts; i++) { parts.push(body.readInt32LE(q)); q += 4; }
      parts.push(numPoints);
      const lines = [];
      for (let i = 0; i < numParts; i++) {
        const line = [];
        for (let j = parts[i]; j < parts[i + 1]; j++) { line.push([body.readDoubleLE(q), body.readDoubleLE(q + 8)]); q += 16; }
        lines.push(line);
      }
      geoms.push({ lines });
    } else {
      geoms.push(null);
    }
    off += 8 + contentLen;
  }
  return geoms;
}

const fold = deg => ((deg % 90) + 90) % 90;
const circDist = (a, b) => { let d = Math.abs(a - b) % 90; return d > 45 ? 90 - d : d; };

// 环的边方向直方图: 返回按长度排序的边列表 [{deg, len}] (折叠角)
function edgeHistogram(ring) {
  const edges = [];
  for (let i = 0; i < ring.length - 1; i++) {
    const [x1, y1] = ring[i], [x2, y2] = ring[i + 1];
    const len = Math.hypot(x2 - x1, y2 - y1);
    if (len < 0.5) continue;
    const deg = fold(Math.atan2(y2 - y1, x2 - x1) * 180 / Math.PI);
    edges.push({ deg, len });
  }
  edges.sort((a, b) => b.len - a.len);
  return edges;
}

// 把边按 8° 聚类成方向簇(简单贪心)
function clusterEdges(edges, tolDeg = 8) {
  const clusters = [];
  for (const e of edges) {
    let best = null;
    for (const c of clusters) {
      if (circDist(e.deg, c.mean) <= tolDeg) { best = c; break; }
    }
    if (best) {
      const w = best.totalLen + e.len;
      best.mean = (best.mean * best.totalLen + e.deg * e.len) / w;
      best.totalLen = w; best.count++;
    } else {
      clusters.push({ mean: e.deg, totalLen: e.len, count: 1 });
    }
  }
  clusters.sort((a, b) => b.totalLen - a.totalLen);
  return clusters;
}

const fmtC = cs => cs.map(c => `${c.mean.toFixed(1)}°(${c.totalLen.toFixed(0)}m/${c.count}条)`).join(' ');

// ---------- 1. 方向检测 + 用户标注 ----------
const dirRows = parseDbf('debug_direction_systems');
const dirGeoms = parseShp('debug_direction_systems');
const byId = new Map();       // id -> 程序检测系统列表
const userMarks = new Map();  // id -> [{rank, angle}]
dirRows.forEach((r, i) => {
  const id = r.src_fid, g = dirGeoms[i];
  const isUser = r.weight.includes('*');
  // 几何反算角度(所有线都算, 用户的只能从几何拿)
  let ang = null;
  const lines = g?.lines;
  if (lines && lines[0]?.length >= 2) {
    const a = lines[0][0], b = lines[0][lines[0].length - 1];
    ang = fold(Math.atan2(b[1] - a[1], b[0] - a[0]) * 180 / Math.PI);
  }
  if (isUser) {
    if (!userMarks.has(id)) userMarks.set(id, []);
    userMarks.get(String(id)).push({ rank: +r.rank, angle: ang });
  } else {
    if (!byId.has(id)) byId.set(id, []);
    byId.get(String(id)).push({ rank: +r.rank, angle: +r.angle_deg, weight: +r.weight, conf: +r.conf, chains: +r.chains, valid: +r.valid, reject: r.reject });
  }
});

// ---------- 2. 初始轮廓 / 输出 ----------
const initRows = parseDbf('initial_building_outline');
const initGeoms = parseShp('initial_building_outline');
const outRows = parseDbf('regularized_building_1');
const outGeoms = parseShp('regularized_building_1');
const initById = new Map(), outById = new Map();
initRows.forEach((r, i) => { if (!initById.has(r.id)) initById.set(r.id, initGeoms[i]); });
outRows.forEach((r, i) => { if (!outById.has(r.id)) outById.set(r.id, outGeoms[i]); });

// ---------- 3. 问题1: 拒绝案例 ----------
const rejected = [...byId.entries()].filter(([, v]) => +v[0].valid === 0).map(([id]) => id);
console.log('===== 问题1: 检测被拒绝的 id =====');
console.log(rejected.join(', '));

// ---------- 4. 逐 id 报告 ----------
function report(id, tag) {
  console.log(`\n===== id=${id} ${tag} =====`);
  const sys = byId.get(String(id)) ?? [];
  console.log('  检测系统: ' + (sys.length ? sys.map(s =>
    `rank${s.rank}=${s.angle.toFixed(1)}° (w=${s.weight.toFixed(0)} conf=${s.conf.toFixed(2)} chains=${s.chains}${s.reject ? ' REJECT:' + s.reject : ''})`).join(' | ') : '无记录'));
  const marks = userMarks.get(String(id));
  if (marks) console.log('  用户标注: ' + marks.map(m => `rank${m.rank}=${m.angle.toFixed(1)}°`).join(' | '));
  const initGeom = initById.get(String(id));
  if (initGeom?.rings) {
    const hist = edgeHistogram(initGeom.rings[0]);
    const cl = clusterEdges(hist);
    console.log(`  初始轮廓边簇(总长${hist.reduce((s, e) => s + e.len, 0).toFixed(0)}m): ` + fmtC(cl.slice(0, 5)));
  }
  const outGeom = outById.get(String(id));
  if (outGeom?.rings) {
    const hist = edgeHistogram(outGeom.rings[0]);
    const cl = clusterEdges(hist);
    console.log(`  输出边簇(总长${hist.reduce((s, e) => s + e.len, 0).toFixed(0)}m): ` + fmtC(cl.slice(0, 5)));
  }
}

for (const id of [43, 78, 85, 96, 124, 132, 135, 141, 188, 80, 152, 257]) report(id, '问题1: 检测被拒绝');
for (const id of [4, 13, 28]) report(id, '问题2: 多方向检测但输出缺第二方向');
for (const id of [140, 130, 142, 209, 234, 159, 103, 45, 48, 300, 274, 254]) report(id, '问题3: 单方向误判为多方向');
for (const id of [291, 120]) report(id, '问题4: 多方向漏检第二方向');
report(120, '问题5: 主方向不准');
