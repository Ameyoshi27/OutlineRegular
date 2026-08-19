// 临时脚本：解析 shapefile（无依赖），输出多边形几何为可读文本
import { readFileSync } from 'fs';

const SHP_TYPE = { 0: 'Null', 1: 'Point', 3: 'Polyline', 5: 'Polygon', 8: 'MultiPoint' };

function parse(path) {
  const buf = readFileSync(path);
  const records = [];
  let off = 100; // skip header
  while (off + 8 <= buf.length) {
    const num = buf.readInt32BE(off);
    const contentLen = buf.readInt32BE(off + 4) * 2; // words -> bytes
    const body = buf.subarray(off + 8, off + 8 + contentLen);
    if (body.length < 4) break;
    const type = body.readInt32LE(0);
    if (contentLen <= 0 || off + 8 + contentLen > buf.length) {
      console.log(`  [warn] record#${num} bad contentLen=${contentLen} at off=${off}, fileLen=${buf.length}; stop`);
      break;
    }
    if (type === 5 || type === 15 || type === 25) { // (Multi)Polygon[M,Z]
      // ESRI 白皮书：Box, NumParts, NumPoints, Parts, Points(XY), 然后才是 Z/M 数组
      let p = 4;
      const box = [];
      for (let i = 0; i < 4; i++) { box.push(body.readDoubleLE(p)); p += 8; }
      const numParts = body.readInt32LE(p); p += 4;
      const numPoints = body.readInt32LE(p); p += 4;
      const parts = [];
      for (let i = 0; i < numParts; i++) { parts.push(body.readInt32LE(p)); p += 4; }
      parts.push(numPoints);
      const rings = [];
      for (let i = 0; i < numParts; i++) {
        const ring = [];
        for (let j = parts[i]; j < parts[i + 1]; j++) {
          const x = body.readDoubleLE(p); p += 8;
          const y = body.readDoubleLE(p); p += 8;
          ring.push([x, y]);
        }
        rings.push(ring);
      }
      records.push({ num, box, rings });
    } else {
      records.push({ num, type: SHP_TYPE[type] || type, skipped: true });
    }
    off += 8 + contentLen;
  }
  return records;
}

const [,, ...paths] = process.argv;
for (const path of paths) {
  console.log('=====', path);
  const recs = parse(path);
  console.log(`features: ${recs.length}`);
  recs.forEach((r, i) => {
    if (r.skipped) { console.log(`[${i}] type=${r.type} (skipped)`); return; }
    console.log(`[${i}] parts=${r.rings.length} box=[${r.box.map(v => v.toFixed(1)).join(', ')}]`);
    r.rings.forEach((ring, k) => {
      // 输出每条边的长度与方向角，便于看主方向结构
      let edges = [];
      for (let j = 0; j + 1 < ring.length; j++) {
        const [x1, y1] = ring[j], [x2, y2] = ring[j + 1];
        const dx = x2 - x1, dy = y2 - y1;
        const len = Math.hypot(dx, dy);
        let ang = Math.atan2(dy, dx) * 180 / Math.PI;
        if (ang < 0) ang += 180; // 无向边
        edges.push({ len, ang });
      }
      edges = edges.filter(e => e.len > 0.05);
      const total = edges.reduce((s, e) => s + e.len, 0);
      console.log(`  ring${k}: ${ring.length} pts, ${edges.length} edges (>${0.05}m), perim=${total.toFixed(1)}m`);
      // 按方向聚类统计（10°桶）
      const buckets = new Map();
      for (const e of edges) {
        const b = Math.round(e.ang / 10) * 10 % 180;
        const key = `${b}°`;
        buckets.set(key, (buckets.get(key) || 0) + e.len);
      }
      const sorted = [...buckets.entries()].sort((a, b) => b[1] - a[1]);
      console.log(`  方向直方图(按边长): ` + sorted.slice(0, 6).map(([k, v]) => `${k}:${v.toFixed(1)}m(${(100 * v / total).toFixed(0)}%)`).join('  '));
      // 打印所有 >1m 的边
      edges.filter(e => e.len > 1.0).forEach(e => console.log(`    edge len=${e.len.toFixed(2)}m ang=${e.ang.toFixed(1)}°`));
    });
  });
}
