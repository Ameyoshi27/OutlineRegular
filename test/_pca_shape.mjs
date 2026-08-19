// 临时分析:初始轮廓整体主轴(PCA) + 加权方向拟合,判断真实建筑朝向
import { readFileSync } from 'fs';

function parseShp(path) {
  const buf = readFileSync(path);
  const rings = [];
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
      for (let i = 0; i < numParts; i++) {
        const ring = [];
        for (let j = parts[i]; j < parts[i + 1]; j++) { ring.push([body.readDoubleLE(p), body.readDoubleLE(p + 8)]); p += 16; }
        if (ring.length > 1) {
          const a = ring[0], b = ring[ring.length - 1];
          if (Math.abs(a[0]-b[0]) < 1e-9 && Math.abs(a[1]-b[1]) < 1e-9) ring.pop();
        }
        if (ring.length >= 3) rings.push(ring);
      }
    }
    off += 8 + contentLen;
  }
  return rings;
}

for (const path of process.argv.slice(2)) {
  for (const ring of parseShp(path)) {
    // PCA 主轴
    let mx = 0, my = 0;
    for (const [x, y] of ring) { mx += x; my += y; }
    mx /= ring.length; my /= ring.length;
    let sxx = 0, syy = 0, sxy = 0;
    for (const [x, y] of ring) { const dx = x - mx, dy = y - my; sxx += dx * dx; syy += dy * dy; sxy += dx * dy; }
    const theta = 0.5 * Math.atan2(2 * sxy, sxx - syy) * 180 / Math.PI;
    const eigRatio = Math.sqrt(Math.max(0, (sxx + syy - Math.hypot(sxx - syy, 2 * sxy)) / 2)) /
      Math.sqrt((sxx + syy + Math.hypot(sxx - syy, 2 * sxy)) / 2);
    let t = ((theta % 90) + 90) % 90;
    // 楼梯消除后的"整体边方向": 长度加权,先把每条边按角度模90聚合,再用高斯平滑看峰
    const bins = new Float64Array(90);
    for (let i = 0; i < ring.length; i++) {
      const [x1, y1] = ring[i], [x2, y2] = ring[(i + 1) % ring.length];
      const len = Math.hypot(x2 - x1, y2 - y1);
      if (len < 0.01) continue;
      let ang = ((Math.atan2(y2 - y1, x2 - x1) * 180 / Math.PI) % 90 + 90) % 90;
      bins[Math.min(89, Math.floor(ang))] += len;
    }
    const smooth = new Float64Array(90);
    for (let i = 0; i < 90; i++) for (let d = -3; d <= 3; d++) smooth[i] += bins[((i + d) % 90 + 90) % 90];
    let best = 0, bi = 0;
    smooth.forEach((v, i) => { if (v > best) { best = v; bi = i; } });
    // 累计转角(整体形状是否平滑旋转)
    console.log(`${path}: verts=${ring.length} PCA主轴=${t.toFixed(1)}°(mod90) 长轴/短轴=${(1 / Math.max(eigRatio, 1e-9)).toFixed(2)} 边直方图峰=${bi}°`);
  }
}
