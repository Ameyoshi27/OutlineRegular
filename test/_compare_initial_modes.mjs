// Compare two initial-outline Shapefiles in their common spatial extent.
// Dependency-free: node test/_compare_initial_modes.mjs old.shp new.shp
import { readFileSync } from 'fs';
import { basename, dirname, join } from 'path';

function parseShp(path) {
  const buffer = readFileSync(path);
  const records = [];
  let offset = 100;
  while (offset + 8 <= buffer.length) {
    const bytes = buffer.readInt32BE(offset + 4) * 2;
    if (bytes <= 0 || offset + 8 + bytes > buffer.length) break;
    const body = buffer.subarray(offset + 8, offset + 8 + bytes);
    const type = body.readInt32LE(0);
    let rings = [];
    if (type === 5 || type === 15 || type === 25) {
      let p = 36;
      const partCount = body.readInt32LE(p); p += 4;
      const pointCount = body.readInt32LE(p); p += 4;
      const starts = [];
      for (let i = 0; i < partCount; ++i) { starts.push(body.readInt32LE(p)); p += 4; }
      starts.push(pointCount);
      for (let part = 0; part < partCount; ++part) {
        const ring = [];
        for (let i = starts[part]; i < starts[part + 1]; ++i) {
          ring.push([body.readDoubleLE(p), body.readDoubleLE(p + 8)]);
          p += 16;
        }
        if (ring.length > 1 && Math.hypot(
          ring[0][0] - ring.at(-1)[0], ring[0][1] - ring.at(-1)[1]) < 1e-8) ring.pop();
        if (ring.length >= 3) rings.push(ring);
      }
    }
    records.push({ rings });
    offset += 8 + bytes;
  }
  return records;
}

function parseDbf(path) {
  const buffer = readFileSync(path);
  const count = buffer.readUInt32LE(4);
  const headerSize = buffer.readUInt16LE(8);
  const recordSize = buffer.readUInt16LE(10);
  const fields = [];
  for (let p = 32; p + 32 <= headerSize && buffer[p] !== 0x0d; p += 32) {
    fields.push({
      name: buffer.toString('ascii', p, p + 11).replace(/\0.*$/, '').trim(),
      type: String.fromCharCode(buffer[p + 11]),
      length: buffer[p + 16]
    });
  }
  const records = [];
  for (let i = 0; i < count; ++i) {
    const base = headerSize + i * recordSize;
    if (base + recordSize > buffer.length) break;
    const row = { _deleted: buffer[base] === 0x2a };
    let p = base + 1;
    for (const field of fields) {
      const raw = buffer.toString('latin1', p, p + field.length).trim();
      row[field.name] = field.type === 'N' || field.type === 'F'
        ? (raw ? Number(raw) : null) : raw;
      p += field.length;
    }
    records.push(row);
  }
  return records;
}

function load(path) {
  const shapes = parseShp(path);
  const dbf = parseDbf(join(dirname(path), `${basename(path, '.shp')}.dbf`));
  const features = [];
  for (let i = 0; i < Math.min(shapes.length, dbf.length); ++i) {
    if (dbf[i]._deleted || !shapes[i].rings.length) continue;
    const ring = shapes[i].rings[0];
    features.push({ record: i, ring, attrs: dbf[i], box: bbox(ring), area: area(ring) });
  }
  return features;
}

function area(ring) {
  return Math.abs(ring.reduce((sum, p, i) => {
    const q = ring[(i + 1) % ring.length];
    return sum + p[0] * q[1] - q[0] * p[1];
  }, 0) / 2);
}

function perimeter(ring) {
  return ring.reduce((sum, p, i) => {
    const q = ring[(i + 1) % ring.length];
    return sum + Math.hypot(q[0] - p[0], q[1] - p[1]);
  }, 0);
}

function bbox(ring) {
  const box = { minX: Infinity, minY: Infinity, maxX: -Infinity, maxY: -Infinity };
  for (const [x, y] of ring) {
    box.minX = Math.min(box.minX, x); box.minY = Math.min(box.minY, y);
    box.maxX = Math.max(box.maxX, x); box.maxY = Math.max(box.maxY, y);
  }
  return box;
}

const intersects = (a, b, pad = 0) =>
  a.maxX + pad >= b.minX && b.maxX + pad >= a.minX &&
  a.maxY + pad >= b.minY && b.maxY + pad >= a.minY;

function pointInRing([x, y], ring) {
  let inside = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const [xi, yi] = ring[i], [xj, yj] = ring[j];
    if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

function approximateIntersection(a, b, step = 0.5) {
  const box = {
    minX: Math.max(a.box.minX, b.box.minX), minY: Math.max(a.box.minY, b.box.minY),
    maxX: Math.min(a.box.maxX, b.box.maxX), maxY: Math.min(a.box.maxY, b.box.maxY)
  };
  if (box.minX >= box.maxX || box.minY >= box.maxY) return 0;
  let count = 0;
  for (let y = box.minY + step / 2; y < box.maxY; y += step) {
    for (let x = box.minX + step / 2; x < box.maxX; x += step) {
      if (pointInRing([x, y], a.ring) && pointInRing([x, y], b.ring)) ++count;
    }
  }
  return count * step * step;
}

function segmentProjectionOverlap(a0, a1, b0, b1, tolerance = 0.65) {
  const adx = a1[0] - a0[0], ady = a1[1] - a0[1];
  const bdx = b1[0] - b0[0], bdy = b1[1] - b0[1];
  const al = Math.hypot(adx, ady), bl = Math.hypot(bdx, bdy);
  if (al < 1e-8 || bl < 1e-8) return 0;
  const ux = adx / al, uy = ady / al;
  if (Math.abs((adx * bdx + ady * bdy) / (al * bl)) < Math.cos(8 * Math.PI / 180)) return 0;
  const lineDistance = p => Math.abs((p[0] - a0[0]) * uy - (p[1] - a0[1]) * ux);
  if (Math.min(lineDistance(b0), lineDistance(b1)) > tolerance) return 0;
  const t0 = (b0[0] - a0[0]) * ux + (b0[1] - a0[1]) * uy;
  const t1 = (b1[0] - a0[0]) * ux + (b1[1] - a0[1]) * uy;
  return Math.max(0, Math.min(al, Math.max(t0, t1)) - Math.max(0, Math.min(t0, t1)));
}

function robustSharedLength(a, b) {
  let length = 0;
  for (let i = 0; i < a.ring.length; ++i) {
    const a0 = a.ring[i], a1 = a.ring[(i + 1) % a.ring.length];
    for (let j = 0; j < b.ring.length; ++j) {
      length += segmentProjectionOverlap(a0, a1, b.ring[j], b.ring[(j + 1) % b.ring.length]);
    }
  }
  return length;
}

function quantile(values, q) {
  if (!values.length) return NaN;
  const sorted = [...values].sort((a, b) => a - b);
  return sorted[Math.min(sorted.length - 1, Math.floor(q * sorted.length))];
}

const [,, oldPath, currentPath] = process.argv;
if (!oldPath || !currentPath) throw new Error('usage: node _compare_initial_modes.mjs old.shp current.shp');
const oldAll = load(oldPath);
const currentAll = load(currentPath);
const extent = currentAll.reduce((box, f) => ({
  minX: Math.min(box.minX, f.box.minX), minY: Math.min(box.minY, f.box.minY),
  maxX: Math.max(box.maxX, f.box.maxX), maxY: Math.max(box.maxY, f.box.maxY)
}), { minX: Infinity, minY: Infinity, maxX: -Infinity, maxY: -Infinity });
const isNearBlackParent = feature => {
  const value = Number(feature.attrs.parent);
  if (!Number.isFinite(value)) return false;
  return ((value >> 16) & 255) <= 2 && ((value >> 8) & 255) <= 2 && (value & 255) <= 2;
};
const nearBlack = currentAll.filter(isNearBlackParent);
const current = currentAll.filter(feature => !isNearBlackParent(feature));
const old = oldAll.filter(f => {
  const x = (f.box.minX + f.box.maxX) / 2, y = (f.box.minY + f.box.maxY) / 2;
  return x >= extent.minX && x <= extent.maxX && y >= extent.minY && y <= extent.maxY;
});

console.log(`old_all=${oldAll.length} old_common=${old.length} current_raw=${currentAll.length} current_without_near_black=${current.length}`);
console.log(`current_near_black_parent_features=${nearBlack.length} area_sum=${nearBlack.reduce((sum, f) => sum + f.area, 0).toFixed(1)}`);
console.log(`common_extent=${extent.minX.toFixed(2)},${extent.minY.toFixed(2)} .. ${extent.maxX.toFixed(2)},${extent.maxY.toFixed(2)}`);

const pairs = [];
for (let i = 0; i < old.length; ++i) for (let j = 0; j < current.length; ++j) {
  if (!intersects(old[i].box, current[j].box)) continue;
  const intersection = approximateIntersection(old[i], current[j]);
  if (intersection < 1) continue;
  const iou = intersection / Math.max(1e-9, old[i].area + current[j].area - intersection);
  const smallerCoverage = intersection / Math.max(1e-9, Math.min(old[i].area, current[j].area));
  pairs.push({ i, j, intersection, iou, smallerCoverage });
}

const greedy = [...pairs].sort((a, b) => b.iou - a.iou);
const usedOld = new Set(), usedCurrent = new Set(), matches = [];
for (const pair of greedy) {
  if (usedOld.has(pair.i) || usedCurrent.has(pair.j) || pair.iou < 0.05) continue;
  usedOld.add(pair.i); usedCurrent.add(pair.j); matches.push(pair);
}
console.log(`greedy_matches=${matches.length} iou_median=${quantile(matches.map(p => p.iou), 0.5).toFixed(3)} iou_p25=${quantile(matches.map(p => p.iou), 0.25).toFixed(3)}`);
console.log(`matches_iou_ge_0.8=${matches.filter(p => p.iou >= 0.8).length}`);

const oldLinks = Array.from({ length: old.length }, () => []);
const currentLinks = Array.from({ length: current.length }, () => []);
for (const pair of pairs) if (pair.smallerCoverage >= 0.35) {
  oldLinks[pair.i].push(pair);
  currentLinks[pair.j].push(pair);
}
console.log(`old_to_multiple_current=${oldLinks.filter(v => v.length >= 2).length}`);
console.log(`current_to_multiple_old=${currentLinks.filter(v => v.length >= 2).length}`);
console.log(`old_without_35pct_overlap=${oldLinks.filter(v => v.length === 0).length}`);
console.log(`current_without_35pct_overlap=${currentLinks.filter(v => v.length === 0).length}`);

let adjacent = 0, sameParent = 0, longSeam = 0, shapeGate = 0, ambiguous = 0;
for (let i = 0; i < current.length; ++i) for (let j = i + 1; j < current.length; ++j) {
  // Match the production merge prefilter: envelopes must touch/overlap before
  // the 0.65 m robust shared-boundary tolerance is evaluated.
  if (!intersects(current[i].box, current[j].box)) continue;
  const shared = robustSharedLength(current[i], current[j]);
  if (shared < 0.6) continue;
  ++adjacent;
  const parentA = Number(current[i].attrs.parent ?? 0);
  const parentB = Number(current[j].attrs.parent ?? 0);
  const parentMatch = parentA > 0 && parentA === parentB;
  if (parentMatch) ++sameParent;
  if (shared >= 3) ++longSeam;
  const sharedRatio = shared / Math.max(1e-9, Math.min(perimeter(current[i].ring), perimeter(current[j].ring)));
  const passesShapeEntry = shared >= 2 || sharedRatio >= 0.015;
  if (passesShapeEntry) ++shapeGate;
  if (!parentMatch && shared >= 3 && Math.min(current[i].area, current[j].area) <= 1500) ++ambiguous;
}
console.log(`current_adjacent_ge_0.6m=${adjacent}`);
console.log(`current_same_parent_adjacent=${sameParent}`);
console.log(`current_shared_ge_3m=${longSeam}`);
console.log(`current_shape_entry_gate=${shapeGate}`);
console.log(`current_osgb_ambiguous_pairs=${ambiguous}`);
