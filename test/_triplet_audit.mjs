// Compare mask outline, algorithm output and a manually corrected reference.
// Dependency-free so it can be rerun with the project's existing Node setup.
import { readFileSync, existsSync } from 'fs';
import { dirname, basename, join } from 'path';

function parseShp(path) {
  const buf = readFileSync(path);
  const records = [];
  let off = 100;
  while (off + 8 <= buf.length) {
    const bytes = buf.readInt32BE(off + 4) * 2;
    if (bytes <= 0 || off + 8 + bytes > buf.length) break;
    const body = buf.subarray(off + 8, off + 8 + bytes);
    const type = body.readInt32LE(0);
    let rings = null;
    if (type === 5 || type === 15 || type === 25) {
      let p = 36;
      const partCount = body.readInt32LE(p); p += 4;
      const pointCount = body.readInt32LE(p); p += 4;
      const starts = [];
      for (let i = 0; i < partCount; ++i) { starts.push(body.readInt32LE(p)); p += 4; }
      starts.push(pointCount);
      rings = [];
      for (let part = 0; part < partCount; ++part) {
        const ring = [];
        for (let i = starts[part]; i < starts[part + 1]; ++i) {
          ring.push([body.readDoubleLE(p), body.readDoubleLE(p + 8)]);
          p += 16;
        }
        if (ring.length > 1) {
          const a = ring[0], b = ring[ring.length - 1];
          if (Math.hypot(a[0] - b[0], a[1] - b[1]) < 1e-8) ring.pop();
        }
        if (ring.length >= 3) rings.push(ring);
      }
    }
    records.push({ rings });
    off += 8 + bytes;
  }
  return records;
}

function parseDbf(path) {
  const buf = readFileSync(path);
  const count = buf.readUInt32LE(4);
  const headerSize = buf.readUInt16LE(8);
  const recordSize = buf.readUInt16LE(10);
  const fields = [];
  let p = 32;
  while (p + 32 <= headerSize && buf[p] !== 0x0d) {
    fields.push({
      name: buf.toString('ascii', p, p + 11).replace(/\0.*$/, '').trim(),
      type: String.fromCharCode(buf[p + 11]),
      length: buf[p + 16],
      decimals: buf[p + 17]
    });
    p += 32;
  }
  const records = [];
  for (let i = 0; i < count; ++i) {
    const base = headerSize + i * recordSize;
    if (base + recordSize > buf.length) break;
    const row = { _deleted: buf[base] === 0x2a };
    let q = base + 1;
    for (const field of fields) {
      const raw = buf.toString('latin1', q, q + field.length).trim();
      row[field.name] = field.type === 'N' || field.type === 'F'
        ? (raw === '' ? null : Number(raw))
        : raw;
      q += field.length;
    }
    records.push(row);
  }
  return { fields, records };
}

function loadDataset(shpPath, label) {
  const dbfPath = join(dirname(shpPath), basename(shpPath, '.shp') + '.dbf');
  const shapes = parseShp(shpPath);
  const dbf = parseDbf(dbfPath);
  const features = [];
  for (let i = 0; i < Math.min(shapes.length, dbf.records.length); ++i) {
    if (!shapes[i].rings?.length || dbf.records[i]._deleted) continue;
    const attrs = dbf.records[i];
    const idField = ['id', 'ID', 'Id', 'mask', 'MASK', 'source_fid', 'SOURCE_FID']
      .find(name => attrs[name] !== undefined && attrs[name] !== null && attrs[name] !== '');
    features.push({
      record: i,
      id: idField ? String(attrs[idField]) : String(i),
      idField: idField ?? 'record',
      attrs,
      rings: shapes[i].rings,
      ring: shapes[i].rings[0]
    });
  }
  console.log(`${label}: features=${features.length}, fields=${dbf.fields.map(f => f.name).join(',')}, id_field=${features[0]?.idField ?? 'none'}`);
  return features;
}

const signedArea = ring => ring.reduce((sum, p, i) => {
  const q = ring[(i + 1) % ring.length];
  return sum + p[0] * q[1] - q[0] * p[1];
}, 0) / 2;
const area = ring => Math.abs(signedArea(ring));
const perimeter = ring => ring.reduce((sum, p, i) => {
  const q = ring[(i + 1) % ring.length];
  return sum + Math.hypot(q[0] - p[0], q[1] - p[1]);
}, 0);

function bbox(ring) {
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  for (const [x, y] of ring) {
    minX = Math.min(minX, x); minY = Math.min(minY, y);
    maxX = Math.max(maxX, x); maxY = Math.max(maxY, y);
  }
  return { minX, minY, maxX, maxY };
}

function pointInRing([x, y], ring) {
  let inside = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const [xi, yi] = ring[i], [xj, yj] = ring[j];
    if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

function pointSegmentDistance(p, a, b) {
  const dx = b[0] - a[0], dy = b[1] - a[1];
  const d2 = dx * dx + dy * dy;
  const t = d2 > 1e-12
    ? Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / d2))
    : 0;
  return Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy);
}

function distanceToRing(p, ring) {
  let best = Infinity;
  for (let i = 0; i < ring.length; ++i) {
    best = Math.min(best, pointSegmentDistance(p, ring[i], ring[(i + 1) % ring.length]));
  }
  return best;
}

function boundarySamples(ring, spacing = 0.5) {
  const out = [];
  for (let i = 0; i < ring.length; ++i) {
    const a = ring[i], b = ring[(i + 1) % ring.length];
    const len = Math.hypot(b[0] - a[0], b[1] - a[1]);
    const count = Math.max(1, Math.ceil(len / spacing));
    for (let k = 0; k < count; ++k) {
      const t = k / count;
      out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
    }
  }
  return out;
}

function directedBoundaryStats(from, to) {
  const distances = boundarySamples(from).map(p => distanceToRing(p, to)).sort((a, b) => a - b);
  const mean = distances.reduce((sum, value) => sum + value, 0) / Math.max(1, distances.length);
  const at = q => distances[Math.min(distances.length - 1, Math.floor(q * distances.length))] ?? 0;
  return { mean, q90: at(0.9), q95: at(0.95), max: distances.at(-1) ?? 0 };
}

function approximateIou(a, b, step = 0.30) {
  const ba = bbox(a), bb = bbox(b);
  const minX = Math.min(ba.minX, bb.minX), minY = Math.min(ba.minY, bb.minY);
  const maxX = Math.max(ba.maxX, bb.maxX), maxY = Math.max(ba.maxY, bb.maxY);
  let intersection = 0, union = 0;
  for (let y = minY + step / 2; y < maxY; y += step) {
    for (let x = minX + step / 2; x < maxX; x += step) {
      const inA = pointInRing([x, y], a), inB = pointInRing([x, y], b);
      if (inA || inB) ++union;
      if (inA && inB) ++intersection;
    }
  }
  return union ? intersection / union : 0;
}

function lineAngle90(a, b) {
  const degrees = Math.atan2(b[1] - a[1], b[0] - a[0]) * 180 / Math.PI;
  return ((degrees % 90) + 90) % 90;
}

function circularDifference90(a, b) {
  const d = Math.abs(a - b);
  return Math.min(d, 90 - d);
}

function directionPeaks(ring) {
  const bins = new Array(90).fill(0);
  let total = 0;
  for (let i = 0; i < ring.length; ++i) {
    const a = ring[i], b = ring[(i + 1) % ring.length];
    const len = Math.hypot(b[0] - a[0], b[1] - a[1]);
    if (len < 1.0) continue;
    bins[Math.floor(lineAngle90(a, b))] += len;
    total += len;
  }
  const smooth = bins.map((_, center) => {
    let value = 0;
    for (let d = -4; d <= 4; ++d) value += bins[(center + d + 90) % 90];
    return value;
  });
  const peaks = [];
  for (let rank = 0; rank < 3; ++rank) {
    let best = -1, bestValue = 0;
    for (let i = 0; i < 90; ++i) {
      if (peaks.some(p => circularDifference90(i, p.angle) < 12)) continue;
      if (smooth[i] > bestValue) { best = i; bestValue = smooth[i]; }
    }
    if (best < 0 || bestValue < 0.08 * total) break;
    peaks.push({ angle: best, ratio: total ? bestValue / total : 0 });
  }
  return peaks;
}

function vertexStats(ring) {
  let short = 0, nearlyCollinear = 0, sharp = 0, reflex = 0;
  const orientation = Math.sign(signedArea(ring)) || 1;
  for (let i = 0; i < ring.length; ++i) {
    const a = ring[(i + ring.length - 1) % ring.length];
    const b = ring[i], c = ring[(i + 1) % ring.length];
    const l1 = Math.hypot(b[0] - a[0], b[1] - a[1]);
    const l2 = Math.hypot(c[0] - b[0], c[1] - b[1]);
    if (l2 < 2.0) ++short;
    const ux = (a[0] - b[0]) / Math.max(l1, 1e-12);
    const uy = (a[1] - b[1]) / Math.max(l1, 1e-12);
    const vx = (c[0] - b[0]) / Math.max(l2, 1e-12);
    const vy = (c[1] - b[1]) / Math.max(l2, 1e-12);
    const angle = Math.acos(Math.max(-1, Math.min(1, ux * vx + uy * vy))) * 180 / Math.PI;
    if (angle > 165) ++nearlyCollinear;
    if (angle < 35) ++sharp;
    const cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
    if (cross * orientation < 0) ++reflex;
  }
  return { short, nearlyCollinear, sharp, reflex };
}

function compare(from, target) {
  const a = area(from), b = area(target);
  const forward = directedBoundaryStats(from, target);
  const reverse = directedBoundaryStats(target, from);
  return {
    iou: approximateIou(from, target),
    areaRatio: a / Math.max(b, 1e-9),
    mean: Math.max(forward.mean, reverse.mean),
    q90: Math.max(forward.q90, reverse.q90),
    hausdorff: Math.max(forward.max, reverse.max)
  };
}

function segmentsProperlyIntersect(a, b, c, d) {
  const cross = (p, q, r) => (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
  const x1 = cross(a, b, c), x2 = cross(a, b, d), x3 = cross(c, d, a), x4 = cross(c, d, b);
  return x1 * x2 < -1e-10 && x3 * x4 < -1e-10;
}

function overlapAreaApprox(a, b, step = 0.25) {
  const ba = bbox(a), bb = bbox(b);
  const minX = Math.max(ba.minX, bb.minX), minY = Math.max(ba.minY, bb.minY);
  const maxX = Math.min(ba.maxX, bb.maxX), maxY = Math.min(ba.maxY, bb.maxY);
  if (minX >= maxX || minY >= maxY) return 0;
  let count = 0;
  for (let y = minY + step / 2; y < maxY; y += step) {
    for (let x = minX + step / 2; x < maxX; x += step) {
      if (pointInRing([x, y], a) && pointInRing([x, y], b)) ++count;
    }
  }
  return count * step * step;
}

function overlaps(features) {
  const result = [];
  for (let i = 0; i < features.length; ++i) {
    const a = features[i].ring, ba = bbox(a);
    for (let j = i + 1; j < features.length; ++j) {
      const b = features[j].ring, bb = bbox(b);
      if (ba.maxX <= bb.minX || bb.maxX <= ba.minX || ba.maxY <= bb.minY || bb.maxY <= ba.minY) continue;
      let crosses = false;
      for (let x = 0; x < a.length && !crosses; ++x) {
        for (let y = 0; y < b.length; ++y) {
          if (segmentsProperlyIntersect(a[x], a[(x + 1) % a.length], b[y], b[(y + 1) % b.length])) {
            crosses = true; break;
          }
        }
      }
      const contained = pointInRing(a[0], b) || pointInRing(b[0], a);
      if (crosses || contained) {
        const overlap = overlapAreaApprox(a, b);
        if (overlap > 0.05) result.push({ a: features[i].id, b: features[j].id, overlap });
      }
    }
  }
  return result;
}

function extractLogBlocks(logPath, ids) {
  if (!logPath || !existsSync(logPath)) return new Map();
  const lines = readFileSync(logPath, 'utf8').split(/\r?\n/);
  const starts = [];
  for (let i = 0; i < lines.length; ++i) {
    const m = lines[i].match(/\[SupportEvidence\] source_fid=(-?\d+)/);
    if (m) starts.push({ index: i, id: m[1] });
  }
  const wanted = new Set(ids);
  const blocks = new Map();
  for (let i = 0; i < starts.length; ++i) {
    if (!wanted.has(starts[i].id)) continue;
    const end = i + 1 < starts.length ? starts[i + 1].index : lines.length;
    const selected = lines.slice(starts[i].index, end).filter(line =>
      /\[SupportEvidence\]|\[SupportFilter\]|\[RegTime\] ===|\[DirectionDecision\]|\[SingleFirst\]|\[BuildingMode\]|\[HypothesisRepair\]|\[MultiChain\]|\[RegularizationGuard\]|\[Fallback\]/.test(line));
    blocks.set(starts[i].id, selected);
  }
  return blocks;
}

const initialPath = process.argv[2];
const algorithmPath = process.argv[3];
const currentPath = process.argv[4];
const logPath = process.argv[5];
if (!initialPath || !algorithmPath || !currentPath) {
  console.error('usage: node test/_triplet_audit.mjs initial.shp algorithm.shp current.shp [console_output.txt]');
  process.exit(2);
}

const initial = loadDataset(initialPath, 'initial');
const algorithm = loadDataset(algorithmPath, 'algorithm');
const current = loadDataset(currentPath, 'current');
const runDir = logPath ? dirname(logPath) : '';
const hypothesisPath = runDir ? join(runDir, 'debug_best_hypothesis.shp') : '';
const fullInitialPath = runDir ? join(runDir, 'initial_building_outline.shp') : '';
const hypothesis = hypothesisPath && existsSync(hypothesisPath)
  ? loadDataset(hypothesisPath, 'hypothesis') : [];
const fullInitial = fullInitialPath && existsSync(fullInitialPath)
  ? loadDataset(fullInitialPath, 'full_initial') : [];
const groupById = values => {
  const groups = new Map();
  for (const value of values) {
    if (!groups.has(value.id)) groups.set(value.id, []);
    groups.get(value.id).push(value);
  }
  return groups;
};
const chooseForReference = (groups, id, reference) => {
  const candidates = groups.get(id) ?? [];
  if (candidates.length <= 1 || !reference) return candidates[0];
  let best = candidates[0], bestIou = -1;
  for (const candidate of candidates) {
    const value = approximateIou(candidate.ring, reference.ring, 0.5);
    if (value > bestIou) { best = candidate; bestIou = value; }
  }
  return best;
};
const initGroups = groupById(initial), curGroups = groupById(current);
const hypGroups = groupById(hypothesis), algGroups = groupById(algorithm);
const fullInitGroups = groupById(fullInitial);
const initById = new Map([...initGroups].map(([id, values]) => [id, values[0]]));
const curById = new Map([...curGroups].map(([id, values]) => [id, values[0]]));
// The algorithm input may be the full run while initial/current are a selected
// case set. Restrict the report to that case set.
const ids = [...new Set([...initById.keys(), ...curById.keys()])]
  .sort((a, b) => Number(a) - Number(b));
const hypById = new Map(ids.map(id => [id, chooseForReference(hypGroups, id, curById.get(id))]).filter(([, value]) => value));
const algById = new Map(ids.map(id => [id, chooseForReference(algGroups, id, curById.get(id))]).filter(([, value]) => value));
const fullInitById = new Map(ids.map(id => [id, chooseForReference(fullInitGroups, id, initById.get(id))]).filter(([, value]) => value));

const duplicateIds = groups => [...groups].filter(([, values]) => values.length > 1)
  .map(([id, values]) => `${id}x${values.length}`);
console.log(`duplicate ids: initial=${duplicateIds(initGroups).join(',') || 'none'}; algorithm=${duplicateIds(algGroups).join(',') || 'none'}; current=${duplicateIds(curGroups).join(',') || 'none'}; hypothesis=${duplicateIds(hypGroups).join(',') || 'none'}`);

console.log('\nPer-feature comparison (I=initial, H=best hypothesis, A=algorithm, C=corrected):');
console.log('id | verts I/H/A/C | area I/H/A/C | IoU I/C | H/C | A/C | area A/C | bQ90 A/C | Haus A/C | dirs A -> C | reflex A/C');
for (const id of ids) {
  const i = initById.get(id), h = hypById.get(id), a = algById.get(id), c = curById.get(id);
  if (!i || !c) {
    console.log(`${id} | missing initial/current`);
    continue;
  }
  const ic = compare(i.ring, c.ring);
  const hc = h ? compare(h.ring, c.ring) : null;
  if (!a) {
    console.log(`${id} | ${i.ring.length}/${h?.ring.length ?? 'MISS'}/MISS/${c.ring.length} | ${area(i.ring).toFixed(0)}/${h ? area(h.ring).toFixed(0) : 'MISS'}/MISS/${area(c.ring).toFixed(0)} | ${ic.iou.toFixed(3)} | ${hc?.iou.toFixed(3) ?? 'MISS'} | MISS | MISS | MISS | MISS | MISS | MISS`);
    continue;
  }
  const ac = compare(a.ring, c.ring);
  const dirsA = directionPeaks(a.ring).map(p => `${p.angle}(${Math.round(100 * p.ratio)}%)`).join(',');
  const dirsC = directionPeaks(c.ring).map(p => `${p.angle}(${Math.round(100 * p.ratio)}%)`).join(',');
  const va = vertexStats(a.ring), vc = vertexStats(c.ring);
  console.log(`${id} | ${i.ring.length}/${h?.ring.length ?? 'MISS'}/${a.ring.length}/${c.ring.length} | ${area(i.ring).toFixed(0)}/${h ? area(h.ring).toFixed(0) : 'MISS'}/${area(a.ring).toFixed(0)}/${area(c.ring).toFixed(0)} | ${ic.iou.toFixed(3)} | ${hc?.iou.toFixed(3) ?? 'MISS'} | ${ac.iou.toFixed(3)} | ${ac.areaRatio.toFixed(3)} | ${ac.q90.toFixed(2)} | ${ac.hausdorff.toFixed(2)} | ${dirsA} -> ${dirsC} | ${va.reflex}/${vc.reflex}`);
}

const matched = ids.filter(id => algById.has(id) && curById.has(id));
const metrics = matched.map(id => compare(algById.get(id).ring, curById.get(id).ring));
const sorted = (values, selector) => values.map(selector).sort((a, b) => a - b);
const median = values => values[Math.floor(values.length / 2)] ?? 0;
console.log('\nSummary:');
console.log(`matched=${matched.length}, missing_algorithm=${ids.filter(id => initById.has(id) && curById.has(id) && !algById.has(id)).join(',') || 'none'}`);
console.log(`IoU median=${median(sorted(metrics, x => x.iou)).toFixed(3)}, <0.80=${metrics.filter(x => x.iou < 0.8).length}, <0.70=${metrics.filter(x => x.iou < 0.7).length}`);
console.log(`boundary q90 median=${median(sorted(metrics, x => x.q90)).toFixed(2)}m, Hausdorff median=${median(sorted(metrics, x => x.hausdorff)).toFixed(2)}m`);
console.log(`algorithm/current matched area=${matched.reduce((s, id) => s + area(algById.get(id).ring), 0).toFixed(0)}/${matched.reduce((s, id) => s + area(curById.get(id).ring), 0).toFixed(0)}m2`);

const selectedAlgorithm = [...algById.values()];
for (const [name, features] of [['initial', initial], ['algorithm', selectedAlgorithm], ['current', current]]) {
  const pairs = overlaps(features);
  console.log(`${name} overlap pairs >0.05m2=${pairs.length}${pairs.length ? ': ' + pairs.map(p => `${p.a}x${p.b}=${p.overlap.toFixed(2)}`).join(', ') : ''}`);
}

const blocks = extractLogBlocks(logPath, ids);
if (logPath && existsSync(logPath) && fullInitial.length) {
  const logText = readFileSync(logPath, 'utf8');
  console.log('\nInput/output identity map:');
  for (const id of ids) {
    const full = fullInitById.get(id);
    if (!full) continue;
    const skipped = logText.includes(`[Output] skip empty geometry for fid=${full.record}`);
    console.log(`id=${id} input_ogr_fid=${full.record} algorithm_present=${algById.has(id) ? 1 : 0} skip_empty=${skipped ? 1 : 0}`);
  }
}
if (blocks.size) {
  console.log('\nLog decisions:');
  for (const id of ids) {
    const lines = blocks.get(id);
    if (!lines) continue;
    console.log(`--- id=${id}`);
    for (const line of lines) console.log(line);
  }
}
