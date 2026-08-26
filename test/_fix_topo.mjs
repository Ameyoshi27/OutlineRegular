// 重写 TopologyPreservingRegularize 函数
import { readFileSync, writeFileSync } from 'fs';
const lines = readFileSync('src/outlineRegular.cpp', 'utf8').split('\r\n');
// 找到函数起止(4125行到4393行, 0基4124-4392)
const start = 4124;
const end = 4392;
// 确认起止
console.log('start:', lines[start].slice(0,50));
console.log('end:', lines[end].slice(0,50));
// 替换
const newFn = readFileSync('test/_topo_v2.txt', 'utf8').split('\r\n');
lines.splice(start, end - start + 1, ...newFn);
writeFileSync('src/outlineRegular.cpp', lines.join('\r\n'), 'utf8');
console.log('replaced', end - start + 1, 'old lines with', newFn.length, 'new lines');
