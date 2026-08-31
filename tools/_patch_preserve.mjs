import { readFileSync, writeFileSync } from 'fs';
let src = readFileSync('src/outlineRegular.cpp', 'utf8');

const oldBlock = `        // 传递给下一次迭代
        // 注: 平行相邻边(台阶 jog)求交 NaN 丢顶点、退化合并是优化器清理
        // 台阶伪影的自然机制, 不强行锁顶点数(锁住会保留中间态斜边);
        // 输出层面的方向一致性由拓扑通道的 rogue 检查兜底。
        const bool new_polygon_valid =
            new_polygon.size() >= 3 &&
            isSimplePolygon2D(new_polygon) &&
            polygonIoU2D(new_polygon, backup_polygon) >= (allow_diagonal_edges ? 0.55 : 0.65);

        if (!new_polygon_valid) {
            std::cerr << "[Ceres] reject invalid iteration polygon; keep last valid polygon"
                      << std::endl;
            current_polygon = last_valid_polygon;
            break;
        }

        if (new_polygon.size() == n && new_polygon.size() >= 3) {
            current_polygon = new_polygon;
            last_valid_polygon = current_polygon;
        }
        else if (new_polygon.size() >= 3) {
            // 顶点数变化但多边形仍有效，继续迭代
            current_polygon = new_polygon;
            last_valid_polygon = current_polygon;
        }
        else {
            // 如果求交失败导致拓扑破坏，终止迭代
            break;`;

const newBlock = `        // 传递给下一次迭代
        // preserve_topology: 顶点数必须守恒——近平行求交 NaN 丢点属于
        // 隐式删点(违反"删点只能由约简操作显式完成"的契约),
        // 该轮整体回滚到 last_valid_polygon
        const bool new_polygon_valid =
            new_polygon.size() >= 3 &&
            isSimplePolygon2D(new_polygon) &&
            polygonIoU2D(new_polygon, backup_polygon) >= (allow_diagonal_edges ? 0.55 : 0.65) &&
            (!preserve_topology || new_polygon.size() == n);

        if (!new_polygon_valid) {
            std::cerr << "[Ceres] reject invalid iteration polygon; keep last valid polygon"
                      << std::endl;
            current_polygon = last_valid_polygon;
            break;
        }

        current_polygon = new_polygon;
        last_valid_polygon = current_polygon;`;

if (!src.includes(oldBlock)) { console.log('OLD NOT FOUND'); process.exit(1); }
src = src.replace(oldBlock, newBlock);
writeFileSync('src/outlineRegular.cpp', src);
console.log('preserve_topology enforced');
