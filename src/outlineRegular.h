/*
 * @Descripttion: LOD轮廓规则化
 * @version:
 * @Author: JiangTao
 * @Date: 2025/08/24, 16:45:54
 * @LastEditors: JiangTao
 * @LastEditTime: 2025/09/24, 16:46:49
 */
// =============================================================================
// outlineRegular.h
// 作用：声明"建筑物轮廓规则化"的核心类 outlineRegular，以及形状检测器 ShapeDetector。
//       outlineRegular 的目标：把不规则的建筑底图多边形纠正成主要由正交(平行/垂直)边
//       组成的规整轮廓。内部综合运用：多边形假设生成、能量最小化(data/model/regularity/
//       support/dlg)、Ceres 优化、CGAL 轮廓规则化、以及各种拓扑修复(短边/尖刺/自交等)。
//       典型用法：构造 outlineRegular(原始多边形, 支撑点云[, DLG先验]) 后调用 regular_Contour()，
//       规则化结果存于成员 final_points。
// =============================================================================

#pragma once
#include "MyCloud.h"
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/centroid.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include"modularFunction.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

//using namespace std;
//using Kernel_r = CGAL::Exact_predicates_inexact_constructions_kernel;
//using FT = typename Kernel_r::FT;
//using Point_2 = typename Kernel_r::Point_2;
//using Contour = std::vector<Point_2>;
//using Contour_directions = CGAL::Shape_regularization::Contours::Multiple_directions_2<Kernel_r, Contour>;

using std::string;
using std::vector;

//用于存储单条边的优化信息
struct EdgeRegularInfo {
	int id;
	int type;             // 0: 平行于基准, 1: 垂直于基准, -1: 自由边
	double theta_offset;  // 相对基准角度的偏移 (0, PI/2, PI, etc.)
	double d;             // 待优化的距离参数 (Hessian Normal Form: x*cos - y*sin - d = 0)
	double current_theta; // 当前角度（用于自由边或初始计算）
	std::vector<pcl::PointXYZ> associated_points; // 属于这条边的 RANSAC 内点
	std::vector<double> associated_weights;        // 与 associated_points 一一对应的支撑点可信度
};

// 定义识别结果结构体
struct ShapeResult {
	enum Type { NONE, CIRCLE, ELLIPSE };
	Type type;
	double rmse;          // 拟合误差
	double axis_a;        // 长轴
	double axis_b;        // 短轴
	double center_x;
	double center_y;
	double angle_deg;     // 旋转角度
	double ratio;         // 长短轴之比 (a/b)
};

// outlineRegular：轮廓规则化器。输入原始多边形 + RANSAC 支撑点云(+可选 DLG 先验)，
// 调用 regular_Contour() 完成规则化，结果存于 final_points。
class outlineRegular
{

public:
	struct SortableElement {
		XG::modularFunction::Points_dp point;  // 直接复用已有结构体
		double score;                         // 排序用的second值
		int source;                           // 0=sp_p1, 1=sp_p2（sp_p1优先）
		size_t index;                         // 原向量中的索引

		// 自定义比较规则（严格遵循需求）
		bool operator<(const SortableElement& other) const {
			// 规则1：score大的优先（降序）
			if (score != other.score) {
				return score > other.score;
			}
			// 规则2：score相同时，sp_p1（source=0）优先于sp_p2（source=1）
			if (source != other.source) {
				return source < other.source;
			}
			// 规则3：来源和score都相同时，索引小的优先
			return index < other.index;
		}
	};

	// 辅助结构体：记录顶点信息及其对形状的贡献（面积）
	struct VertexNode {
		int id;               // 原始索引
		double area;          // 该点与其左右邻居构成的三角形面积
		int prev, next;       // 双向链表索引，维持多边形拓扑
		bool is_deleted;      // 标记是否已被删除

		// 最小堆比较算子
		bool operator>(const VertexNode& other) const {
			return area > other.area;
		}
	};

public:
	struct ReferenceWallLine {
		double normal_theta = 0.0;
		double d = 0.0;
		double length = 0.0;
	};

	outlineRegular(vector<pcl::PointXYZ> p, pcl::PointCloud<pcl::PointXYZ>::Ptr ransac_inner_cloud);
	outlineRegular(vector<pcl::PointXYZ> p,
		pcl::PointCloud<pcl::PointXYZ>::Ptr ransac_inner_cloud,
		const std::vector<double>& support_weights);
	outlineRegular(vector<pcl::PointXYZ> p,
		pcl::PointCloud<pcl::PointXYZ>::Ptr ransac_inner_cloud,
		const std::vector<pcl::PointXYZ>& dlg_polygon);
	outlineRegular();
	~outlineRegular();

	void regular_Contour();
	void setSupportDirectionHint(double angle, double peakRatio, std::size_t pairCount);
	void setSourceFeatureId(long long fid) { source_feature_id_ = fid; }
// Mask-only 模式关闭曲线恢复：无正射证据仲裁时，
// 曲线检测器会把栅格楼梯拟合成伪曲线。
	void setCurveRestorationEnabled(bool enabled) { curve_restoration_enabled_ = enabled; }
// 打印结构感知假设修复的运行级汇总。
	static void PrintHypothesisRepairSummary();
	static bool estimateSupportDirection2D(
		const pcl::PointCloud<pcl::PointXYZ>::Ptr& support,
		double& angle,
		double& peakRatio,
		std::size_t& pairCount);
	static bool estimatePcaDirection2D(
		const pcl::PointCloud<pcl::PointXYZ>::Ptr& points,
		double& angle,
		double& axisRatio);
	const std::vector<pcl::PointXYZ>& getBestEnergyHypothesis() const {
		return best_energy_hypothesis_;
	}
	void setInterFloorRegularizationContext(
		const std::vector<double>& building_line_angles,
		const std::vector<ReferenceWallLine>& reference_walls,
		double wall_snap_weight = 0.25,
		double max_wall_snap_distance = 1.2,
		double max_wall_snap_angle_deg = 12.0);
	double GraphZ;//当前规则化底图的z值
	pcl::PointCloud<pcl::PointXYZ>::Ptr ransac_inner_cloud;//一个图中所有线段的ransac拟合内点组成的点云，用于垂直优化过程中计算距离误差
	std::vector<double> support_weights_;

	vector<pcl::PointXYZ> original_points;//纠正前的点
	pcl::PointCloud<pcl::PointXYZ>::Ptr regularzed_points;//纠正后的点
	pcl::PointCloud<pcl::PointXYZ>::Ptr final_points;//平行纠正后的点
	std::vector<pcl::PointXYZ> regularizeBuildingFootprint(std::vector<pcl::PointXYZ> input_points, float angle_remove_threshold = 3.0f,
		float angle_adjust_threshold = 30.0f, float distance_threshold = 10.0f);
	void removeCollinearPoints(std::vector<pcl::PointXYZ>& polygon, double threshold_deg);
	static void regularize_cgal(const std::vector<pcl::PointXYZ>& input_vertices,  // 输入顶点
			//const pcl::PointCloud<pcl::PointXYZ>::Ptr& point_cloud,  // 点云，用于偏差控制/后处理
			std::vector<pcl::PointXYZ>& output_vertices,  // 输出规则化顶点
			double min_length = 3,  // 最小边长（基于点云密度）
			double max_angle_deg = 10.0,  // 最大角度偏差（度）
			double max_offset = 2  // 最大偏移（基于点云分辨率）
		);
	void VDPEnergySimplify(const std::vector<pcl::PointXYZ>& input_vertices,std::vector<pcl::PointXYZ>& output_vertices,double& resolution);
	static void pruneShortEdges(const std::vector<pcl::PointXYZ>& input, std::vector<pcl::PointXYZ>& output, double min_edge_length);
private:

	void generatePolygonalHypotheses(const std::vector<pcl::PointXYZ>& original_points, std::vector<std::vector<pcl::PointXYZ>>& hypotheses);//生成所有假设的多边形
	void generateHypothesesByAreaSubtraction(const std::vector<pcl::PointXYZ>& original_points,std::vector<std::vector<pcl::PointXYZ>>& hypotheses);
	void vertexDrivenDouglasPeucker(const std::vector<XG::modularFunction::Points_dp>& points, double& epsilon,
		std::vector<XG::modularFunction::Points_dp>& out, size_t target_vertices);//通过道格拉斯普克生成指定顶点数的多边形
	double computePointToSegmentDistance(const pcl::PointXYZ& point, const pcl::PointXYZ& seg_start, const pcl::PointXYZ& seg_end);//计算点到线距离
	double computeDataEnergy(const std::vector<pcl::PointXYZ>& hypothesis, const std::vector<pcl::PointXYZ>& original_points);//计算数据能量
	double computeTotalEnergy(const std::vector<pcl::PointXYZ>& hypothesis, const std::vector<pcl::PointXYZ>& original_points, const double lambda, bool verbose = true);//计算总能量
	bool refineThinFeatureHypothesis(std::vector<pcl::PointXYZ>& hypothesis,
		const std::vector<pcl::PointXYZ>& original_points,
		double lambda,
		double resolution);
	double computeOBBArea(const std::vector<pcl::PointXYZ>& points);
	void computeConstraintPairs(const std::vector<pcl::PointXYZ>& poly,std::vector<std::pair<int, int>>& perp, std::vector<std::pair<int, int>>& parallelPairs);
	void ceresOptimize(const std::vector<pcl::PointXYZ>& poly,const pcl::PointCloud<pcl::PointXYZ>::Ptr& points
		,const std::vector<std::pair<int, int>>& perp, const std::vector<std::pair<int, int>>& parallelPairs);
	void updateFinalPointsFromLines();

	std::vector<Eigen::Vector3d> optimized_params; // 优化后的直线参数 [a,b,c]
	double t_threshold = 20.0 * M_PI / 180.0; // 垂直角度阈值（弧度）
	double p_threshold = 20.0 * M_PI / 180.0; // 平行角度阈值（弧度）

	// ransac_inner_cloud的二维包围盒（已扩展20单位），用于最终结果越界检测
	double ransac_min_x, ransac_max_x, ransac_min_y, ransac_max_y;
	//double lambda_g = 100.0; // 约束权重

	
	void resolveShortEdgeIntersections(std::vector<pcl::PointXYZ>& pts, double short_edge_threshold);
	vector<double> data_energy_vec, model_energy_vec;

	void saveAllHypotheses(const std::vector<std::vector<pcl::PointXYZ>>& hypotheses);

	std::vector<XG::modularFunction::Points_dp> mergeAndSort(
		const std::vector<std::pair<XG::modularFunction::Points_dp, double>>& sp_p1,
		const std::vector<std::pair<XG::modularFunction::Points_dp, double>>& sp_p2);

	std::vector<pcl::PointXYZ> best_hypothesis; // 最优假设
	std::vector<pcl::PointXYZ> best_energy_hypothesis_; // VDP+energy 最小化后的原始候选
	std::vector<pcl::PointXYZ> dlg_polygon_;
	double dlg_confidence_ = 0.0;
	bool use_dlg_direction_ = false;
	bool use_dlg_position_ = false;
	void optimizeWithHardConstraints(const std::vector<pcl::PointXYZ>& hypothesis_raw,
		const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
		bool allow_diagonal_edges,
		const std::vector<double>& preferred_line_angles = std::vector<double>());
	double computeDLGPriorEnergy(const std::vector<pcl::PointXYZ>& hypothesis);
	void initializeDLGPrior();

	std::vector<double> building_line_angles_;
	std::vector<ReferenceWallLine> reference_walls_;
	double inter_floor_wall_snap_weight_ = 0.0;
	double inter_floor_max_wall_snap_distance_ = 1.2;
	double inter_floor_max_wall_snap_angle_ = 12.0 * M_PI / 180.0;
	bool has_support_direction_hint_ = false;
	double support_direction_hint_ = 0.0;
	double support_direction_peak_ratio_ = 0.0;
	std::size_t support_direction_pair_count_ = 0;
	long long source_feature_id_ = -1;   // for [HypothesisRepair] logging only
	bool curve_restoration_enabled_ = true;

	double computeAdaptiveLambda(double resolution, const std::vector<pcl::PointXYZ>& points);

	// 基于加权直方图计算全局主方向
	double calculateGlobalDominantAngle(const std::vector<pcl::PointXYZ>& hypothesis);

	

	void OptimizeFinal_points(pcl::PointCloud<pcl::PointXYZ>::Ptr& input);

	// 处理短斜边连接形成的异常（不区分角度类型）
	void resolveOrthogonalSpikes(std::vector<pcl::PointXYZ>& pts, double max_spike_length);

	void cgalContourRegularize(std::vector<pcl::PointXYZ>& pts, double distance_tolerance = 0.5);

	void regularizeBuildingOutline(std::vector<pcl::PointXYZ>& pts, double tolerance);

	void resolveParallelStep(std::vector<pcl::PointXYZ>& pts, double step_threshold);

	// 处理由两条连续短边造成的直角拐角破缺（将 b, c, d 替换为 ab 和 de 的交点）
	void resolveDoubleShortEdgeSpike(std::vector<pcl::PointXYZ>& pts, double threshold);

	// 通用短边折叠：短边夹在两条长边之间时直接延长求交，不判断角度类型
 	void resolveGenericShortEdge(std::vector<pcl::PointXYZ>& pts, double threshold);

	void resolveShortEdgeCluster(std::vector<pcl::PointXYZ>& pts, double threshold);

	// 基于角度和局部形态清理尖角/倒角/缺口，补足短边阈值覆盖不到的错误
	void resolveAngleAnomalies(std::vector<pcl::PointXYZ>& pts, double threshold);

	// 批量评分排序+冲突解决框架，替代原有的贪心顺序修复
	void repairTopologyBatch(std::vector<pcl::PointXYZ>& pts, double threshold);

private:
	pcl::PointXYZ computeIntersection(const Eigen::Vector3d& l1, const Eigen::Vector3d& l2, float z_value);
	// 安全获取得分：越界返回极小值
	double safeGetScore(const std::vector<double>& scores, size_t idx) const {
		return (idx < scores.size()) ? scores[idx] : -1e18;
	}

	// 提取单层中间元素并追加到结果
	void appendLayerMidElements(
		const std::vector<XG::modularFunction::Points_dp>& layer,
		std::vector<XG::modularFunction::Points_dp>& result) const {
		if (layer.size() >= 2) {
			result.insert(result.end(), layer.begin() + 1, layer.end() - 1);
		}
	}

	// 这里保证vector<std::vector<XG::modularFunction::Points_dp>>的一个标准形态
	void  juged_sytle(vector<std::vector<XG::modularFunction::Points_dp>>& points)
	{
		// 遍历每个子向量，确保不越界
		for (int i = 0; i < static_cast<int>(points.size()); ++i) {
			int target_size = i + 3;
			auto& current_vec = points[i];

			// 1. 长度超过目标：截断到目标长度
			if (current_vec.size() > target_size) {
				current_vec.resize(target_size); // 直接resize更高效，替代while循环
				continue;
			}

			// 2. 长度等于目标：无需处理
			if (current_vec.size() == target_size) {
				continue;
			}

			// 3. 长度不足：从后续子向量找可补充的元素（核心修复部分）
			int deficit = target_size - static_cast<int>(current_vec.size());
			int iter = 1;
			int supply_idx = -1;

			// 找后续有多余元素的子向量（避免死循环：限制iter范围）
			while (i + iter < static_cast<int>(points.size())) {
				int supply_target = (i + iter) + 3;
				// 只有后续子向量长度超过其自身目标长度时，才有多余元素可借
				if (points[i + iter].size() > supply_target) {
					supply_idx = i + iter;
					break;
				}
				++iter;
			}

			// 容错：如果没找到可补充的子向量，用默认构造的元素填充（避免崩溃）
			if (supply_idx == -1) {
				current_vec.resize(target_size); // 用默认构造的Points_dp填充
				continue;
			}

			// 从找到的子向量中取元素补充（确保不越界）
			auto& supply_vec = points[supply_idx];
			int supply_start = static_cast<int>(current_vec.size());
			// 计算可取的最大元素数：不超过赤字，且不超过供给向量的多余元素数
			int supply_target = supply_idx + 3;
			int max_take = std::min(deficit, static_cast<int>(supply_vec.size()) - supply_target);

			// 补充元素（优先取供给向量的多余部分）
			if (max_take > 0) {
				current_vec.insert(
					current_vec.end(),
					supply_vec.begin() + supply_target, // 从供给向量的多余部分开始取
					supply_vec.begin() + supply_target + max_take
				);
				// 供给向量同步删除借出的元素（保持其长度为目标长度）
				supply_vec.erase(
					supply_vec.begin() + supply_target,
					supply_vec.begin() + supply_target + max_take
				);
			}

			// 若仍有缺口，用默认构造的元素填充（兜底）
			if (current_vec.size() < target_size) {
				current_vec.resize(target_size);
			}
		}
	}
	 
	// 1. 数据关联：将 RANSAC 点云分配给具体的每一条边
	void assignPointsToEdges(const std::vector<pcl::PointXYZ>& hypothesis,
		const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
		std::vector<EdgeRegularInfo>& edge_infos);

	// 2. 几何分析：确定哪些边是平行组，哪些是垂直组
	void classifyEdgeConstraints(const std::vector<pcl::PointXYZ>& hypothesis,
		std::vector<EdgeRegularInfo>& edge_infos,
		double& initial_base_theta);

	// 3. 硬约束优化：替代原有的 ceresOptimize
	void optimizeWithHardConstraints(std::vector<EdgeRegularInfo>& edge_infos,
		double& base_theta);

	// 4. 重建：根据优化后的 theta 和 d 计算新顶点
	std::vector<pcl::PointXYZ> reconstructPolygon(const std::vector<EdgeRegularInfo>& edge_infos,
		double base_theta);

	// 辅助计算交点
	pcl::PointXYZ computeIntersection(double theta1, double d1, double theta2, double d2);

};

class ShapeDetector {
public:
	/**
	 * @brief 判断点云形状
	 * @param cloud 输入点云 (ransac_inner_cloud)
	 * @param rmse_threshold RMSE阈值，建议设为分辨率的 1.5~2.0 倍 (例如 0.1m)
	 * @param circle_tolerance 圆的容差率，默认 0.1 (即长短轴差异在 10% 以内视为圆)
	 */
	static ShapeResult analyzeCloud(const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
		double rmse_threshold = 0.1,
		double circle_tolerance = 0.1)
	{
		ShapeResult result = { ShapeResult::NONE, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

		if (cloud->points.size() < 6) {
			std::cerr << "[ShapeDetector] 点数太少，无法拟合。" << std::endl;
			return result;
		}

		// 1. 数据准备：构建设计矩阵 D
		// 方程: Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
		int n = cloud->points.size();
		Eigen::MatrixXd D(n, 6);

		// 为了数值稳定性，先去中心化
		Eigen::Vector4f centroid;
		pcl::compute3DCentroid(*cloud, centroid);
		double cx = centroid[0];
		double cy = centroid[1];

		for (int i = 0; i < n; ++i) {
			double x = cloud->points[i].x - cx;
			double y = cloud->points[i].y - cy;
			D.row(i) << x * x, x* y, y* y, x, y, 1.0;
		}

		// 2. 求解最小二乘解 (SVD分解)
		// 求解 D * params = 0，解是 D^T*D 最小特征值对应的特征向量
		// 或者直接对 D 进行 SVD，取右奇异矩阵 V 的最后一列
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(D, Eigen::ComputeThinV);
		Eigen::VectorXd params = svd.matrixV().col(5);

		double A = params(0);
		double B = params(1);
		double C = params(2);
		double D_param = params(3);
		double E = params(4);
		double F = params(5);

		// 3. 判别式检查 (B^2 - 4AC)
		// 椭圆必须满足 B^2 - 4AC < 0
		double det = B * B - 4 * A * C;
		if (det >= 0) {
			// std::cout << "拟合结果是双曲线或抛物线，非闭合形状。" << std::endl;
			return result;
		}

		// 4. 将代数参数转换为几何参数 (中心, 长短轴, 角度)
		// 参考公式: MathWorld - Ellipse
		double q = 64 * (F * (4 * A * C - B * B) - A * E * E + B * D_param * E - C * D_param * D_param) / pow(4 * A * C - B * B, 2);
		double s = 0.25 * sqrt(pow(A - C, 2) + B * B);

		// 几何中心（相对于去中心化坐标系）
		double dx = (2 * C * D_param - B * E) / (B * B - 4 * A * C);
		double dy = (2 * A * E - B * D_param) / (B * B - 4 * A * C);

		// 恢复绝对坐标
		result.center_x = dx + cx;
		result.center_y = dy + cy;

		// 计算长短轴
		// 这里的公式比较复杂，使用特征值法简化理解：
		// 椭圆矩阵 M = [A B/2; B/2 C] 的特征值 lambda1, lambda2
		// 轴长对应 1/sqrt(lambda) (归一化后)

		// 简化计算长短轴的方法：
		double tmp1 = (A * E * E + C * D_param * D_param - B * D_param * E + (B * B - 4 * A * C) * F);
		double tmp2 = (A + C);
		double tmp3 = sqrt((A - C) * (A - C) + B * B);
		double term1 = -sqrt(2 * tmp1 * (tmp2 + tmp3)) / (B * B - 4 * A * C);
		double term2 = -sqrt(2 * tmp1 * (tmp2 - tmp3)) / (B * B - 4 * A * C);

		result.axis_a = std::max(term1, term2);
		result.axis_b = std::min(term1, term2);

		// 计算旋转角
		if (B != 0) {
			result.angle_deg = 0.5 * std::atan2(B, A - C) * 180.0 / M_PI;
		}
		else {
			result.angle_deg = (A < C) ? 0 : 90;
		}

		// 5. 计算 RMSE (评估拟合好坏)
		// 由于计算点到椭圆的几何距离很复杂，这里使用代数距离的近似或采样计算
		// 简单方法：计算点到椭圆方程值的残差，归一化梯度
		double total_error = 0.0;
		for (int i = 0; i < n; ++i) {
			double x = cloud->points[i].x - result.center_x; // 这里要用绝对坐标减绝对中心
			double y = cloud->points[i].y - result.center_y;

			// 旋转回无角度坐标系计算距离偏差 (更准确的几何近似)
			double rad = -result.angle_deg * M_PI / 180.0;
			double xr = x * cos(rad) - y * sin(rad);
			double yr = x * sin(rad) + y * cos(rad);

			// 椭圆方程 (x/a)^2 + (y/b)^2 = 1
			// 距离约为 |r - r_ellipse|
			double r_point = sqrt(xr * xr + yr * yr);
			double angle = atan2(yr, xr);
			double r_ellipse = (result.axis_a * result.axis_b) /
				sqrt(pow(result.axis_b * cos(angle), 2) + pow(result.axis_a * sin(angle), 2));

			total_error += pow(r_point - r_ellipse, 2);
		}
		result.rmse = sqrt(total_error / n);

		// 6. 最终判定
		if (result.rmse > rmse_threshold) {
			result.type = ShapeResult::NONE; // 误差太大，不是椭圆
			// std::cout << "RMSE过大: " << result.rmse << std::endl;
		}
		else {
			result.ratio = result.axis_a / result.axis_b;
			if (std::abs(result.ratio - 1.0) < circle_tolerance) {
				result.type = ShapeResult::CIRCLE;
			}
			else {
				result.type = ShapeResult::ELLIPSE;
			}
		}

		return result;
	}
};
