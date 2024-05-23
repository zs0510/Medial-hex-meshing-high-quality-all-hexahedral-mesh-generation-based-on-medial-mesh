#pragma once
#include <set>
#include <omp.h>
#include "Kernel/Mesh.h"
#include "Tools/meshmath.h"
#include "Streamline.h"

class QuadMesher {
public:
	// 输入：三角形网格和特征边
	QuadMesher(MeshKernel::SurfaceMesh _trimesh, std::vector<EH> _feature_ehs): trimesh(_trimesh), feature_ehs(_feature_ehs){};
	
	bool meshing(MeshKernel::SurfaceMesh& quadmesh);

	void get_streamlines(std::vector<Streamline>& res, std::vector<SingularPoint>& points) {// 返回流线
		res = streamlines;
		points = singular_points;
	}

	void get_cross_field(std::vector<std::vector<Vec>>& _cross_field) {
		_cross_field = crossfield;
	}

	void get_streamlines_bad(std::vector<Streamline>& res) {
		res = bad_streamlines;
	}

private:

#pragma region data

	MeshKernel::SurfaceMesh trimesh;// input parameter
	std::vector<EH> feature_ehs;

	std::vector<SingularPoint> singular_points;// 所有的奇异点
	std::vector<Streamline> streamlines;// 所有有效的流线
	std::vector<std::vector<Vec>> crossfield;// 所有顶点的标架信息
	std::vector<EH> boundary_ehs;
	std::unordered_map<EH, int> singular_neighbor_ehs;// eh 到 该边上奇异点id 的映射
	std::unordered_map<FH, int> singular_neighbor_fhs;// fh 到 该面上奇异点id 的映射

	std::vector<Streamline> bad_streamlines;// 无效的流线，调试用

	// 三角形网格 handle 到其对应奇异点 id 的映射
	std::unordered_map<VH, int> vh_to_spid;
	std::unordered_map<EH, int> eh_to_spid;
	std::unordered_map<FH, int> fh_to_spid;

#pragma endregion

#pragma region pipeline
	bool preprocessing();

	void calculate_cross_field();// 计算标架场 by Su
	void generate_initial_streamlines();// 生成流线
	void simplify_streamlines();// 流线的简化
	void quad_meshing();

#pragma endregion

	void init_singular_neighbor_handles();

	/* 新方法生成流线 */
	void lengthen_initial_streamlines(const MarchNode&, Streamline&, std::unordered_map<FH, bool>& visited, int bdy_face_cnt = 0);// 不使用上个面片的标架信息

	void merge_neighbor_edgehandle();// 简化前的预处理，将相邻的终点边合并

	bool get_first_point_streamline(Vec v, EH eh, FH fh, std::vector<Vec>& march_dir, std::vector<double>& weight_vh1);// 检测是否有从奇异点 v 出发射向 eh 的流线
	std::vector<Vec> get_cross_dirs(const Vec&, EH eh, FH fh);// 返回该点在边上插值得到的标架
	std::vector<Vec> get_cross_dirs(VH vh, FH fh);// 返回该点在该三角形面片上的标架
	std::vector<Vec> get_cross_dirs(const Vec&);// 通过输入的一个标架场分量生成四个完整的标架
	bool get_intersection_line2(const Vex& source, const Vec& dir, EH eh, Vex& intersection);// 返回二维射线与线段 eh 的交点, true 表示相交, false 表示未相交
	Vec get_most_similar_dir(const Vec& ref_dir, const std::vector<Vec>& dirs);// 返回 dirs 中与 ref_dir 方向最接近的方向向量
	std::unordered_set<EH> get_smallest_bounding_ehs(const EH& eh);// 返回这条边的最小多边形的边
	std::vector<std::pair<EH, FH>> get_samllest_bounding_ehs_fhs(const EH& eh);// 返回包围这条边的最小多边形的边以及面
	void get_streamline_itrp(const Streamline& sl0, const Streamline& sl1, Streamline& sl_new);// 由 sl0 和 sl1 插值得到新流线, 注意: sl0 须与 sl1 同向 
	void get_among_ehs(const EH& beg, const EH& end, std::set<EH>& res);// 返回连接 边beg 和 边end 之间的最小个数边界边集合(含 beg 和 end)
	void generate_initial_streamlines_v(const SingularPoint& point);// 生成流线
	void generate_initial_streamlines_e(const SingularPoint& point);// 生成流线
	void generate_initial_streamlines_f(const SingularPoint& point);// 生成流线
	
};