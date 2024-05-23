#pragma once
#include <unordered_map>
#include <vector>
#include <Tools/Math_PCA.h>
#include "Tools/BVH.h"
#include "Skeletal Mesh/SkeletalMesh.h"
#include "Kernel/Mesh.h"
#include "Kernel/IO.h"
#include "Subdivision_Edge.h"
#include <glm/gtx/quaternion.hpp>

// 为中轴网格生成其匹配的六面体网格
// 
// 任务1: 处理 线面分支 节点
//        主要是将与多个曲线骨架相交的面细分, 应保证输出六面体的每个面只与一条骨架线相交
// 
// 任务2: 处理 线线分支 节点
//        在每个线线分支节点上放置合适数目、合适朝向的六面体, 应保证输出六面体的每个面只与一条骨架线相交
// 
// 任务3: 处理 线面分支 节点
//        两个在不同体且相邻的分支面, 中间应该插入一个面作为缓冲地带, 这样在后面的拟合可以取得更好的效果

using namespace std;

class HexMeshing_Branch_Processor {

public:
	HexMeshing_Branch_Processor(SkeletalMesh& skl_mesh, MeshKernel::VolumeMesh& hex_mesh)
		:sklmesh(skl_mesh), hexmesh(hex_mesh) {}

	void optimize(VH& vh, bool only_place__not_split = false);

private:
	SkeletalMesh& sklmesh;
	MeshKernel::VolumeMesh& hexmesh;

	// 分支结点
	unordered_map<EH, bool> curve_skel_ehs;// 纯曲线骨架
	unordered_map<VH, bool> branch_curvh;
	BVH_Tree* bvh_tree = nullptr;
	vector<Ray> rays;
	unordered_map<VH, vector<int>> sklvh_to_ray;// 记录每个顶点和与它相关的射线
	unordered_map<FH, vector<int>> hexfh_to_ray;// 记录每个面与哪些射线相交
	vector<FH> rayid_to_hexfh;// 一条光线只记录最近的那一个面
	vector<Vex> rayid_to_intersection;// 一条光线只有一个交点

	double radius_scale = 0.75;
	
	void init_nodes_map();

	void optimize_curve_branch();// 在每个线线分支结点上放置六面体

	void optimize_face_branch();// 检测每个被多条射线相交的面, 并贪心地细分, 必须保证每个面只被一条射线相交

	void rebuild_bvh_tree();// 使用当前的六面体重新生成 BVH 树
	void recast_rays();// 重新更新所有射线的相交信息
	void recast_ray(int rayid);// 重新更新该条射线的相交信息
	double dist_point_line_segment(Vex& v, Vex& v0, Vex& v1, Vex& nearest_vertex);// 点到直线的最近距离

	

	vector<FH> check_unique_intersect();

	void place_branch_hex(VH vh, vector<pair<VH, Vec>>& adj_vecs);

	vector<Eigen::Vector3d> optimize_orientation(Eigen::Vector3d vec_pca, Eigen::Vector3d vec_base, vector<Eigen::Vector3d>& dirs);

	vector<Eigen::Vector3d> optimize_orientation_zmy(Eigen::Vector3d vec_pca, vector<Eigen::Vector3d>& dirs);

	bool is_corner(VH sklvh);
	glm::quat inAngleAxis(glm::vec3 RotationAxis, double RotationAngle);
	void HexMeshing_Branch_Processor::rotateByQuat(const glm::quat& q, Eigen::Vector3d& in);

};
