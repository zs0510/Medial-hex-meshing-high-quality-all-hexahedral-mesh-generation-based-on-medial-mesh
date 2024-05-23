#pragma once
#include "SkeletalMesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <queue>
#include <vector>

struct V_SEQ {
	Eigen::Matrix4d A = Eigen::Matrix4d::Zero();
	Eigen::Vector4d b = Eigen::Vector4d::Zero();
	double c = 0;
	V_SEQ() {};
};

struct Skl_CollapseNode {
	bool is_valid = true;
	EH eh = MeshKernel::iGameEdgeHandle(-1);
	Eigen::Vector4d sphere = Eigen::Vector4d::Zero();
	double cost = std::numeric_limits<double>::max();
	Skl_CollapseNode() {};
};

struct Skl_CollapseNode_CMP {
	bool operator() (Skl_CollapseNode* n1, Skl_CollapseNode* n2) {
		return n1->cost > n2->cost;
	}
};


class Q_MAT {
public:

	Q_MAT(SkeletalMesh& _mesh) : mesh(_mesh){

	}
	void simplification(int target_vertices);
	

private:

	SkeletalMesh& mesh;
	std::priority_queue<Skl_CollapseNode*, std::vector<Skl_CollapseNode*>, Skl_CollapseNode_CMP> pri_que;
	int check_topology_threshold = 200;
	double hard_k = 0.00001 * 1000 * 1000;// hard_k Խ���Ƴ�ë���ٶ�Խ�죬����ֲ�Խ�����ȣ�Ĭ��ֵ 1E-5F when bounding box diagonal length = 1
	double soft_k = 0.1;// soft_k Խ�󣬱����߽�Ч��Խ�ã���ֵӦ�� 0 - 1 ֮��
	std::vector<V_SEQ> v_sqe;// ����ÿ���������ڵ���ƽ��
	std::vector<std::vector<Skl_CollapseNode*>> v_collapse;// ����ÿ��������ص�����������

	void calc_total_sqe();// Slab Quadratic Error
	void calc_face_sqe(FH);
	void calc_edge_sqe(EH);
	void calc_vertex_sqe(VH);

	void calc_collapse_cost(EH);

	double get_edge_stability(const MeshKernel::iGameEdge& edge);
	double get_radius_by_interpolation(const MeshKernel::iGameVertex&, VH , VH);

	bool get_triangle_from_three_spheres(Eigen::Vector3d c0, double r0, Eigen::Vector3d c1, double r1, Eigen::Vector3d c2, double r2, 
		Eigen::Vector3d& n0, Eigen::Vector3d& n1);
	bool distance_to_line(const Eigen::Vector3d& p, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, double& dist, Eigen::Vector3d& fp);
	
};