#pragma once
#include <vector>
#include <utility>
#include <Kernel/Mesh.h>
#include <Skeletal Mesh/SkeletalMesh.h>
#include <Tools/BVH.h>
#include "Volume_Algorithm/Add_Cell.h"
#include "Subdivision_Edge.h"

struct JointPlane {
	double a, b, c, d;
	JointPlane(){};
	JointPlane(double _a, double _b, double _c, double _d):
		a(_a), b(_b), c(_c), d(_d) {
	}
};

class HexMeshing_CurveSkeleton {

public:
	HexMeshing_CurveSkeleton(SkeletalMesh& _sklmesh, MeshKernel::VolumeMesh& _hexmesh): 
		sklmesh(_sklmesh), hexmesh(_hexmesh) {

	}
	void hexmeshing();

private:
	SkeletalMesh& sklmesh;
	MeshKernel::VolumeMesh& hexmesh;

	std::unordered_map<VH, bool> branching_nodes;
	std::unordered_map<VH, bool> joint_nodes;
	std::unordered_map<VH, bool> end_nodes;

	std::unordered_map<EH, bool> curve_skel_ehs;// 骨架的曲线边
	
	std::unordered_map<VH, std::unordered_map<EH, FH>> branchingnode_to_fh;// 骨架VH 到 六面体CH 的映射
	std::unordered_map<VH, FH> jointnode_to_fh;// 骨架VH 到 六面体FH 的映射
	std::unordered_map<VH, FH> endnode_to_fh;// 骨架VH 到 六面体FH 的映射

	std::unordered_map<VH, JointPlane> jointnode_to_jp;// 骨架VH 到 平面的映射
	std::unordered_map<VH, JointPlane> endnode_to_jp;

	double radius_scale = 0.75;

	void init_nodes_map();// 初始化骨架的曲线
	void set_joint_plane();
	void set_end_plane();
	void set_branch_face();
	void lanuch_from_branching_node();
	
	bool dfs_curve_vhs(std::vector<VH>& vhs);
	void extrude_curve_vhs(std::vector<VH>& vhs);
	
	void cone_detection();
	bool is_corner(VH sklvh);

};