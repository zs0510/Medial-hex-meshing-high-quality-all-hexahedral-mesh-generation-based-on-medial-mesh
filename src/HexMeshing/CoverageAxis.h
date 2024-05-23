#pragma once
#include <unordered_set>
#include <queue>
#include <omp.h>
#include <gurobi_c++.h>
#include "Kernel/Mesh.h"
#include "Tools/VoronoiDiagram.h"
#include "Tools/AABB_Tree.h"
#include "Tools/BVH.h"
#include "Skeletal Mesh/SkeletalMesh.h"
#include "Skeletal Mesh/VoronoiMesh.h"

using namespace std;

// 计算一个三角形网格的覆盖轴网格（中轴）

class CoverageAxis {

	struct InnerPoint;
	struct SurfacePoint;
	struct CoverageNode;

	struct InnerPoint {
		int id = -1;
		Vector3d pos;
		double radius = 0;
		CoverageNode* coverage_node = nullptr;
		InnerPoint() {}
	};

	struct SurfacePoint {
		int id = -1;
		Vector3d pos;
		vector<CoverageNode*> coverage_nodes;
		SurfacePoint() {}
	};

	struct CoverageNode {
		bool is_valid = true;
		InnerPoint* centroid = nullptr;
		unordered_set<SurfacePoint*> coverage_points;
		CoverageNode() {}
	};

	struct CoverageNode_CMP {
		bool operator()(CoverageNode* cn1, CoverageNode* cn2) {
			return cn1->coverage_points.size() < cn2->coverage_points.size();// 大顶堆
		}
	};

public:
	CoverageAxis(MeshKernel::SurfaceMesh& _mesh);

	void compute_coverage_axis();

	void get_inner_points_selected(vector<Vector3d>& inner_points_pos);

	void get_medial_mesh(SkeletalMesh&);

private:
	MeshKernel::SurfaceMesh& mesh;

	bool is_computed = false;

	SkeletalMesh medial_mesh;

	unordered_map<InnerPoint*, VH> ip_to_vh;
	unordered_map<VH, InnerPoint*> vh_to_ip;

	AABB_Tree* ab_tree = nullptr;
	BVH_Tree* bvh_tree = nullptr;

	vector<SurfacePoint*> surface_points;// 输入网格表面上的采样点
	vector<InnerPoint*> inner_points;

	unordered_set<SurfacePoint*> surface_points_coveraged;
	unordered_set<InnerPoint*> inner_points_selected;
	unordered_set<VH> vhs_selected;

	//const double scale_radius = 1.5;
	double offset_radius = 0.02;

	priority_queue<CoverageNode*, vector<CoverageNode*>, CoverageNode_CMP> heap;

	void init();

	void inner_points_generation();

	void point_selection_based_on_set_coverage_heap();

	void point_selection_based_on_set_coverage_MILP();

	void connection_establishment();
	void connection_establishment(vector<bool>& is_selected);// 留下选择的顶点

	void hole_filling();
	bool find_boundary_edges_loop(EH eh, vector<VH>&, vector<EH>&);
	bool neighbor_with_hanging_edge(VH vh);
	bool is_all_hanging(vector<EH> ehs);

	void postprocessing();

};

class CoverageAxis_Vorn {

private:
	struct InnerPoint;
	struct SurfacePoint;
	struct CoverageNode;

	struct InnerPoint {
		int id = -1;
		Vector3d pos;
		double radius = 0;
		CoverageNode* coverage_node = nullptr;
		InnerPoint() {}
	};

	struct SurfacePoint {
		int id = -1;
		Vector3d pos;
		vector<CoverageNode*> coverage_nodes;
		SurfacePoint() {}
	};

	struct CoverageNode {
		bool is_valid = true;
		InnerPoint* centroid = nullptr;
		unordered_set<SurfacePoint*> coverage_points;
		CoverageNode() {}
	};

public:

	CoverageAxis_Vorn(MeshKernel::SurfaceMesh& trimesh, SkeletalMesh& medial_mesh);

	void compute_coverage_axis();

private:
	MeshKernel::SurfaceMesh& mesh;
	SkeletalMesh& medial_mesh;

	AABB_Tree* ab_tree = nullptr;
	BVH_Tree* bvh_tree = nullptr;
	//const double scale_radius = 1.5;
	double offset_radius = 0.02;
	vector<SurfacePoint*> surface_points;// 输入网格表面上的采样点
	vector<InnerPoint*> inner_points;

	void inner_points_generation();
	void point_selection_based_on_set_coverage_MILP();
	void connection_establishment(vector<bool>& is_selected);// 留下选择的顶点

	void postprocessing();

};