#pragma once
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <time.h>
#include "Kernel/IO.h"
#include "Tools/AABB_Tree.h"
#include "Tools/KdTree.h"
#include "Tools/KDTree6d.h"
#include "Surface_Algorithm/QuadMeshSubdivision_CC.h"
#include "Volume_Algorithm/Smoothing.h"
#include "Volume_Algorithm/Subdivision_Linear.h"

// Auther: Zhang Sheng
// Date: 2022/7/20
// Reference: Quality guaranteed all-hex mesh generation by a constrained volume iterative fitting algorithm_lin2015

// Improve: add feature-preserving at 2023/6/18 by Zhang Sheng 

using namespace std;

class VolumeFitting_CVIF {
public:
	VolumeFitting_CVIF(MeshKernel::VolumeMesh& _hexmesh, MeshKernel::SurfaceMesh& _refmesh): 
		hexmesh(_hexmesh), refmesh(_refmesh) {
		refmesh.genAllVerticesNormal();
	}
	void volume_fitting(int subdivision_times = 1, 
		bool project_to_surface_flag = false, 
		int iter_times = 100, 
		unordered_map<VH, Vex> corner_features_hex = {},
		unordered_map<VH, int> line_features_hex = {},
		vector<vector<EH>> line_features_ref = {}
		);

private:
	MeshKernel::VolumeMesh& hexmesh;
	MeshKernel::SurfaceMesh& refmesh;
	MeshKernel::SurfaceMesh control_mesh;
	MeshKernel::SurfaceMesh subdivision_mesh;
	MeshKernel::SurfaceMesh subdivision_tri_mesh;

	unordered_map<VH, VH> hex2ctrl;// 输入六面体 VH 与其控制网格 VH 的映射关系
	unordered_map<VH, VH> ctrl2hex;
	unordered_map<FH, FH> trisub2quadsub;

	unordered_map<VH, unordered_map<VH, double>> control_weights;// 细分曲面 VH 与 控制曲面 VH 之间的权值关系

	MeshKernel::iGameKdTree* kd_tree = nullptr;
	KDTree6d* kdtree6d = nullptr;

	AABB_Tree* aabb_tree = nullptr;

	unordered_map<VH, vector<pair<double, Vec>>> movements;// 控制曲面的移动向量

	int ctrl_subdivision_iter = 1;
	bool use_kdtree6d_flag = false;
	bool is_converged = false;
	double epsilon_error = 0.001;
	double rms_error = -1;// root mean square (RMS) fitting error of the kth iteration

	unordered_map<VH, Vex> corner_features_hex;
	unordered_map<VH, int> line_features_hex;
	vector<vector<EH>> line_features_ref;

	bool use_refmesh_face = true;
	vector<Vex> refmesh_points;
	vector<Vector3d> refmesh_normals;
	vector<VH> nearest_subvhs;
	vector<Vec> nearest_deltas;

	void init_refmesh_data();
	void init_control_mesh();
	void init_kdtree();
	void calc_movement();
	void update_control_mesh();
	void update_hexmesh();

	void project_to_surface();

	void feature_recovery();

	void calc_movement_face_center();
	double dist_point_line_segment(iGameVertex& v, iGameVertex& v0, iGameVertex& v1, iGameVertex& nearest_vertex);

};

class Coarse_CVIF {

public:
	Coarse_CVIF(MeshKernel::VolumeMesh& _hexmesh, MeshKernel::SurfaceMesh& _refmesh) : hexmesh(_hexmesh), trimesh(_refmesh) {

	}

	void cvif_coarse();

private:

	MeshKernel::VolumeMesh& hexmesh;
	MeshKernel::SurfaceMesh& trimesh;

};