#pragma once
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <time.h>
#include "Kernel/Mesh.h"
#include "Surface_Algorithm/QuadMeshSubdivision_CC.h"
#include "Volume_Algorithm/Smoothing.h"
#include "Volume_Algorithm/Subdivision_Linear.h"
#include "Tools/KdTree.h"
#include "Tools/KDTree6d.h"

using namespace std;

class CVIF_Bezier_Unsmoothed {

public:
	CVIF_Bezier_Unsmoothed(MeshKernel::VolumeMesh& hexmesh, MeshKernel::SurfaceMesh& target_mesh);

	void getControlMesh(MeshKernel::SurfaceMesh& ctrlMesh) {
		ctrlMesh = control_mesh;
	}

	void cvif_bezier(vector<vector<vector<Vector3d>>>& MutiBezierPatch);

private:
	MeshKernel::SurfaceMesh control_mesh;
	MeshKernel::SurfaceMesh refmesh;

	MeshKernel::SurfaceMesh subdivision_mesh;
	MeshKernel::SurfaceMesh subdivision_tri_mesh;
	unordered_map<VH, unordered_map<VH, double>> control_weights;// 细分曲面 VH 与 控制曲面 VH 之间的权值关系
	unordered_map<FH, FH> trisub2quadsub;

	MeshKernel::iGameKdTree* kd_tree = nullptr;
	KDTree6d* kdtree6d = nullptr;

	unordered_map<VH, vector<pair<double, Vec>>> movements;// 控制曲面的移动向量

	int ctrl_subdivision_iter = 1;
	bool use_kdtree6d_flag = false;
	bool is_converged = false;
	double epsilon_error = 0.001;
	double rms_error = -1;// root mean square (RMS) fitting error of the kth iteration

	vector<VH> nearest_subvhs;
	vector<Vec> nearest_deltas;

	vector<vector<vector<VH>>> mutil_bezier_patch;

	void generate_control_mesh(MeshKernel::VolumeMesh& hexmesh);

	void sub_to_bezier_ctrlmesh();

	void init_kdtree();
	void calc_movement();
	void update_control_mesh();
	void extract_bezier_surface(vector<vector<vector<Vector3d>>>& MutiBezierPatch);

};
