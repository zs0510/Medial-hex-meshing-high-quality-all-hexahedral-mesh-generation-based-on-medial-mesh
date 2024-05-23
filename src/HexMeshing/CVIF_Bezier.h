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

class CVIF_Bezier {

public:
	CVIF_Bezier(MeshKernel::VolumeMesh& hexmesh, MeshKernel::SurfaceMesh& target_mesh);

	void getSurfaceMesh(MeshKernel::SurfaceMesh& surfaceMesh) {
		surfaceMesh = surface_mesh;
	}

	void getBezierMesh(MeshKernel::SurfaceMesh& bezierMesh) {
		bezierMesh = bezier_mesh;
	}

	void getBezierSurface(vector<vector<vector<Vector3d>>>& MutiBezierPatch) {
		generate_bezier_mesh(true);
		extract_bezier_surface(MutiBezierPatch);
	}

	void cvif_bezier(vector<vector<vector<Vector3d>>>& MutiBezierPatch);

private:
	MeshKernel::SurfaceMesh surface_mesh;// 从六面体表面提取出来的曲面网格
	MeshKernel::SurfaceMesh refmesh;

	MeshKernel::SurfaceMesh bezier_mesh;// 曲面网格对应的 Bezier 网格
	MeshKernel::SurfaceMesh subdivision_mesh;
	MeshKernel::SurfaceMesh subdivision_tri_mesh;
	unordered_map<FH, FH> trisub2quadsub;

	// 计算 Bezier 曲面
	unordered_map<VH, VH> v_to_corner_point;// 曲面网格顶点 到 Bezier曲面顶点的映射
	unordered_map<FH, vector<VH>> fv_to_inter_point;// 曲面网格面与顶点 到 Bezier曲面顶点的映射
	unordered_map<EH, vector<VH>> ev_to_edge_point;// 曲面网格边与顶点 到 Bezier曲面顶点的映射

	unordered_map<VH, unordered_map<VH, double>> bezier_to_control_weights;// Bezier曲面 VH 与 控制曲面 VH 之间的权值关系
	unordered_map<VH, unordered_map<VH, double>> subdivision_to_bezier_weights;// 细分曲面 VH 与 Bezier曲面 VH 之间的权值关系
	unordered_map<VH, vector<pair<double, Vec>>> bezier_movements;// Bezier曲面的移动向量
	unordered_map<VH, vector<pair<double, Vec>>> control_movements;// 控制曲面的移动向量
	
	MeshKernel::iGameKdTree* kd_tree = nullptr;
	KDTree6d* kdtree6d = nullptr;

	bool use_kdtree6d_flag = false;// 一般情况下不使用 6-D 的 KD树
	bool is_converged = false;
	double epsilon_error = 0.001;
	double rms_error = -1;// root mean square (RMS) fitting error of the kth iteration

	vector<VH> nearest_subvhs;
	vector<Vec> nearest_deltas;

	vector<vector<vector<VH>>> mutil_bezier_patch;

	void generate_surface_mesh(MeshKernel::VolumeMesh& hexmesh);// 提取六面体的表面作为四边形曲面网格

	void generate_bezier_mesh(bool extract_patch = false);// 根据四边形曲面网格计算相应 Bezier 网格的控制顶点

	void init_kdtree();
	void calc_movement();
	void update_surface_mesh();
	void extract_bezier_surface(vector<vector<vector<Vector3d>>>& MutiBezierPatch);

	void calc_corner_point(VH);
	void calc_inter_point(FH, VH);
	void calc_edge_point(EH, VH);

	void calc_movement_face_center();

};