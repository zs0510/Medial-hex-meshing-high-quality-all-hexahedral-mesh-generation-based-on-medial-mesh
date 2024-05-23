#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Skeletal Mesh/SkeletalMesh.h>

class CAD_Reconstruction {

public:

	CAD_Reconstruction() {};

	void get_bezier_control_mesh(MeshKernel::VolumeMesh& _hexmesh, vector<vector<vector<Vector3d>>>& MutiBezierPatch);

	void get_bezier_control_mesh(MeshKernel::SurfaceMesh& _quadmesh, vector<vector<vector<Vector3d>>>& MutiBezierPatch);

	void get_bezier_control_mesh(MeshKernel::SurfaceMesh& _quadmesh, vector<vector<vector<Vector3d>>>& MutiBezierPatch, 
		MeshKernel::SurfaceMesh& _refmesh, unordered_set<EH>& feature_ehs, unordered_set<VH>& feature_vhs);

	void solve_aux_surface_mesh(MeshKernel::VolumeMesh& _hexmesh, MeshKernel::SurfaceMesh& result_mesh);
	void solve_aux_surface_mesh(MeshKernel::SurfaceMesh& _quadmesh, MeshKernel::SurfaceMesh& result_mesh);

private:
	MeshKernel::SurfaceMesh surface_mesh;
	MeshKernel::SurfaceMesh bezier_mesh;// ���������Ӧ�� Bezier ����
	MeshKernel::SurfaceMesh aux_surface_mesh;// ���� aux_surface_mesh ����� Bezier �������񶥵�ֱ�������� surface_mesh ��
	// ���� Bezier ����
	unordered_map<VH, VH> v_to_corner_point;// �������񶥵� �� Bezier���涥���ӳ��
	unordered_map<FH, vector<VH>> fv_to_inter_point;// �����������붥�� �� Bezier���涥���ӳ��
	unordered_map<EH, vector<VH>> ev_to_edge_point;// ����������붥�� �� Bezier���涥���ӳ��
	vector<vector<vector<VH>>> mutil_bezier_patch;

	vector<Eigen::Triplet<double>> triplets;// �� aux_surface_mesh ��
	vector<double> corner_weights_sum;

	void get_bezier_control_mesh(vector<vector<vector<Vector3d>>>& MutiBezierPatch);
	void get_bezier_control_mesh(vector<vector<vector<Vector3d>>>& MutiBezierPatch, 
		MeshKernel::SurfaceMesh& _refmesh, unordered_set<EH>& feature_ehs, unordered_set<VH>& feature_vhs);

	void generate_surface_mesh(MeshKernel::VolumeMesh& hexmesh);
	void calc_corner_point(VH);
	void calc_inter_point(FH, VH);
	void calc_edge_point(EH, VH);

	Vex calc_corner_point_aux(VH);
	double dist_point_line_segment(MeshKernel::iGameVertex& v, MeshKernel::iGameVertex& v0, MeshKernel::iGameVertex& v1, 
		MeshKernel::iGameVertex& nearest_vertex);

};