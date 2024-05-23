#pragma once
#include <gurobi_c++.h>
#include "Kernel/Mesh.h"
#include "Kernel/IO.h"
#include "Skeletal Mesh/SkeletalMesh.h"
#include "Project_KD_Tree.h"
#include "Tools/BVH.h"
#include "Volume_Algorithm/VolumeEditer.h"
//#include "Volume_Algorithm/Add_Cell.h"
#include "Subdivision_Edge.h"

using namespace std;

struct CurveLine {
	VH beg_vh, end_vh;// ����ǼܵĶ���
	EH beg_eh, end_eh;// ����Ǽܵı�
	vector<VH> vertices;
	vector<VH> hexvhs_src, hexvhs_end;
	CurveLine() {}
};

struct CurveConstraint {
	EH hexeh_p[2];
	bool is_b2b = true;
	CurveConstraint() {}
};

struct RayConstraint {
	int rid = 0, segment = 1;
	FH hexfh;
	EH hexeh;
	RayConstraint() {}
};

class HexMeshing_ConformingTessellation {

public:
	HexMeshing_ConformingTessellation(SkeletalMesh& _sklmesh, MeshKernel::VolumeMesh& _hexmesh) : sklmesh(_sklmesh), hexmesh(_hexmesh) {
		
	}

	void conforming_split(VH);

	vector<pair<Vex, Vex>> hexeh_pairs;
	vector<pair<FH, FH>> hexfh_pairs;
	vector<vector<Vex>> path_dense_vex;

	unordered_map<FH, vector<vector<FH>>> hexfh_origin_to_current;
	//unordered_map<FH, vector<vector<VH>>> hexfh_origin_to_current_vhs;
	unordered_map<FH, vector<VH>> hexfh_origin_to_sub_order;
	vector<vector<vector<FH>>> rayid_to_sub_fhs;// ��¼�����ÿ�����ߵ�ϸ�����������

	vector<int> debug_should_be_com_faces;
	vector<pair<Vex, Vex>> matched_vex;

private:
	SkeletalMesh& sklmesh;
	MeshKernel::VolumeMesh& hexmesh;
	MeshKernel::VolumeMesh hexmesh_origin;
	
	BVH_Tree* bvh_tree = nullptr;

	unordered_set<VH> branch_sklvh;
	unordered_set<VH> joint_sklvh;
	unordered_set<VH> end_sklvh;

	unordered_set<EH> curve_skel_ehs;
	unordered_map<VH, vector<int>> sklvh_to_rayid;
	unordered_map<VH, unordered_map<EH, int>> sklvh_skleh_to_rayid;

	vector<Ray> rays;
	vector<FH> rayid_to_hexfh;
	vector<Vex> rayid_to_intersection;
	vector<pair<VH, EH>> rayid_to_sklvh_skleh;// ��¼ÿ�����߶�Ӧ������Ǽܵ� �� ��
	

	vector<RayConstraint> ray_constraints;
	vector<int> constraint_equal_rcs;

	unordered_map<FH, vector<int>> hexfh_to_rayid;// ��¼ÿ��������Щ�����ཻ
	
	unordered_map<EH, int> hexeh_to_minsegnum;
	unordered_map<FH, pair<EH, EH>> hexfh_to_hexeh;// ��¼ÿ������������ ��Ҫ����ϸ�ֶ����� ��
	unordered_map<FH, EH> hexfh_to_subeh;// ��¼ÿ������������ ������Ҫϸ�ֵ� ��

	vector<CurveLine> curve_lines_b2b, curve_lines_b2e;
	int constraint_count_b2b = 0;

	unordered_set<EH> visited;

	unordered_map<EH, unordered_set<EH>> constraint_equal_ehs;// һ��һ���Լ��: ��eh Ӧ���� ��constraint1[eh] ��������һ��

	unordered_map<EH, vector<vector<int>>> hexeh_to_rcid;
	unordered_map<EH, unordered_map<int, int>> hexeh_rayid_to_rcid;

	/* ������ */

	void init_nodes_map(VH);

	void rebulid_bvh();
	void recast_rays();

	void calc_segnum();

	void calc_constraint();

	void opt_conforming_split();

	void assign_sweep_face();

	void sweep_curve_skel();

	void cone_detection();

	/* �Ӻ��� */
	void assign_face_dense_to_origin();
	void assign_face_to_ray();

	void recast_ray(int rayid, int sklvh = -1, int skleh = -1);

	bool dfs_curve_vhs_between_branch(std::vector<VH>& vhs);

	bool is_corner(VH vh);

	void min_dist(vector<Vex>, vector<VH>&);
	void min_dir(vector<Vex>, vector<VH>&, Vec);

	FH getFaceHandle(MeshKernel::VolumeMesh& meah, EH, EH);

	FH getFaceHandle(MeshKernel::VolumeMesh& mesh, VH, VH, VH);

};