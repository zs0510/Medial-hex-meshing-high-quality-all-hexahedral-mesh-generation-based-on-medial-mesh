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

// Ϊ��������������ƥ�������������
// 
// ����1: ���� �����֧ �ڵ�
//        ��Ҫ�ǽ��������߹Ǽ��ཻ����ϸ��, Ӧ��֤����������ÿ����ֻ��һ���Ǽ����ཻ
// 
// ����2: ���� ���߷�֧ �ڵ�
//        ��ÿ�����߷�֧�ڵ��Ϸ��ú�����Ŀ�����ʳ����������, Ӧ��֤����������ÿ����ֻ��һ���Ǽ����ཻ
// 
// ����3: ���� �����֧ �ڵ�
//        �����ڲ�ͬ�������ڵķ�֧��, �м�Ӧ�ò���һ������Ϊ����ش�, �����ں������Ͽ���ȡ�ø��õ�Ч��

using namespace std;

class HexMeshing_Branch_Processor {

public:
	HexMeshing_Branch_Processor(SkeletalMesh& skl_mesh, MeshKernel::VolumeMesh& hex_mesh)
		:sklmesh(skl_mesh), hexmesh(hex_mesh) {}

	void optimize(VH& vh, bool only_place__not_split = false);

private:
	SkeletalMesh& sklmesh;
	MeshKernel::VolumeMesh& hexmesh;

	// ��֧���
	unordered_map<EH, bool> curve_skel_ehs;// �����߹Ǽ�
	unordered_map<VH, bool> branch_curvh;
	BVH_Tree* bvh_tree = nullptr;
	vector<Ray> rays;
	unordered_map<VH, vector<int>> sklvh_to_ray;// ��¼ÿ�������������ص�����
	unordered_map<FH, vector<int>> hexfh_to_ray;// ��¼ÿ��������Щ�����ཻ
	vector<FH> rayid_to_hexfh;// һ������ֻ��¼�������һ����
	vector<Vex> rayid_to_intersection;// һ������ֻ��һ������

	double radius_scale = 0.75;
	
	void init_nodes_map();

	void optimize_curve_branch();// ��ÿ�����߷�֧����Ϸ���������

	void optimize_face_branch();// ���ÿ�������������ཻ����, ��̰�ĵ�ϸ��, ���뱣֤ÿ����ֻ��һ�������ཻ

	void rebuild_bvh_tree();// ʹ�õ�ǰ���������������� BVH ��
	void recast_rays();// ���¸����������ߵ��ཻ��Ϣ
	void recast_ray(int rayid);// ���¸��¸������ߵ��ཻ��Ϣ
	double dist_point_line_segment(Vex& v, Vex& v0, Vex& v1, Vex& nearest_vertex);// �㵽ֱ�ߵ��������

	

	vector<FH> check_unique_intersect();

	void place_branch_hex(VH vh, vector<pair<VH, Vec>>& adj_vecs);

	vector<Eigen::Vector3d> optimize_orientation(Eigen::Vector3d vec_pca, Eigen::Vector3d vec_base, vector<Eigen::Vector3d>& dirs);

	vector<Eigen::Vector3d> optimize_orientation_zmy(Eigen::Vector3d vec_pca, vector<Eigen::Vector3d>& dirs);

	bool is_corner(VH sklvh);
	glm::quat inAngleAxis(glm::vec3 RotationAxis, double RotationAngle);
	void HexMeshing_Branch_Processor::rotateByQuat(const glm::quat& q, Eigen::Vector3d& in);

};
