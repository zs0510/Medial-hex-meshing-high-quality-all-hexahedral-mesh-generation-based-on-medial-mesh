#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Skeletal Mesh/SkeletalMesh.h>
#include <Kernel/Mesh.h>
#include <Kernel/IO.h>

#include "DomainProcessor.h"
#include "Remeshing_2.h"
#include "Constrained_Delaunay_Triangulation_2.h"
#include "HexMeshing_FaceSkeleton.h"
#include "HexMeshing_CurveSkeleton.h"

using namespace std;

class HexMeshing_based_on_SkelMesh {

public:
	HexMeshing_based_on_SkelMesh() {};
	void hexmeshing(SkeletalMesh* _skeletal_mesh, MeshKernel::HexMesh& _hex_mesh);

private:
	SkeletalMesh* sklMesh = nullptr;// ����
	MeshKernel::HexMesh hexMesh;// ���

	unordered_map<VH, CH> sklvh_to_hexch;// �Ǽ�����VH ���Ӧ�� ����������CH

	
	vector<vector<VH>> skl_faces_bdy;// �Ǽ��������Ƭ��
	vector<Eigen::Matrix4d> mat_affine_3d_to_2d;// ����Ƭʹ�õķ������
	unordered_map<VH, Eigen::Vector2d> sklvh_to_vex2;// �Ǽ�����VH ���ά����Ķ�Ӧ��ϵ
	 


#pragma region function

	void hexmeshing_face();

	void hexmeshing_curve();

	void detect_triangles();



#pragma endregion

};
