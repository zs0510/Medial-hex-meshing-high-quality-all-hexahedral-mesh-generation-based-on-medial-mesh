#pragma once
#include <QtWidgets/qfiledialog.h>
#include <QtCore/QtCore>
#include <QtWidgets/QLayout>
#include <QtWidgets/QMessageBox>
#include <vector>
#include <queue>
#include <utility>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "Skeletal Mesh/SkeletalMesh.h"
#include "Surface_Algorithm/Parameterization_ARAP.h"
#include "Kernel/IO.h"
#include "Remeshing_2.h"
#include "Tools/MeshMath.h"

using namespace std;

class DomainProcessor {
public:

	DomainProcessor(){}

	void detect_faces(SkeletalMesh& mesh, vector<vector<VH>>& domains);// 返回若干圈边界边
	void get_input_domains(SkeletalMesh& mesh, vector<vector<VH>>& domains, vector<vector<Vex>>& expanded_domains);// 将边界边扩张
	void remeshing_3d(SkeletalMesh& mesh, vector<vector<VH>>& domains_vhs, vector<vector<Vex>>& expanded_domains, vector<MeshKernel::SurfaceMesh>& trimeshs);

	// 改进的解法: 仍不支持非流形
	bool is_manifold(SkeletalMesh& mesh);
	void expand_boundary(SkeletalMesh& mesh, vector<vector<Vex>>& expanded_boundary);
	void decomposition(SkeletalMesh& mesh);

private:
	
	double radius_scale = 0.75;

};

