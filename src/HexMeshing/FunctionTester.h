#pragma once
#include "Kernel/Mesh.h"
#include "Tools/KDTree6d.h"
#include "Tools/MeshMath.h"
#include "Skeletal Mesh/SkeletalMesh.h"

namespace FunctionTester {

	void test_kdtree6d(MeshKernel::SurfaceMesh& trimesh, std::vector<int>& fhs);

	void test_pca(SkeletalMesh& sklmesh, std::vector<std::vector<Vex>>& pca_res);

	void flip_same_orientation_faces(MeshKernel::SurfaceMesh& mesh, FH fh_src, vector<int>& fhs);

};