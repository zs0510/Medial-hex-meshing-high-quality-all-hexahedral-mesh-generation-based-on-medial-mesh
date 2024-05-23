#pragma once
#include "Kernel/Mesh.h"
#include "Skeletal Mesh/SkeletalMesh.h"
#include "Project_KD_Tree.h"

namespace BoundarySmoother {


	void smoothing_2D(MeshKernel::SurfaceMesh& mesh, int iters = 30, double move_speed = 0.1);

	void smoothing_2D(SkeletalMesh& mesh, std::vector<int>& non_mainfold_vhs, int iters = 30, double move_speed = 0.1);

	//void branch_smoothing(SkeletalMesh& mesh, double move_ratio = 0.5);

	bool distance_to_line(const Vex& p, const Vex& v0, const Vex& v1, double& dist, Vex& np);

	void smoothing_branch(SkeletalMesh& mesh);

};