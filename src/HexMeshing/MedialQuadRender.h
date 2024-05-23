#pragma once
#include <QtWidgets/qinputdialog.h>
#include <omp.h>
#include <pmp/algorithms/DistancePointTriangle.h>
#include <pmp/algorithms/BarycentricCoordinates.h>
#include "Kernel/Mesh.h"
#include "Kernel/IO.h"
#include "Skeletal Mesh/SkeletalMesh.h"
#include "Project_KD_Tree.h"
#include "Tools/BVH.h"
#include "Tools/AABB_Tree.h"

namespace MedialQuadRender {

	using namespace std;

	void get_streamlines(SkeletalMesh& skel_mesh, vector<vector<Vex>>& res);

}