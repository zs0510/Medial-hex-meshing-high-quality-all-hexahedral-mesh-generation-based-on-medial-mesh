#pragma once
#include <QtWidgets/qinputdialog.h>
#include <omp.h>
#include <pmp/algorithms/DistancePointTriangle.h>
#include <pmp/algorithms/BarycentricCoordinates.h>
#include "Kernel/Mesh.h"
#include "Kernel/IO.h"
#include "Skeletal Mesh/SkeletalMesh.h"
#include "Tools/KdTreeMedialMesh.h"
#include "Tools/BVH.h"
#include "Tools/AABB_Tree.h"

class HexMeshing_FaceSkeleton_Nonmanifold {// 通过直接读入四边形网格来生成面上的六面体
public:
	HexMeshing_FaceSkeleton_Nonmanifold(SkeletalMesh& _sklmesh, MeshKernel::VolumeMesh& _hexmesh): hexmesh(_hexmesh), sklmesh(_sklmesh) {}

	void hexmeshing();

private:
	MeshKernel::VolumeMesh& hexmesh;
	SkeletalMesh& sklmesh;
	double radius_scale = 0.85;
	void normal_filtering(MeshKernel::SurfaceMesh& quad_mesh, vector<Vec>& normals);

};
