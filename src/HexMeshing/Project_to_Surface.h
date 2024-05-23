#pragma once
#include <omp.h>
#include "Kernel/Mesh.h"
#include "Kernel/IO.h"
#include "Tools/AABB_Tree.h"
#include "Tools/BVH.h"
#include "Volume_Algorithm/Smoothing.h"

class Project_to_Surface {
public:
	Project_to_Surface(MeshKernel::SurfaceMesh& _ref_mesh, MeshKernel::VolumeMesh& _hex_mesh)
		: refmesh(_ref_mesh), hexmesh(_hex_mesh) {

	}

	void project_to_nearest_point();

private:
	MeshKernel::SurfaceMesh& refmesh;
	MeshKernel::VolumeMesh& hexmesh;
	AABB_Tree* ab_tree = nullptr;

	void project_exterior();
	
};
