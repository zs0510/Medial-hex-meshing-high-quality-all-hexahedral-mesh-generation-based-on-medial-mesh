#pragma once
#include "Kernel/Mesh.h"
#include "Kernel/IO.h"

class SingularGraph {
public:
	SingularGraph(MeshKernel::VolumeMesh& _hexmesh): hexmesh(_hexmesh) {

	}

	void get_singular_branches(std::vector<std::vector<EH>>& singular_branches);

private:

	MeshKernel::VolumeMesh& hexmesh;
	std::unordered_set<VH> singular_vhs;
	std::unordered_set<EH> singular_ehs;
	std::unordered_set<FH> boundary_fhs;

};