#pragma once
#include "Kernel/Mesh.h"
#include "Volume_Algorithm/Smoothing.h"
//#include "cinolib/padding.h"

using namespace std;

namespace HexMesh_TopologicalOperators {

	void padding(MeshKernel::VolumeMesh& hexmesh, bool add_inner_vertices = true);

	void dicing(MeshKernel::VolumeMesh& hexmesh, EH input_eh, int num_of_segments);

	void get_dicing_vhs(MeshKernel::VolumeMesh& hexmesh, EH input_eh, std::unordered_map<VH, VH>& left2right);

	bool is_dicing_ok(MeshKernel::VolumeMesh& hexmesh, EH input_eh);

};