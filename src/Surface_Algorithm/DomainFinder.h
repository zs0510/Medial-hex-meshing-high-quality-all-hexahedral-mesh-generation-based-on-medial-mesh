#pragma once
#include "Kernel/Mesh.h"

namespace DomainFinder {

	void get_decomposition_from_quad(MeshKernel::SurfaceMesh& quad_mesh, std::vector<std::pair<Vex, Vex>>& edges);

}