#pragma once
#include "Kernel/Mesh.h"

#ifndef D_PI
#define D_PI 6.2831852f
#endif // !D_PI

namespace Analysis_HexMesh {

	void checkTwoFacesOnBoundary(MeshKernel::VolumeMesh& hexmesh, std::unordered_map<CH, bool>& bad_cells, std::unordered_map<FH, bool>& bad_faces);

};
