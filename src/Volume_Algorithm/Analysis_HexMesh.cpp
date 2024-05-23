#include "Analysis_HexMesh.h"

void Analysis_HexMesh::checkTwoFacesOnBoundary(MeshKernel::VolumeMesh& hexmesh, 
	std::unordered_map<CH, bool>& bad_cells, std::unordered_map<FH, bool>& bad_faces) {

	bad_cells.clear();
	bad_faces.clear();

	for (auto& cp : hexmesh.allcells()) {

		if (!hexmesh.isOnBoundary(cp.first)) continue;

		const auto& fhs = cp.second.getFaceHandle();

		int surface_count = 0;
		for (auto& fh : fhs) {
			if (hexmesh.isOnBoundary(fh)) {
				surface_count++;
			}
		}

		if (surface_count > 1) {
			bad_cells[cp.first] = true;
			for (auto& fh : fhs) {
				bad_faces[fh] = true;
			}
		}

	}

	//std::cout << "There ara " << bad_cells.size() << " boundary elements have more than one facet exposed on the surface.\n";

}