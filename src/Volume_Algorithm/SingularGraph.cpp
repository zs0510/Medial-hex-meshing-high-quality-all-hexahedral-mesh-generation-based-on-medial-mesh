#include "SingularGraph.h"

void SingularGraph::get_singular_branches(std::vector<std::vector<EH>>& singular_branches) {

	for (auto& fp : hexmesh.allfaces()) {
		if (hexmesh.NeighborCh(fp.first).size() == 1) {
			boundary_fhs.insert(fp.first);
		}
	}

	auto vh_is_on_boundary = [](MeshKernel::VolumeMesh& hexmesh, std::unordered_set<FH>& boundary_fhs, VH vh) {
		for (auto& fh : hexmesh.NeighborFh(vh)) {
			if (boundary_fhs.count(fh)) {
				return true;
			}
		}
		return false;
	};

	auto eh_is_on_boundary = [](MeshKernel::VolumeMesh & hexmesh, std::unordered_set<FH>&boundary_fhs, EH eh) {
		for (auto& fh : hexmesh.NeighborFh(eh)) {
			if (boundary_fhs.count(fh)) {
				return true;
			}
		}
		return false;
	};

	for (auto& vp : hexmesh.allvertices()) {
		int regular = vh_is_on_boundary(hexmesh, boundary_fhs, vp.first) ? 4 : 6;
		if (hexmesh.NeighborEh(vp.first).size() != regular) {
			singular_vhs.insert(vp.first);
		}
	}
	std::cout << "The size of singular_vhs is " << singular_vhs.size() << std::endl;

	for (auto& ep : hexmesh.alledges()) {
		if (eh_is_on_boundary(hexmesh, boundary_fhs, ep.first)) continue;
		auto& edge = hexmesh.edges(ep.first);
		if (singular_vhs.count(edge.vh1()) && singular_vhs.count(edge.vh2())) {
			singular_ehs.insert(ep.first);
		}
	}
	std::cout << "The size of singular_ehs is " << singular_ehs.size() << std::endl;

	// for test
	std::vector<EH> singular_ehs_vec(singular_ehs.begin(), singular_ehs.end());
	singular_branches.push_back(singular_ehs_vec);

}