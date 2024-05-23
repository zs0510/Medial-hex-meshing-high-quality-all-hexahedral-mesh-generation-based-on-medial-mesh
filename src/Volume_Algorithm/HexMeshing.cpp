#include "HexMeshing.h"

void HexMeshing::hexmeshing_based_on_skeleton() {


}

void HexMeshing::classify_skl_vertices() {
	if (!sklMesh) return;

	

}

void HexMeshing::resolution_control() {



}

void HexMeshing::find_partition_boundary(std::vector<vector<EH>>& ehs, std::vector<FH>& fhs) {

	if (sklMesh == nullptr) return;

	std::unordered_map<EH, bool> is_boundary;
	for (auto& ep : sklMesh->alledges()) {
		if (sklMesh->NeighborFh(ep.first).size() == 1) {
			is_boundary[ep.first] = true;
			//ehs.push_back({ep.first});
		}
	}

	// find and remove non-mainfold edges

	ehs.clear();
	for (auto& ep : is_boundary) {
		if (!ep.second) continue;// has been visited
		
		std::vector<EH> loop;
		std::unordered_map<EH, bool> visited;
		
		bool is_loop = find_partition_boundary_recursive(ep.first, loop, visited, is_boundary);
		
		if (is_loop) {
			ehs.push_back(loop);
			for (auto& eh : loop) {
				is_boundary[eh] = false;// markd has been visited
			}
		}
	}



	std::cout << "Find partition boundary success. Domain size = " << ehs.size() << endl;

}

bool HexMeshing::find_partition_boundary_recursive(EH eh, std::vector<EH>& loop, std::unordered_map<EH, bool>& stored, 
	std::unordered_map<EH, bool>& is_on_boundary) {

	if (!loop.empty() && eh == loop[0]) return true;
	if (stored.count(eh)) return false;// ÓÐÐ¡Ñ­»·

	loop.push_back(eh);
	stored[eh] = true;

	for (auto& adjeh : sklMesh->NeighborEh(eh)) {
		if (loop.size() >= 2 && adjeh == loop[loop.size() - 2]) continue;
		if (is_on_boundary.count(adjeh) && find_partition_boundary_recursive(adjeh, loop, stored, is_on_boundary)) {
			return true;
		}
	}

	loop.pop_back();
	stored.erase(eh);
	return false;
}