#include "DomainFinder.h"

void DomainFinder::get_decomposition_from_quad(MeshKernel::SurfaceMesh& quad_mesh, std::vector<std::pair<Vex, Vex>>& edges) {

	std::unordered_set<EH> lines;
	for (auto& ep : quad_mesh.alledges()) {
		if (quad_mesh.isOnBoundary(ep.first)) {
			lines.insert(ep.first);
		}
	}

	std::unordered_set<VH> singular_points;
	for (auto& vp : quad_mesh.allvertices()) {
		if (quad_mesh.isOnBoundary(vp.first)) {
			if (quad_mesh.NeighborVh(vp.first).size() == 4) {
				singular_points.insert(vp.first);
			}
		} else if (quad_mesh.NeighborVh(vp.first).size() != 4) {
			singular_points.insert(vp.first);
		}
	}

	for (auto& svh : singular_points) {
		for (auto& beg_eh : quad_mesh.NeighborEh(svh)) {
			VH vh_pre = svh;
			EH eh = beg_eh;
			while (!lines.count(eh)) {
				lines.insert(eh);
				auto& edge = quad_mesh.edges(eh);
				VH vh_next(edge.vh1() + edge.vh2() - vh_pre);
				if (singular_points.count(vh_next) || quad_mesh.isOnBoundary(vh_next)) break;
				// 选择下一条边
				EH eh_next(-1);
				std::unordered_set<EH> ehs_in_same_faces;
				for (auto& adjfh : quad_mesh.NeighborFh(eh)) {
					auto& face = quad_mesh.faces(adjfh);
					auto fehs = face.getEdgeHandle();
					ehs_in_same_faces.insert(fehs.begin(), fehs.end());
				}
				for (auto& adjeh : quad_mesh.NeighborEh(vh_next)) {
					if (ehs_in_same_faces.count(adjeh) || lines.count(adjeh)) {
						continue;
					}
					eh_next = adjeh;
					break;
				}
				if (eh_next != -1) {
					eh = eh_next;
					vh_pre = vh_next;
				}
			}
		}
	}

	for (auto& eh : lines) {

		auto& edge = quad_mesh.edges(eh);
		edges.emplace_back(quad_mesh.vertices(edge.vh1()), quad_mesh.vertices(edge.vh2()));

	}

}