#include "BoundaryVertexMerger.h"

void BoundaryVertexMerger::merge() {

	std::unordered_map<VH, bool> boundary_vertices;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			boundary_vertices[vp.first] = true;
		}
	}

	double len_avg = 0;
	mesh.genAllEdgesLength();
	for (auto& ep : mesh.alledges()) {
		len_avg += ep.second.getLength();
	}
	len_avg /= mesh.esize();

	for (auto& vp0 : boundary_vertices) {
		if (vp0.second == false) continue;
		vp0.second = false;
		VH vh0 = vp0.first, vh1(-1);
		auto& v0 = mesh.vertices(vh0);
		double dis_min = 999999999;
		for (auto& vp1 : boundary_vertices) {
			if (vp1.second == false) continue;
			auto& v1 = mesh.vertices(vp1.first);
			double dis = (v0 - v1).norm();
			if (dis < dis_min) {
				dis_min = dis;
				vh1 = vp1.first;
				//std::cout << dis << std::endl;
			}
		}
		std::cout << "vh0 = " << vh0 << ", vh1 = " << vh1 << std::endl;
		if (!mesh.isValid(vh1)) {
			std::cerr << "Error: get a invalid vertex handle.\n";
			continue;
		} else if (mesh.isConnected(vh0, vh1)) {
			std::cout << "Warning: get a pair of connected vertex pair.\n";
			continue;
		} else if (dis_min < len_avg) {
			std::cout << "Dis is too big.\n";
			continue;
		} else {
			std::cout << "We merge a pair of vertices. Dis = " << dis_min << "\n" ;
		}
		boundary_vertices[vh1] = false;
		std::vector<std::vector<VH>> newfaces;
		for (auto& adjfh : mesh.NeighborFh(vh1)) {
			auto& adjf =mesh.faces(adjfh);
			auto vhs = adjf.getVertexHandle();
			for (auto& vh : vhs) {
				if (vh == vh1) {
					vh = vh0;
					break;
				}
			}
			newfaces.push_back(vhs);
		}

		mesh.DeleteVertex(vh1);

		for (auto& f : newfaces) {
			mesh.AddFace(f);
		}

	}

}