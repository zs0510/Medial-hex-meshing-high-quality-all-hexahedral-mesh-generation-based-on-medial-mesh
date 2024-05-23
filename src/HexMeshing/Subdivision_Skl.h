#pragma once
#include "Skeletal Mesh/SkeletalMesh.h"

class Subdivision_Skl {
public:
	Subdivision_Skl(SkeletalMesh* _mesh) : mesh(_mesh) {};
	void subdivision(int iter);

private:
	SkeletalMesh* mesh;

};

void Subdivision_Skl::subdivision(int iter) {

	for (int it = 0; it < iter; ++it) {

		auto faces_original = mesh->allfaces();
		std::vector<EH> face_edges;
		std::unordered_map<EH, VH> eh2vh;
		for (auto& ep : mesh->alledges()) {
			if (!mesh->NeighborFh(ep.first).empty()) {
				auto& e = mesh->edges(ep.first);
				auto ev = mesh->getEdgeMidpoint(ep.first);
				auto evh = mesh->AddVertex(ev);
				eh2vh[ep.first] = evh;
				mesh->radius[evh] = (mesh->radius[e.vh1()] + mesh->radius[e.vh2()]) * 0.5;
				face_edges.push_back(ep.first);
			}
		}

		for (auto& fp : faces_original) {
			auto& face = mesh->faces(fp.first);
			auto ehs = face.getEdgeHandle();
			auto vhs = face.getVertexHandle();
			assert(ehs.size() == 3 && vhs.size() == 3);
			mesh->AddFace({ VH(eh2vh[ehs[0]]), VH(eh2vh[ehs[1]]), VH(eh2vh[ehs[2]]) });
			for (int i = 0; i < 3; ++i) {
				auto cur = vhs[i];
				auto pre = vhs[(i + 2) % 3];
				auto next = vhs[(i + 1) % 3];
				VH new_pre(-1), new_next(-1);
				for (auto eh : ehs) {
					const auto& e = mesh->edges(eh);
					if (e.vh1() != cur && e.vh2() != cur) continue;
					int vex = e.vh1() + e.vh2() - cur;
					if (new_pre == -1 && vex == pre) {
						new_pre = eh2vh[eh];
					}
					if (new_next == -1 && vex == next) {
						new_next = eh2vh[eh];
					}
				}
				assert(new_pre != -1 && new_next != -1);
				mesh->AddFace({ new_pre, cur, new_next });
			}
			//std::cout << "Subdivision Skl: #F " << fp.first << " subdivision success.\n";
		}

		for (auto& eh : face_edges) {
			if (mesh->isValid(eh)) {
				mesh->DeleteEdge(eh);
			}
		}

	}

	mesh->updateAllHandles();

}