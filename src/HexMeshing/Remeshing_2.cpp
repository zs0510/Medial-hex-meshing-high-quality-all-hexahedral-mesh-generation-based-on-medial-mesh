#include "Remeshing_2.h"

void Remeshing_2::execute_igame(double ratio) {

	double len = 0;
	mesh.genAllEdgesLength();
	for (auto& ep : mesh.alledges()) {
		auto& e = mesh.edges(ep.first);
		len += e.getLength();
	}

	len /= mesh.esize();
	len *= ratio;

	for (int i = 0; i < 10; ++i) {

		split_long_edges(len);

		collapse_short_edges(len);

		equalize_valence();

		relaxation();

	}

	mesh.updateAllHandles();
	//relaxation(10);

}

void Remeshing_2::split_long_edges(double target_len) {

	double upper_length = target_len * 4 / 3;

	int split_cnt = 0;
	auto oldEhs = mesh.alledges();
	for (auto& ep : oldEhs) {
		if (!mesh.isValid(ep.first)) continue;
		auto& e = mesh.edges(ep.first);
		auto vh1 = ep.second.vh1();
		auto vh2 = ep.second.vh2();
		auto& v1 = mesh.vertices(vh1);
		auto& v2 = mesh.vertices(vh2);
		double len = (v1 - v2).norm();
		if (len < upper_length) continue;
		MeshKernel::iGameVertex new_v = (v1 + v2) / 2;
		auto adjfhs = mesh.NeighborFh(ep.first);
		auto new_vh = mesh.AddVertex(new_v);
		for (auto& fh : adjfhs) {// adjacent triangles are bisected
			auto newFace = mesh.faces(fh).getVertexHandle();
			if (newFace[0] == vh1) newFace[0] = new_vh;
			else if (newFace[1] == vh1) newFace[1] = new_vh;
			else if (newFace[2] == vh1) newFace[2] = new_vh;
			mesh.AddFace(newFace);
			newFace = mesh.faces(fh).getVertexHandle();
			if (newFace[0] == vh2) newFace[0] = new_vh;
			else if (newFace[1] == vh2) newFace[1] = new_vh;
			else if (newFace[2] == vh2) newFace[2] = new_vh;
			mesh.AddFace(newFace);
		}
		mesh.DeleteEdge(ep.first);
		split_cnt++;

	}

}

void Remeshing_2::collapse_short_edges(double target_len) {

	double lower_length = target_len * 0.8;
	double upper_length = target_len * 4 / 3;

	int collapse_cnt = 0;
	auto oldEhs = mesh.alledges();
	for (auto& ep : oldEhs) {
		if (!mesh.isValid(ep.first) || mesh.isOnBoundary(ep.first) || mesh.isOnBoundary(ep.second.vh1()) || mesh.isOnBoundary(ep.second.vh2())) continue;
		auto& e = mesh.edges(ep.first);
		auto vh1 = ep.second.vh1();
		auto vh2 = ep.second.vh2();
		auto& v1 = mesh.vertices(vh1);
		auto& v2 = mesh.vertices(vh2);
		double len = (v1 - v2).norm();
		if (len > lower_length) continue;
		auto vh_keep = vh1;
		auto vh_remove = vh2;
		auto& vex_keep = mesh.vertices(vh_keep);
		auto& vex_remove = mesh.vertices(vh_remove);
		// preserve boundary
		if (mesh.isOnBoundary(vh_keep)) {
			// do nothing
		} else if (mesh.isOnBoundary(vh_remove)) {
			vex_keep = vex_remove;
		} else {
			vex_keep = (vex_keep + vex_remove) / 2;
		}

		std::vector<std::vector<MeshKernel::iGameVertexHandle>> new_triangles;
		for (auto& fh : mesh.NeighborFh(vh_remove)) {
			auto vex = (mesh.faces(fh)).getVertexHandle();
			assert(vex.size() > 2);
			if (std::find(vex.begin(), vex.end(), vh_keep) != vex.end()) continue;
			for (int i = 0; i < vex.size(); ++i) {
				if (vex[i] == vh_remove) {
					vex[i] = vh_keep;
					new_triangles.push_back(vex);
					break;
				}
			}
		}

		mesh.DeleteVertex(vh_remove);

		// Ìí¼ÓÈý½ÇÐÎ
		for (auto& tri : new_triangles) {
			mesh.AddFace(tri);
		}

		collapse_cnt++;

	}

}

void Remeshing_2::equalize_valence() {

	int flip_cnt = 0;
	auto oldEdges = mesh.alledges();
	for (auto& ep : oldEdges) {
		auto eh = ep.first;
		if (!mesh.isValid(eh) || mesh.isOnBoundary(eh)) continue;
		auto& e = mesh.edges(eh);

		//if (!is_flip_ok(eh)) continue;
		auto vh1 = e.vh1();
		auto vh2 = e.vh2();
		MeshKernel::iGameVertexHandle vh3(-1), vh4(-1);
		auto adjfhs = mesh.NeighborFh(eh);
		//if (adjfhs.size() != 2 || !mesh.isValid(*adjfhs.begin()) || !mesh.isValid(*(adjfhs.begin()++))) continue;
		//assert(adjfhs.size() == 2 && !mesh.isValid(*adjfhs.begin()) && !mesh.isValid(*(adjfhs.begin()++)));
		for (auto fh : adjfhs) {
			auto vex = (mesh.faces(fh)).getVertexHandle();
			assert(vex.size() > 2);
			int tmp = vex[0] + vex[1] + vex[2] - vh1 - vh2;
			if (vh3 == -1) vh3 = MeshKernel::iGameVertexHandle(tmp);
			else vh4 = MeshKernel::iGameVertexHandle(tmp);
		}
		//if (vh3 == vh4 || vh3 == -1 || vh4 == -1) continue;
		int target_degree1 = (mesh.isOnBoundary(vh1)) ? 4 : 6;
		int target_degree2 = (mesh.isOnBoundary(vh2)) ? 4 : 6;
		int target_degree3 = (mesh.isOnBoundary(vh3)) ? 4 : 6;
		int target_degree4 = (mesh.isOnBoundary(vh4)) ? 4 : 6;
		int diff1 = mesh.NeighborVh(vh1).size() - target_degree1;
		int diff2 = mesh.NeighborVh(vh2).size() - target_degree2;
		int diff3 = mesh.NeighborVh(vh3).size() - target_degree3;
		int diff4 = mesh.NeighborVh(vh4).size() - target_degree4;
		int old_diff = diff1 * diff1 + diff2 * diff2 + diff3 * diff3 + diff4 * diff4;
		diff1--, diff2--;
		diff3++, diff4++;
		int new_diff = diff1 * diff1 + diff2 * diff2 + diff3 * diff3 + diff4 * diff4;
		if (new_diff < old_diff) {
			flip_cnt++;
			flip(eh);
		}
	}

}

void Remeshing_2::flip(EH eh) {
	if (!mesh.isValid(eh)) return;
	auto& e = mesh.edges(eh);
	auto vh1 = e.vh1();
	auto vh2 = e.vh2();
	MeshKernel::iGameVertexHandle vh3(-1), vh4(-1);
	auto adjfhs = mesh.NeighborFh(eh);
	//if (adjfhs.size() != 2 || !mesh.isValid(*adjfhs.begin()) || !mesh.isValid(*(adjfhs.begin()++))) continue;
	//assert(adjfhs.size() == 2 && !mesh.isValid(*adjfhs.begin()) && !mesh.isValid(*(adjfhs.begin()++)));
	for (auto fh : adjfhs) {
		auto vex = (mesh.faces(fh)).getVertexHandle();
		assert(vex.size() > 2);
		int tmp = vex[0] + vex[1] + vex[2] - vh1 - vh2;
		if (vh3 == -1) vh3 = MeshKernel::iGameVertexHandle(tmp);
		else vh4 = MeshKernel::iGameVertexHandle(tmp);
	}
	for (auto& fh : adjfhs) {
		auto vhs = (mesh.faces(fh)).getVertexHandle();
		int i;
		for (i = 0; i < 3; ++i) {
			int j = (i + 1) % 3;
			if ((vhs[i] == vh1 && vhs[j] == vh2) || (vhs[i] == vh2 && vhs[j] == vh1)) {
				int tmp = vhs[0] + vhs[1] + vhs[2] - vh1 - vh2;
				if (tmp == vh3) vhs[i] = vh4;
				else vhs[i] = vh3;
				assert(vhs[0] != vhs[1] && vhs[1] != vhs[2] && vhs[2] != vhs[0]);
				break;
			}
		}
		assert(i < 4);
		auto new_fh = mesh.AddFace(vhs);
	}
	mesh.DeleteEdge(eh);
}

void Remeshing_2::relaxation(int iter) {

	for (int it = 0; it < iter; ++it) {

		std::vector<VH> allvhs;
		allvhs.reserve(mesh.vsize());
		for (auto& vp : mesh.allvertices()) {
			allvhs.push_back(vp.first);
		}
		std::vector<Vex> newpos(allvhs.size());

#pragma omp parallel for
		for (int i = 0; i < allvhs.size(); ++i) {
			VH vh = allvhs[i];
			if (mesh.isOnBoundary(vh)) continue;
			newpos[i] = Vex(0, 0, 0);
			const auto& adjvhs = mesh.NeighborVh(vh);
			if (adjvhs.empty()) continue;
			for (auto& adjvh : adjvhs) {
				auto& adjv = mesh.vertices(adjvh);
				newpos[i] += adjv;
			}
			newpos[i] /= adjvhs.size();
		}

#pragma omp parallel for
		for (int i = 0; i < allvhs.size(); ++i) {
			VH vh = allvhs[i];
			if (mesh.isOnBoundary(vh)) continue;
			auto& v = mesh.vertices(vh);
			v.setPosition(newpos[i].x(), newpos[i].y(), newpos[i].z());
		}


	}

}