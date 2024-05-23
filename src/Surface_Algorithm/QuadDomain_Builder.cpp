#include "QuadDomain_Builder.h"

void QuadDomain_Builder::get_quad_domains(vector<int>& ehs) {

	for (auto& vp : singular_vhs) {
		VH vh = vp.first;
		const auto& adjehs = mesh.NeighborEh(vh);
		for (auto& adjeh : adjehs) {
			auto& edge = mesh.edges(adjeh);
			VH adjvh(edge.vh1() + edge.vh2() - vh);
			if (singular_visited[vh].count(adjvh)) continue;
			singular_visited[vh][adjvh] = true;
			vector<EH> quad_edge;
			unordered_map<EH, bool> visited;
			march_inr(adjeh, vh, adjvh, quad_edge, visited);

			quad_edges.push_back(quad_edge);
			quad_edges_srcvh.push_back(vh);
		}
	}

	for (auto& vp : bdy_singular_vhs) {
		VH vh = vp.first;
		const auto& adjehs = mesh.NeighborEh(vh);
		for (auto& adjeh : adjehs) {
			if (!mesh.isOnBoundary(adjeh)) continue;
			auto& edge = mesh.edges(adjeh);
			VH adjvh(edge.vh1() + edge.vh2() - vh);
			if (singular_visited[vh].count(adjvh)) continue;
			singular_visited[vh][adjvh] = true;
			vector<EH> quad_edge;
			unordered_map<EH, bool> visited;
			march_bdy(adjeh, vh, adjvh, quad_edge, visited);
			quad_edges.push_back(quad_edge);
			quad_edges_srcvh.push_back(vh);
		}
	}

	ehs.clear();
	for (auto& quad_edge : quad_edges) {
		for (auto& eh : quad_edge) {
			ehs.push_back(eh);
		}
	}

}

void QuadDomain_Builder::march_inr(EH eh, VH src, VH tgt, vector<EH>& quad_edge, unordered_map<EH, bool>& visited) {

	if (visited.count(eh)) return;

	quad_edge.push_back(eh);
	visited[eh] = true;

	if (singular_vhs.count(tgt)) {
		singular_visited[tgt][src] = true;
		return;
	}

	if (boundary_vhs.count(tgt)) {
		bdy_singular_vhs[tgt] = true;
		return;
	}

	EH next_eh = get_next_eh(eh, src, tgt);

	if (next_eh == -1) return;

	auto& edge = mesh.edges(next_eh);
	src = tgt;
	tgt = VH(edge.vh1() + edge.vh2() - src);

	march_inr(next_eh, src, tgt, quad_edge, visited);

}

void QuadDomain_Builder::march_bdy(EH eh, VH src, VH tgt, vector<EH>& quad_edge, unordered_map<EH, bool>& visited) {

	if (visited.count(eh)) return;
	quad_edge.push_back(eh);
	visited[eh] = true;

	if (bdy_singular_vhs.count(tgt)) {
		singular_visited[tgt][src] = true;
		return;
	}

	const auto& adjehs = mesh.NeighborEh(tgt);
	EH next_eh(-1);
	for (auto& adjeh : adjehs) {
		if (mesh.isOnBoundary(adjeh) && adjeh != eh) {
			next_eh = adjeh;
			break;
		}
	}

	if (next_eh == -1) return;
	auto& edge = mesh.edges(next_eh);
	src = tgt;
	tgt = VH(edge.vh1() + edge.vh2() - src);

	march_bdy(next_eh, src, tgt, quad_edge, visited);

}

EH QuadDomain_Builder::get_next_eh(EH eh, VH src, VH tgt) {

	const auto& e_fhs = mesh.NeighborFh(eh);
	unordered_set<EH> src_ehs;
	for (auto& e_fh : e_fhs) {
		auto& face = mesh.faces(e_fh);
		const auto& f_ehs = face.getEdgeHandle();
		for (auto& f_eh : f_ehs) {
			src_ehs.insert(f_eh);
		}
	}

	const auto& tgt_ehs = mesh.NeighborEh(tgt);
	EH res(-1);
	for (auto& t_eh : tgt_ehs) {
		if (src_ehs.count(t_eh)) continue;
		res = t_eh;
		break;
	}

	return res;

}

void QuadDomain_Builder::generate_quad_mesh(MeshKernel::SurfaceMesh& quad_mesh, vector<int>& vhs) {

	if (quad_edges.empty()) {
		vector<int> tmp;
		get_quad_domains(tmp);
	}

	unordered_map<VH, int> valences;
	for (auto& qe : quad_edges) {
		for (auto& eh : qe) {
			auto& edge = mesh.edges(eh);
			valences[edge.vh1()]++;
			valences[edge.vh2()]++;
		}
	}

	unordered_map<VH, VH> tri2quad;
	for (auto& vp : valences) {
		if (vp.second > 2) {
			vhs.push_back(vp.first);
			auto& v = mesh.vertices(vp.first);
			tri2quad[vp.first] = quad_mesh.AddVertex(v);
		}
	}

	/*for (auto& quad_edge : quad_edges) {
		auto edge = mesh.edges(quad_edge.front());
		VH svh_pre = tri2quad.count(edge.vh1()) ? edge.vh1() : edge.vh2();
		VH vh_pre = svh_pre;
		for (auto& eh : quad_edge) {
			edge = mesh.edges(eh);
			VH vh_cur(edge.vh1() + edge.vh2() - vh_pre);
			if (tri2quad.count(vh_cur)) {
				quad_mesh.AddEdge(tri2quad[svh_pre], tri2quad[vh_cur]);
				svh_pre = vh_cur;
			}
			vh_pre = vh_cur;
		}
	}*/

	for (int i = 0; i < quad_edges.size(); ++i) {
		auto& quad_ehs = quad_edges[i];
		VH vh_pre = quad_edges_srcvh[i];
		vector<VH> quad_vhs = { vh_pre };
		for (auto& eh : quad_ehs) {
			auto& edge = mesh.edges(eh);
			VH vh_cur(edge.vh1() + edge.vh2() - vh_pre);
			if (tri2quad.count(vh_cur)) {
				quad_vhs.push_back(vh_cur);
			}
			vh_pre = vh_cur;
		}

		for (int j = 1; j < quad_vhs.size(); ++j) {
			quad_mesh.AddEdge(tri2quad[quad_vhs[j - 1]], tri2quad[quad_vhs[j]]);
		}

	}

	/*for (auto& vp : quad_mesh.allvertices()) {
		VH vh = vp.first;
		const auto& adjvhs = quad_mesh.NeighborVh(vh);
		bool is_minial_vh = true;
		for (auto& adjvh : adjvhs) {
			if (adjvh < vh) {
				is_minial_vh = false;
				break;
			}
		}
		if (!is_minial_vh) continue;

		for (auto& vh_i : adjvhs) {
			for (auto& vh_j : adjvhs) {
				if (vh_j == vh_i) continue;
				const auto& adjvhs_i = quad_mesh.NeighborVh(vh_i);
				const auto& adjvhs_j = quad_mesh.NeighborVh(vh_j);
				VH vh_common(-1);
				for (auto& adjvh_i : adjvhs_i) {
					if (adjvh_i == vh) continue;
					if (adjvhs_j.count(adjvh_i)) {
						vh_common = adjvh_i;
						break;
					}
				}
				if (vh_common != -1) {
					quad_mesh.AddFace({ vh, vh_i, vh_common, vh_j });
				}
			}
		}
	}*/

	//quad_mesh.eraseComplicatedEdge();


	std::cout << "[Quad Mesh]: edges size = " << quad_mesh.esize() << ", vertices size = " << quad_mesh.vsize() << std::endl;

}