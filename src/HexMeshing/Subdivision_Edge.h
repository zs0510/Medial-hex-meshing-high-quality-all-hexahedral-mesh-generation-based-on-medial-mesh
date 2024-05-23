#pragma once
#include "Kernel/Mesh.h"

namespace Subdivision_Edge {

	void get_total_edges(MeshKernel::VolumeMesh& hexmesh, EH eh, std::vector<int>& res) {
		res.clear();
		if (!hexmesh.isValid(eh)) return;
		std::unordered_map<FH, bool> visited_fh;
		std::unordered_map<EH, bool> visited_eh;
		std::queue<EH> que;
		que.emplace(eh);
		res.emplace_back(eh);
		visited_eh[eh] = true;
		while (!que.empty()) {
			auto cur_eh = que.front();
			que.pop();
			auto& edge = hexmesh.edges(cur_eh);
			VH vh1 = edge.vh1(), vh2 = edge.vh2();
			const auto& adjfhs = hexmesh.NeighborFh(cur_eh);
			for (auto& adjfh : adjfhs) {
				if (visited_fh.count(adjfh)) continue;
				visited_fh[adjfh] = true;
				auto& face = hexmesh.faces(adjfh);
				const auto& f_ehs = face.getEdgeHandle();
				for (auto& f_eh : f_ehs) {
					if (visited_eh.count(f_eh)) continue;
					auto& f_edge = hexmesh.edges(f_eh);
					if (f_edge.vh1() == vh1 || f_edge.vh2() == vh1 
						|| f_edge.vh1() == vh2 || f_edge.vh2() == vh2) continue;
					que.emplace(f_eh);
					res.emplace_back(f_eh);
					visited_eh[f_eh] = true;
				}
			}
		}
		printf("Edges size = %d\n", (int)res.size());
	}

	void subdivision_edges(MeshKernel::VolumeMesh& hexmesh, EH input_eh, int num_of_segments);// 线性细分

	void subdivision_edges_weighted(MeshKernel::VolumeMesh& hexmesh, EH input_eh, double weight_vh1 = 0.5);// 非线性细分, 按给定权值细分

};