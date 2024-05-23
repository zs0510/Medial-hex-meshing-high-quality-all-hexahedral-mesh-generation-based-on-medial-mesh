#include "Subdivision_Butterfly.h"

void Subdivision_Butterfly::subdivision(int iter_times) {

	init_weight_tables();

	for (int it = 0; it < iter_times; ++it) {

		generate_edge_vertex();

		generate_new_faces();

	}


}

void Subdivision_Butterfly::init_weight_tables() {

	//weight_tables[6] = { 0.5, 1.0 / 16.0, -1.0 / 16.0, 0, -1.0 / 16.0, 1.0 / 16.0 };// sum = 0.5

	//weight_tables[3] = { 5.0 / 12.0, -1.0 / 12.0, -1.0 / 12.0 };// paper: { 5.0 / 12.0, -1.0 / 12, -1.0 / 12 }, sum = 0.25

	//weight_tables[4] = { 3.0 / 8.0, 0.0, -1.0 / 8.0, 0.0 };// paper: { 3.0 / 8.0, 0.0, -1.0 / 8.0, 0.0 }, sum = 0.25

	//for (int i = 5; i < 100; ++i) {
	//	if (i == 6) continue;
	//	std::vector<double> weights;
	//	for (int j = 0; j < i; ++j) {
	//		double w = 0.25 + std::cos(2.0 * M_PI * j / i) + 0.5 * std::cos(4.0 * M_PI * j / i);
	//		w /= i;
	//		weights.push_back(w);
	//	}
	//	weight_tables[i] = weights;
	//}

	//for (auto& wp : weight_tables) {
	//	double weight_sum = std::accumulate(wp.second.begin(), wp.second.end(), 0.0);
	//	std::cout << "K = " << wp.second.size() << ", weight_sum = " << weight_sum << std::endl;
	//}

}

void Subdivision_Butterfly::generate_edge_vertex() {

	

	for (auto& ep : mesh.alledges()) {

		EH eh(ep.first);
		auto& e = mesh.edges(ep.first);
		VH vh0 = e.vh1();
		VH vh1 = e.vh2();
		Vex pos;

		if (mesh.isOnBoundary(eh)) {// 在边界

			VH vh2(-1), vh3(-1);
			for (auto& adjvh : mesh.NeighborVh(vh0)) {
				if (adjvh != vh1 && mesh.isOnBoundary(adjvh)) {
					vh2 = adjvh;
					break;
				}
			}

			for (auto& adjvh : mesh.NeighborVh(vh1)) {
				if (adjvh != vh0 && mesh.isOnBoundary(adjvh)) {
					vh3 = adjvh;
					break;
				}
			}

			if (!mesh.isValid(vh2) || !mesh.isValid(vh3)) {
				std::cerr << "[Error]: Incorrect neighbor vhs size!!!\n";
				continue;
			}

			pos = mesh.vertices(vh0) * nine_sixteen + mesh.vertices(vh1) * nine_sixteen
				+ mesh.vertices(vh2) * -1.0 * one_sixteen + mesh.vertices(vh3) * -1.0 * one_sixteen;

		} else {// 不在边界

			int valence0 = mesh.valence(vh0), valence1 = mesh.valence(vh1);
			if (valence0 != 6 && valence1 == 6) {// 若有度为6的顶点, 则优先将其置于 vh0
				std::swap(vh0, vh1);
				std::swap(valence0, valence1);
			}

			if (valence0 != 6) {// 两个都不是规则点

				Vex pos_0, pos_1;
				get_vex_from_extraordinary_vertices(vh0, vh1, pos_0);
				get_vex_from_extraordinary_vertices(vh1, vh0, pos_1);
				pos = (pos_0 + pos_1) * 0.5;

			} else if (valence1 == 6) {// 两个都是规则点

				auto adjvhs = mesh.NeighborVh(vh0);
				auto adjvhs1 = mesh.NeighborVh(vh1);
				adjvhs.insert(adjvhs1.begin(), adjvhs1.end());
				adjvhs.erase(vh0);
				adjvhs.erase(vh1);
				std::vector<VH> vh_b;
				for (auto& vh : adjvhs) {
					if (mesh.isConnected(vh0, vh) && mesh.isConnected(vh1, vh)) {
						vh_b.push_back(vh);
						if (vh_b.size() == 2) break;
					}
				}
				assert(vh_b.size() == 2);
				adjvhs.erase(vh_b[0]);
				adjvhs.erase(vh_b[1]);
				std::vector<VH> vh_c;
				for (auto& vh : adjvhs) {
					if (mesh.isConnected(vh_b[0], vh) || mesh.isConnected(vh_b[1], vh)) {
						vh_c.push_back(vh);
						if (vh_c.size() == 4) break;
					}
				}
				assert(vh_c.size() == 4);

				pos = mesh.vertices(vh0) * 0.5 + mesh.vertices(vh1) * 0.5 + mesh.vertices(vh_b[0]) * one_eight + mesh.vertices(vh_b[1]) * one_eight
					- mesh.vertices(vh_c[0]) * one_sixteen - mesh.vertices(vh_c[1]) * one_sixteen - mesh.vertices(vh_c[2]) * one_sixteen - mesh.vertices(vh_c[3]) * one_sixteen;
				

			} else {// 一个规则点, 一个不规则点
				
				get_vex_from_extraordinary_vertices(vh1, vh0, pos);

			}

		}// end 不在边界

		edge_vh[eh] = mesh.AddVertex(pos);

	}


}

void Subdivision_Butterfly::generate_new_faces() {

	auto oldfaces = mesh.allfaces();
	for (auto& fp : oldfaces) {

		if (!mesh.isValid(fp.first)) continue;

		const auto& old_ehs = fp.second.getEdgeHandle();
		const auto& old_vhs = fp.second.getVertexHandle();

		if (old_ehs.size() != 3) continue;

		mesh.AddFace({ edge_vh[old_ehs[0]], edge_vh[old_ehs[1]], edge_vh[old_ehs[2]] });

		for (int i = 0; i < 3; ++i) {
			auto cur = old_vhs[i];
			auto pre = old_vhs[(i + 2) % 3];
			auto next = old_vhs[(i + 1) % 3];
			VH new_pre(-1), new_next(-1);
			for (auto eh : old_ehs) {
				auto& e = mesh.edges(eh);
				if (e.vh1() != cur && e.vh2() != cur) continue;
				int vex = e.vh1() + e.vh2() - cur;
				if (new_pre == -1 && vex == pre) {
					new_pre = edge_vh[eh];
				}
				if (new_next == -1 && vex == next) {
					new_next = edge_vh[eh];
				}
			}
			assert(new_pre != -1 && new_next != -1);
			mesh.AddFace({ new_pre, cur, new_next });
		}

	}

}

void Subdivision_Butterfly::get_ordered_vvhs(VH center, VH start, std::vector<VH>& res) {

	auto adjvhs = mesh.NeighborVh(center);
	std::unordered_set<EH> oppo_ehs;
	std::unordered_set<VH> visited;
	for (auto& fh : mesh.NeighborFh(center)) {
		auto& f = mesh.faces(fh);
		const auto& ehs = f.getEdgeHandle();
		for (auto& eh : ehs) {
			auto& e = mesh.edges(eh);
			if (e.vh1() != center && e.vh2() != center) {
				oppo_ehs.insert(eh);
				break;
			}
		}
	}
	VH pre_vh = start, cur_vh(-1);
	res.push_back(start);
	visited.insert(start);
	//adjvhs.erase(start);
	for (auto& vh : adjvhs) {
		if (vh != start && mesh.isConnected(vh, start)) {
			cur_vh = vh;
			break;
		}
	}
	while (cur_vh != start && cur_vh != -1) {
		res.push_back(cur_vh);
		visited.insert(cur_vh);
		//adjvhs.erase(cur_vh);
		VH next_vh(-1);
		for (auto& vh : adjvhs) {
			if (visited.count(vh)) continue;
			if (mesh.isConnected(vh, cur_vh)) {
				if (oppo_ehs.count(mesh.getEdgeHandle(vh, cur_vh))) {
					next_vh = vh;
					break;
				}
			}
		}
		//assert(next_vh != -1);// center 在边界才可能会出现 next_vh = -1
		//std::cout << "pre_vh = " << pre_vh << ", cur_vh = " << cur_vh << ", next_vh = " << next_vh << std::endl;
		pre_vh = cur_vh;
		cur_vh = next_vh;
	}
	if (res.size() != adjvhs.size()) {
		std::cerr << "Center: " << center << 
			", res.size() = " << res.size() << 
			", adjvhs.size() = " << adjvhs.size() << std::endl;
		for (auto& vh : res) {
			std::cout << vh << " -> ";
		}
		std::cout << std::endl;
	}

	assert(res.size() == adjvhs.size());

}

void Subdivision_Butterfly::get_vex_from_extraordinary_vertices(VH center, VH start, Vex& res) {

	int valence = mesh.valence(center);

	if (valence == 3) {

		auto adjvhs = mesh.NeighborVh(center);
		adjvhs.erase(start);
		assert(adjvhs.size() == 2);

		auto it = adjvhs.begin();
		VH vh_left = *it++;
		VH vh_right = *it;

		res = mesh.vertices(start) * five_twelve + mesh.vertices(center) * 0.75 - mesh.vertices(vh_left) * one_twelve
			- mesh.vertices(vh_right) * one_twelve;


	} else if (valence == 4) {

		auto adjvhs = mesh.NeighborVh(center);
		adjvhs.erase(start);
		assert(adjvhs.size() == 3);

		VH vh_rear(-1);
		for (auto& vh : adjvhs) {
			if (mesh.isConnected(start, vh)) continue;
			vh_rear = vh;
			break;
		}

		res = mesh.vertices(start) * three_eight + mesh.vertices(center) * 0.75 - mesh.vertices(vh_rear) * one_eight;
		

	} else {// valence >= 5

		Vex pos(0, 0, 0);
		double weight_center = 1.0;
		std::vector<VH> vvhs;
		get_ordered_vvhs(center, start, vvhs);
		assert(vvhs.size() == valence);
		for (int i = 0; i < valence; ++i) {
			double weight_s = (0.25 + std::cos(2.0 * M_PI * i / valence) + 0.5 * std::cos(4.0 * M_PI * i / valence)) / valence;
			weight_center -= weight_s;
			pos += (mesh.vertices(vvhs[i]) * weight_s);
		}
		pos += (mesh.vertices(center) * weight_center);
		res = pos;

	}

}