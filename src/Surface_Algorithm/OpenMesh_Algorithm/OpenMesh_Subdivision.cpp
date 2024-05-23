#include "OpenMesh_Subdivision.h"

void iGame_OpenMesh_Subdivision::subdivision_butterfly(int iter_times) {

	mesh.add_property(edge_vh);
	mesh.request_vertex_status();
	mesh.request_face_status();

	init_weight_tables();

	for (int it = 0; it < iter_times; ++it) {

		generate_edge_vertex();

		generate_new_faces();

		//mesh.garbage_collection();

	}

	mesh.garbage_collection();

	mesh.remove_property(edge_vh);
	mesh.release_vertex_status();
	mesh.release_face_status();

}

void iGame_OpenMesh_Subdivision::generate_edge_vertex() {

	const double nine_sixteen = 9.0 / 16.0;
	const double negative_one_sixteen = -1.0 / 16.0;

	for (auto eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit) {// 为每条边增加一个边点

		OM_EH eh = *eit;
		OM_VH v0 = eit->v0();
		OM_VH v1 = eit->v1();
		OM_HEH heh = eit->h0();
		OM_HEH op_heh = mesh.opposite_halfedge_handle(heh);

		if (mesh.is_boundary(eh)) {// s[-1] = -1/16, s[0] = s[1] = 9/16, s[2] = -1/16

			OM_VH v2, v3;// v0/v1的另外一个边界点
			for (auto vvit = mesh.vv_begin(v0); vvit.is_valid(); ++vvit) {
				if (mesh.is_boundary(*vvit) && (*vvit) != v1) {
					v2 = *vvit;
					break;
				}
			}
			for (auto vvit = mesh.vv_begin(v1); vvit.is_valid(); ++vvit) {
				if (mesh.is_boundary(*vvit) && (*vvit) != v0) {
					v3 = *vvit;
					break;
				}
			}
			if (!mesh.is_valid_handle(v2) || !mesh.is_valid_handle(v3)) {
				std::cerr << "[Error]: Find invalid handle!!!\n";
				continue;
			}
			
			OM_Pos pos = mesh.point(v0)* nine_sixteen + mesh.point(v1) * nine_sixteen + mesh.point(v2) * negative_one_sixteen + mesh.point(v3) * negative_one_sixteen;
			mesh.property(edge_vh, heh) = mesh.property(edge_vh, op_heh) = mesh.add_vertex(pos);

		} else {

			OM_Pos pos(0.0, 0.0, 0.0);
			std::vector<OM_VH> evhs = { v0, v1 };
			for (int i = 0; i < 2; ++i) {
				OM_VH cur_vh = evhs[i], adj_vh = evhs[(i + 1) % 2];
				int valence = mesh.valence(cur_vh);
				const auto& weights = weight_tables[valence];
				std::vector<OM_VH> neighbor_vhs;
				int beg_idx = 0;
				for (auto vvit = mesh.vv_begin(cur_vh); vvit.is_valid(); ++vvit) {
					if (*vvit == adj_vh) beg_idx = neighbor_vhs.size();
					neighbor_vhs.push_back(*vvit);
				}
				if (neighbor_vhs.size() != valence) {
					std::cerr << "[Error]: Incorrect neighbor vhs size!!!\n";
					continue;
				}
				for (int j = 0; j < valence; ++j) {
					int vj = (beg_idx + j) % valence;
					pos += (mesh.point(neighbor_vhs[vj]) * weights[j]);
				}

				//mesh.property(edge_vh, heh) = mesh.property(edge_vh, op_heh) = mesh.add_vertex(pos);
				mesh.property(edge_vh, heh) = mesh.property(edge_vh, op_heh) = mesh.add_vertex(mesh.calc_edge_midpoint(eh));
			}

		}

	}

}

void iGame_OpenMesh_Subdivision::generate_new_faces() {

	auto old_faces = mesh.faces();

	std::vector<std::vector<OM_VH>> total_new_faces;
	total_new_faces.reserve(mesh.n_faces() * 3);

	for (auto& fh : old_faces) {

		if (!mesh.is_valid_handle(fh)) continue;

		std::vector<std::vector<OM_VH>> new_faces;
		std::vector<OM_VH> edge_face;

		for (auto f_he_it = mesh.fh_begin(fh); f_he_it.is_valid(); ++f_he_it) {

			OM_HEH pre_he = mesh.prev_halfedge_handle(*f_he_it);// 上一条半边
			OM_HEH cur_he = *f_he_it;
			
			OM_VH pre_vh = mesh.property(edge_vh, pre_he);
			OM_VH cur_vh = mesh.from_vertex_handle(cur_he);
			OM_VH next_vh = mesh.property(edge_vh, cur_he);
			if (pre_vh != cur_vh && cur_vh != next_vh && next_vh != pre_vh) {
				new_faces.push_back({ pre_vh, cur_vh, next_vh });
			} else {
				std::cerr << "[Error]: Complex vertex handle!!! ";
				std::cout << "VH: " << pre_vh << " & " << cur_vh << " & " << next_vh << ".\n";
			}
			edge_face.push_back(cur_vh);
			
		}

		mesh.delete_face(fh);

		//new_faces.push_back(edge_face);
		//for (auto& new_face : new_faces) mesh.add_face(new_face);

		total_new_faces.push_back(edge_face);
		total_new_faces.insert(total_new_faces.end(), new_faces.begin(), new_faces.end());

	}

	for (auto& new_face : total_new_faces) mesh.add_face(new_face);

}

void iGame_OpenMesh_Subdivision::init_weight_tables() {

	weight_tables[6] = { 0.5, 1.0 / 16.0, -1.0 / 16.0, 0, -1.0 / 16.0, 1.0 / 16.0 };

	weight_tables[3] = { 5.0 / 12.0, -1.0 / 12, -1.0 / 12 };

	weight_tables[4] = { 3.0 / 8.0, 0.0, -1.0 / 8.0, 0.0 };

	for (int i = 5; i < 20; ++i) {
		if (i == 6) continue;
		std::vector<double> weights;
		for (int j = 0; j < i; ++j) {
			double w = 0.25 + std::cos(2.0 * M_PI * j / i) + 0.5 * std::cos(4.0 * M_PI * j / i);
			w /= i;
			weights.push_back(w);
		}
		weight_tables[i] = weights;
	}

	for (auto& wp : weight_tables) {
		double weight_sum = std::accumulate(wp.second.begin(), wp.second.end(), 0.0);
		std::cout << "K = " << wp.second.size() << ", weight_sum = " << weight_sum << std::endl;
	}

}

