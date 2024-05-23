#include "QuadMeshSubdivision_CC.h"

void QuadMeshSubdivision_CC::subdivision(int iter_times) {

	// 1. 初始化
	subdi_mesh = control_mesh;
	for (auto& vp : control_mesh.allvertices()) {
		weights[vp.first][vp.first] = 1.0;
	}

	for (int it = 0; it < iter_times; ++it) {

		unordered_map<FH, VH> face_vhs;
		unordered_map<EH, VH> edge_vhs;

		auto origin_vertices = subdi_mesh.allvertices();
		auto origin_edges = subdi_mesh.alledges();
		auto origin_faces = subdi_mesh.allfaces();
		//unordered_map<VH, unordered_map<VH, double>> new_weights;

		for (auto& fp : subdi_mesh.allfaces()) {
			Vex face_center = subdi_mesh.getFaceCenter(fp.first);
			VH face_vh = face_vhs[fp.first] = subdi_mesh.AddVertex(face_center);
			// 更新权值
			auto& face = subdi_mesh.faces(fp.first);
			const auto& f_vhs = face.getVertexHandle();
			for (auto& f_vh : f_vhs) {
				for (auto& wp : weights[f_vh]) {
					weights[face_vh][wp.first] += wp.second * 0.25;
				}
			}

		}

		for (auto& ep : subdi_mesh.alledges()) {// 请保证输入曲面是封闭的
			Vex edge_midpoint = subdi_mesh.getEdgeMidpoint(ep.first);
			Vex face_point(0, 0, 0);
			const auto& adjfhs = subdi_mesh.NeighborFh(ep.first);
			vector<VH> weighted_vhs = { ep.second.vh1(), ep.second.vh2() };
			for (auto& adjfh : adjfhs) {
				VH f_vh = face_vhs[adjfh];
				face_point += subdi_mesh.vertices(f_vh);
				weighted_vhs.emplace_back(f_vh);
			}
			if (weighted_vhs.size() != 4) {
				std::cerr << "[Error]: Line 46\n";
				return;
			}
			Vex edge_point = edge_midpoint * 0.5 + face_point * 0.25;
			VH edge_vh = edge_vhs[ep.first] = subdi_mesh.AddVertex(edge_point);
			for (auto& w_vh : weighted_vhs) {
				for (auto& wp : weights[w_vh]) {
					weights[edge_vh][wp.first] += wp.second * 0.25;
				}
			}

		}

		for (auto& vp : origin_vertices) {
			VH vh = vp.first;
			int valence = subdi_mesh.valence(vh);
			if (valence == 0) continue;

			double face_weight = 1.0 / valence;
			double edge_weight = 2.0 / valence;
			double vertex_weight = 1.0 - face_weight - edge_weight;

			Vex face_point(0, 0, 0);
			vector<VH> fw_vhs;
			const auto& adjfhs = subdi_mesh.NeighborFh(vh);
			for (auto& adjfh : adjfhs) {
				face_point += subdi_mesh.vertices(face_vhs[adjfh]);
				fw_vhs.emplace_back(face_vhs[adjfh]);
			}
			face_point /= adjfhs.size();

			Vex edge_point(0, 0, 0);
			vector<VH> ew_vhs;
			const auto& adjehs = subdi_mesh.NeighborEh(vh);
			for (auto& adjeh : adjehs) {
				edge_point += subdi_mesh.vertices(edge_vhs[adjeh]);
				ew_vhs.emplace_back(edge_vhs[adjeh]);
			}
			edge_point /= adjehs.size();

			auto& v = subdi_mesh.vertices(vh);
			v = face_point * face_weight + edge_point * edge_weight + v * vertex_weight;
			/*double weight_sum = face_weight + edge_weight + vertex_weight;
			if (weight_sum < 0.9999 || weight_sum > 1.0001) {
				std::cerr << "[Error]: Exist non normalized weight!\n";
				std::cerr << "Vertex weight = " << weight_sum << std::endl;
				return;
			}*/

			//weights[vh].clear();

			unordered_map<VH, double> new_weight;

			//double v_weight_sum = 0;
			for (auto& wp : weights[vh]) {
				new_weight[wp.first] += wp.second * vertex_weight;
				//v_weight_sum += wp.second * vertex_weight;
			}
			

			//double f_weight_sum = 0;
			for (auto& w_vh : fw_vhs) {
				for (auto& wp : weights[w_vh]) {
					new_weight[wp.first] += wp.second * face_weight / valence;
					//f_weight_sum += wp.second * face_weight / valence;
				}
			}

			//double e_weight_sum = 0;
			for (auto& w_vh : ew_vhs) {
				for (auto& wp : weights[w_vh]) {
					new_weight[wp.first] += wp.second * edge_weight / valence;
					//e_weight_sum += wp.second * edge_weight / valence;
				}
			}
			/*std::cout << "Valence = " << valence << ", v_weight = " << vertex_weight << ", v_weight_sum = " << v_weight_sum
				<< ", f_weight = " << face_weight << ", f_weight_sum = " << f_weight_sum
				<< ", e_weight = " << edge_weight << ", e_weight_sum = " << e_weight_sum << std::endl;*/
			weights[vh] = new_weight;

		}

		//// check weights
		//for (auto& vp : origin_vertices) {
		//	double weight_sum = 0;
		//	for (auto& wp : weights[vp.first]) {
		//		weight_sum += wp.second;
		//	}

		//	if (weight_sum < 0.9999 || weight_sum > 1.0001) {
		//		std::cerr << "[Error]: Exist non normalized weight!\n";
		//		std::cerr << "weight = " << weight_sum << std::endl;
		//		return;
		//	}
		//}
		//std::cout << "Vertex weight check success.\n";

		//for (auto& vp : edge_vhs) {
		//	double weight_sum = 0;
		//	for (auto& wp : weights[vp.second]) {
		//		weight_sum += wp.second;
		//	}
		//	
		//	if (weight_sum < 0.9999 || weight_sum > 1.0001) {
		//		std::cerr << "[Error]: Exist non normalized weight!\n";
		//		std::cerr << "weight = " << weight_sum << std::endl;
		//		return;
		//	}
		//}
		//std::cout << "Edge weight check success.\n";

		//for (auto& vp : face_vhs) {
		//	double weight_sum = 0;
		//	for (auto& wp : weights[vp.second]) {
		//		weight_sum += wp.second;
		//	}
		//	if (weight_sum < 0.9999 || weight_sum > 1.0001) {
		//		std::cerr << "[Error]: Exist non normalized weight!\n";
		//		std::cerr << "weight = " << weight_sum << std::endl;
		//		return;
		//	}
		//}
		//std::cout << "Face weight check success.\n";

		for (auto& fp : origin_faces) {

			FH fh(fp.first);
			const auto& f_vhs = fp.second.getVertexHandle();
			const auto& f_ehs = fp.second.getEdgeHandle();

			for (int i = 0; i < 4; ++i) {
				vector<VH> new_f_vhs = f_vhs;
				for (int j = 0; j < 4; ++j) {
					if (j == i) continue;
					if (subdi_mesh.isConnected(f_vhs[i], f_vhs[j])) {
						EH eh = subdi_mesh.getEdgeHandle(f_vhs[i], f_vhs[j]);
						new_f_vhs[j] = edge_vhs[eh];
					} else {
						new_f_vhs[j] = face_vhs[fh];
					}
				}
				subdi_mesh.AddFace(new_f_vhs);
			}

		}
		

		for (auto& ep : origin_edges) {
			subdi_mesh.DeleteEdge(ep.first);
		}


	}

}