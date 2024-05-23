#include "Add_Cell.h"

void HexAddCell::add_cell_2fh(MeshKernel::VolumeMesh& mesh, FH fh1, FH fh2) {
	if (!mesh.isValid(fh1) || !mesh.isValid(fh2)) {
		return;
	}

	auto chs1 = mesh.NeighborCh(fh1);
	auto chs2 = mesh.NeighborCh(fh2);
	for (auto& ch : chs1) {
		if (std::find(chs2.begin(), chs2.end(), ch) != chs2.end()) {
			std::cout << "Error: exist a cell contain input fhs." << std::endl;
			return;
		}
	}
	auto& face1 = mesh.faces(fh1);
	auto& face2 = mesh.faces(fh2);

	if (mesh.isConnected(fh1, fh2)) {// 有一条公共边（有公共顶点不算邻接）

		auto ehs1 = face1.getEdgeHandle();
		auto ehs2 = face2.getEdgeHandle();
		EH common_eh(-1);
		for (auto& eh : ehs1) {
			if (std::find(ehs2.begin(), ehs2.end(), eh) != ehs2.end()) {
				common_eh = eh;
				break;
			}
		}
		if (common_eh == -1) return;
		std::vector<std::vector<VH>> vhs{ face1.getVertexHandle() , face2.getVertexHandle() };
		auto& common_edge = mesh.edges(common_eh);
		std::vector<VH> common_vhs{ common_edge.vh1() , common_edge.vh2() };
		std::vector<std::vector<VH>> adjvhs(2);
		for (int i = 0; i < 2; ++i) {
			for (auto& vh : vhs[i]) {
				if (vh == common_vhs[0] || vh == common_vhs[1]) continue;
				if (mesh.isConnected(vh, common_vhs[0])) adjvhs[0].push_back(vh);
				if (mesh.isConnected(vh, common_vhs[1])) adjvhs[1].push_back(vh);
			}
		}
		assert(adjvhs[0].size() == 2 && adjvhs[1].size() == 2 && " adjvhs size is wrong! ");
		std::vector<VH> new_vhs(2, VH(-1));
		for (int i = 0; i < 2; ++i) {
			auto& v0 = mesh.vertices(common_vhs[i]);
			auto& v1 = mesh.vertices(adjvhs[i][0]);
			auto& v2 = mesh.vertices(adjvhs[i][1]);
			auto vec = (v1 - v0) + (v2 - v0);
			auto new_v = v0 + vec;
			new_vhs[i] = mesh.AddVertex(new_v);
		}
		std::vector<VH> new_cell = { common_vhs[0], adjvhs[0][0], new_vhs[0], adjvhs[0][1], common_vhs[1], adjvhs[1][0], new_vhs[1], adjvhs[1][1] };
		auto new_ch = mesh.AddCell(new_cell);
		std::cout << "add new cell success. new ch = " << new_ch << std::endl;
	} else {// 无公共边

		std::vector<std::vector<VH>> vhs{ face1.getVertexHandle(), face2.getVertexHandle() };

		//{// 判断是否有顶点相连，若有, 则根据之前相连的拓扑规则推导出新体的拓扑
		//	//int connect_count = 0;
		//	std::vector<VH> connect_vhs1;
		//	std::vector<VH> connect_vhs2;
		//	for (auto& vh1 : vhs[0]) {
		//		for (auto& vh2 : vhs[1]) {
		//			if (mesh.isConnected(vh1, vh2)) {
		//				connect_vhs1.emplace_back(vh1);
		//				connect_vhs2.emplace_back(vh2);
		//				break;
		//			}
		//		}
		//	}
		//	if (connect_vhs1.size() == 2) {
		//		
		//		std::cout << "[Add Cell]: We are using new method!!!\n";
		//		return;
		//	}/* else {
		//		for (auto& vh : connect_vhs) {
		//			std::cout << vh << std::endl;
		//		}
		//		std::cout << "[Add Cell]: We not use new method, beacuse connect vhs size = " << connect_vhs.size() << std::endl;
		//	}*/


		//}

		std::vector<Vex> centroid{ mesh.getFaceCenter(fh1), mesh.getFaceCenter(fh2) };
		Vec vec_c = centroid[1] - centroid[0];
		vec_c.normalize();
		typedef std::pair<std::pair<int, int>, double> Node;
		std::vector<Node> nodes;
		for (auto& vh1 : vhs[0]) {
			auto& v1 = mesh.vertices(vh1);
			for (auto& vh2 : vhs[1]) {
				auto& v2 = mesh.vertices(vh2);
				Vec vec12 = (v2 - v1).normalize();
				double cos_theta = vec_c * vec12;
				nodes.emplace_back(std::pair<int, int>(vh1, vh2), cos_theta);
			}
		}
		std::sort(nodes.begin(), nodes.end(), [](Node& n1, Node& n2) {
			return n1.second > n2.second;
			});
		std::unordered_map<int, int> determinate;

		for (int i = 0; i < nodes.size(); ++i) {
			auto& vh1 = nodes[i].first.first;
			auto& vh2 = nodes[i].first.second;
			auto& cos_theta = nodes[i].second;
			if (determinate.count(vh1) || determinate.count(vh2) || cos_theta < 0) continue;
			determinate[vh1] = vh2;
			determinate[vh2] = vh1;
		}
		assert(determinate.size() == 8 && "Error: vh pairs are invalid.");
		std::vector<VH> new_vhs = vhs[0];
		for (auto& vh : vhs[0]) {
			new_vhs.push_back((VH)determinate[vh]);
		}
		auto new_ch = mesh.AddCell(new_vhs);
		std::cout << "add new cell success. new ch = " << new_ch << std::endl;
	}
}

void HexAddCell::add_cell_3fh(MeshKernel::VolumeMesh& mesh, FH fh1, FH fh2, FH fh3) {

	unordered_map<int, int> vh_count;

	auto& face1 = mesh.faces(fh1), & face2 = mesh.faces(fh2), & face3 = mesh.faces(fh3);
	vector<vector<VH>> vhs = { face1.getVertexHandle() , face2.getVertexHandle(), face3.getVertexHandle() };
	vector<vector<EH>> ehs = { face1.getEdgeHandle() , face2.getEdgeHandle(), face3.getEdgeHandle() };
	vector<EH> common_ehs;
	VH common_vh(-1);
	for (auto& _vhs : vhs) {
		for (auto& _vh : _vhs) {
			vh_count[_vh]++;
			if (vh_count[_vh] == 3) {
				if (common_vh != -1) {
					std::cout << "Error: Support special case only." << std::endl;
					return;
				}
				common_vh = _vh;
			}
		}
	}
	if (common_vh == -1) {
		std::cout << "Error: Support special case only." << std::endl;
		return;
	}
	for (auto& _ehs : ehs) {
		for (auto& _eh : _ehs) {
			auto& _edge = mesh.edges(_eh);
			if (_edge.vh1() == common_vh || _edge.vh2() == common_vh) {
				if (find(common_ehs.begin(), common_ehs.end(), _eh) == common_ehs.end())
					common_ehs.push_back(_eh);
			}
		}
	}
	if (common_ehs.size() != 3) {
		std::cout << "Error: Support special case only." << std::endl;
		return;
	}
	Vec _vec;
	for (auto& eh : common_ehs) {
		auto& edge = mesh.edges(eh);
		auto vh1 = edge.vh1(), vh2 = edge.vh2();
		if (vh1 == common_vh) {
			_vec += (mesh.vertices(vh2) - mesh.vertices(vh1));
		} else {
			_vec -= (mesh.vertices(vh2) - mesh.vertices(vh1));
		}
	}
	auto new_v = mesh.vertices(common_vh) + _vec;
	auto new_vh = mesh.AddVertex(new_v);
	auto new_vhs = vhs[0];
	unordered_map<int, int> bi_map;
	for (int i = 1; i < 3; ++i) {
		auto& _ehs = ehs[i];
		for (auto& _eh : _ehs) {
			auto& _edge = mesh.edges(_eh);
			VH vh1 = _edge.vh1(), vh2 = _edge.vh2();
			bool flag1 = find(new_vhs.begin(), new_vhs.end(), vh1) != new_vhs.end();
			bool flag2 = find(new_vhs.begin(), new_vhs.end(), vh2) != new_vhs.end();
			if (flag1 && !flag2) {
				bi_map[vh1] = vh2;
			} else if (!flag1 && flag2) {
				bi_map[vh2] = vh1;
			}
		}
	}
	if (bi_map.size() != 3) {
		std::cout << "Error: Support special case only." << std::endl;
		return;
	}
	for (int i = 0; i < 4; ++i) {
		auto vh = new_vhs[i];
		if (bi_map.count(vh)) new_vhs.push_back((VH)bi_map[vh]);
		else new_vhs.push_back(new_vh);
	}
	auto new_ch = mesh.AddCell(new_vhs);
	std::cout << "add new cell success. new ch = " << new_ch << std::endl;

}

void HexAddCell::add_cell_eh_fh(MeshKernel::VolumeMesh& mesh, EH eh, FH fh) {
	if (!mesh.isValid(eh) || !mesh.isValid(fh)) {
		return;
	}
	auto& face = mesh.faces(fh);
	auto& edge = mesh.edges(eh);
	auto vhs = face.getVertexHandle();
	auto ehs = face.getEdgeHandle();
	// check topology
	VH on_face_vh(-1), out_face_vh(-1);
	bool find_flag1 = std::find(vhs.begin(), vhs.end(), edge.vh1()) != vhs.end();
	bool find_flag2 = std::find(vhs.begin(), vhs.end(), edge.vh2()) != vhs.end();
	if (find_flag1 && !find_flag2) {
		on_face_vh = edge.vh1();
		out_face_vh = edge.vh2();
	} else if (!find_flag1 && find_flag2) {
		on_face_vh = edge.vh2();
		out_face_vh = edge.vh1();
	} else {
		std::cout << "Please check input parameter." << std::endl;
		return;
	}
	Vec vec_e = mesh.vertices(out_face_vh) - mesh.vertices(on_face_vh);
	std::unordered_map<VH, VH> mp;
	mp[on_face_vh] = out_face_vh;
	for (auto& vh : vhs) {
		if (vh == on_face_vh) continue;
		auto& v_from = mesh.vertices(vh);
		Vex v_to = v_from + vec_e;
		mp[vh] = mesh.AddVertex(v_to);
	}
	assert(mp.size() == 4 && "Error: vh pairs are invalid.");
	std::vector<VH> new_vhs = vhs;
	for (auto& vh : vhs) {
		new_vhs.push_back((VH)mp[vh]);
	}
	auto new_ch = mesh.AddCell(new_vhs);
	std::cout << "add new cell success. new ch = " << new_ch << std::endl;

}

void HexAddCell::add_cell_fh(MeshKernel::VolumeMesh& mesh, FH fh) {

	if (!mesh.isValid(fh)) return;

	if (mesh.NeighborCh(fh).size() != 1) return;

	CH ch = *(mesh.NeighborCh(fh).begin());
	auto& cell = mesh.cells(ch);
	auto& face = mesh.faces(fh);
	auto vhs_f = face.getVertexHandle();
	reverse(vhs_f.begin(), vhs_f.end());
	double avg_len = 0.f;
	for (auto& eh : face.getEdgeHandle()) {
		avg_len += mesh.getLength(eh);
	}
	avg_len /= 4;
	auto centroid_cell = mesh.getCellCenter(ch);
	auto centroid_face = mesh.getFaceCenter(fh);
	auto dir = (centroid_face - centroid_cell).normalize();
	auto move = dir * avg_len;
	vector<VH> new_vhs = vhs_f;
	for (auto& _vh : vhs_f) {
		auto pos = mesh.vertices(_vh);
		pos += move;
		auto new_vh = mesh.AddVertex(pos);
		new_vhs.push_back(new_vh);
	}

	mesh.AddCell(new_vhs);

}