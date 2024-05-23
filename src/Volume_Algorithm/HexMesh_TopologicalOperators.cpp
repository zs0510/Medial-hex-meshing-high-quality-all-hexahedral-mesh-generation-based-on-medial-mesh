#include "HexMesh_TopologicalOperators.h"

void HexMesh_TopologicalOperators::padding(MeshKernel::VolumeMesh& hexmesh, bool add_inner_vertices) {

	int before_cell_size = hexmesh.csize();
	double move_length = std::numeric_limits<double>::max();

	unordered_set<FH> surface_fhs;// padding 前 表面的面
	unordered_set<EH> surface_ehs;// padding 前 表面的边
	unordered_map<VH, VH> vmap;// 旧顶点与新顶点的映射
	unordered_map<FH, FH> fmap;
	unordered_map<EH, FH> emap;
	unordered_map<FH, Vec> faces_normal;// 面的法向量
	unordered_map<VH, Vec> vertices_normal;// 顶点的法向量


	// 生成面的法向量
	for (auto& fp : hexmesh.allfaces()) {
		if (hexmesh.isOnBoundary(fp.first)) {
			const auto& vhs = fp.second.getVertexHandle();
			Vec dir_outer = hexmesh.getFaceCenter(fp.first) - hexmesh.getCellCenter(*(hexmesh.NeighborCh(fp.first).begin()));
			dir_outer.normalize();
			Vec normal(0, 0, 0);
			for (int i = 2; i < vhs.size(); ++i) {
				auto& v0 = hexmesh.vertices(vhs[0]);
				auto& v1 = hexmesh.vertices(vhs[i - 1]);
				auto& v2 = hexmesh.vertices(vhs[i]);
				Vec vec01 = v1 - v0;
				Vec vec02 = v2 - v0;
				Vec N = (vec01.cross(vec02)).normalized();
				if (dir_outer.dot(N) < 0) N *= -1;
				normal += N;
			}
			faces_normal[fp.first] = normal.normalized();
			surface_fhs.insert(fp.first);
		}
	}

	// 偏移取最小边长的 0.05
	double edge_len_max_input = 0.0;
	hexmesh.genAllEdgesLength();
	for (auto& ep : hexmesh.alledges()) {
		if (hexmesh.isOnBoundary(ep.first)) {
			move_length = std::fmin(move_length, ep.second.getLength());
			edge_len_max_input = std::fmax(edge_len_max_input, ep.second.getLength());
			surface_ehs.insert(ep.first);
		}
		
	}
	move_length *= 0.05;
	//std::cout << "[Hex Padding]: Move length = " << move_length << std::endl;

	// 生成顶点法向量
	for (auto& vp : hexmesh.allvertices()) {
		if (hexmesh.isOnBoundary(vp.first)) {
			Vec normal(0, 0, 0);
			const auto& adjfhs = hexmesh.NeighborFh(vp.first);
			for (auto& adjfh : adjfhs) {
				if (faces_normal.count(adjfh)) {
					normal += faces_normal[adjfh];
				}
			}
			vertices_normal[vp.first] = normal.normalized();
			/*double norm2 = vertices_normal[vp.first].norm2();
			if (norm2 > 1.001 || norm2 < 0.999) {
				std::cerr << "Normal norm2 = " << norm2 << std::endl;
			}*/
		}
	}

	// 生成新顶点
	for (auto& vp : hexmesh.allvertices()) {
		if (hexmesh.isOnBoundary(vp.first)) {
			VH vh(vp.first);
			auto& v = hexmesh.vertices(vh);
			Vex vex_new = v;
			if (add_inner_vertices) v = v - vertices_normal[vh] * move_length;// 旧的顶点往里移
			else vex_new = v + vertices_normal[vh] * move_length; // 新的顶点往外移
			vmap[vh] = hexmesh.AddVertex(vex_new);
		}
	}

	// 旧面生成新面
	for (auto& fh : surface_fhs) {
		auto& face = hexmesh.faces(fh);
		auto vhs = face.getVertexHandle();
		for (auto& vh : vhs) {
			vh = vmap[vh];
		}
		fmap[fh] = hexmesh.AddFace(vhs);
	}

	// 旧边生成新面
	for (auto& eh : surface_ehs) {
		auto& edge = hexmesh.edges(eh);
		vector<VH> vhs = {
			edge.vh1(), edge.vh2(), vmap[edge.vh2()], vmap[edge.vh1()]
		};
		emap[eh] = hexmesh.AddFace(vhs);
	}

	// 添加体

	for (auto& fp : fmap) {
		FH fh0 = fp.first, fh1 = fp.second;
		vector<VH> cell_vhs;
		const auto& f_vhs0 = hexmesh.faces(fh0).getVertexHandle();
		const auto& f_vhs1 = hexmesh.faces(fh1).getVertexHandle();
		cell_vhs.insert(cell_vhs.end(), f_vhs0.begin(), f_vhs0.end());
		cell_vhs.insert(cell_vhs.end(), f_vhs1.begin(), f_vhs1.end());
		hexmesh.AddCell(cell_vhs);
	}

	hexmesh.updateAllHandles();

	HexMesh_Smoothing app(hexmesh);
	app.Laplacian_Level(15, false);

	int after_cell_size = hexmesh.csize();

	std::cout << "[Hex Padding]: Before cells count = " << before_cell_size 
		<< ", after cells count = " << after_cell_size << std::endl;

	edge_len_max_input *= 2.0;
	int iter_max = 10;
	while (iter_max--) {
		bool remove_flag = false;
		// 检查有无顶点与邻域边长过大
		for (auto& vp : hexmesh.allvertices()) {
			auto& v = hexmesh.vertices(vp.first);
			const auto& adjvhs = hexmesh.NeighborVh(vp.first);
			int degree = 0;
			Vex v_avg(0, 0, 0);
			for (auto& adjvh : adjvhs) {
				auto& adjv = hexmesh.vertices(adjvh);
				v_avg = v_avg + adjv;
				if ((adjv - v).norm() > edge_len_max_input) {
					degree++;
				}
			}
			if (degree > 1) {
				v = v_avg / adjvhs.size();
				remove_flag = true;
			}

		}
		if (!remove_flag) break;
		
	}

	

}


void HexMesh_TopologicalOperators::dicing(MeshKernel::VolumeMesh& hexmesh, EH input_eh, int num_of_segments) {

	if (!hexmesh.isValid(input_eh) || num_of_segments <= 1) return;

	if (!is_dicing_ok(hexmesh, input_eh)) return;

	std::unordered_map<VH, VH> left2right;// 有方向地记录一条边的两个顶点handle
	std::unordered_set<VH> total_vhs;// 记录所有顶点 handle

	std::unordered_map<FH, bool> visited_fh;
	std::unordered_map<EH, bool> visited_eh;
	std::queue<EH> que;
	que.emplace(input_eh);
	VH vh_left = hexmesh.edges(input_eh).vh1(), vh_right = hexmesh.edges(input_eh).vh2();
	left2right[vh_left] = vh_right;
	visited_eh[input_eh] = true;
	total_vhs.insert(vh_left);
	total_vhs.insert(vh_right);

	while (!que.empty()) {// 使用 BFS, 把所有需要细分的边找到
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
					|| f_edge.vh1() == vh2 || f_edge.vh2() == vh2) continue;// 在同一个面上且不相邻的边

				vh_left = f_edge.vh1(), vh_right = f_edge.vh2();
				bool vh_left_is_ok = false;
				for (auto& adjvh : hexmesh.NeighborVh(vh_left)) {
					if (left2right.count(adjvh)) {
						vh_left_is_ok = true;
						break;
					}
				}
				if (!vh_left_is_ok) std::swap(vh_left, vh_right);// 调整方向, 确保细分时不会出现自交
				que.emplace(f_eh);
				left2right[vh_left] = vh_right;
				visited_eh[f_eh] = true;
				total_vhs.insert(vh_left);
				total_vhs.insert(vh_right);
			}
		}
	}

	std::unordered_map<VH, std::vector<VH>> left2vhs;// 记录细分体需要用到的顶点 handle
	for (auto& vp : left2right) {
		auto& v_left = hexmesh.vertices(vp.first);
		auto& v_right = hexmesh.vertices(vp.second);
		Vec vec = v_right - v_left;

		left2vhs[vp.first].emplace_back(vp.first);
		for (int i = 1; i < num_of_segments; ++i) {
			left2vhs[vp.first].emplace_back(hexmesh.AddVertex(v_left + vec * (i * 1.0 / num_of_segments)));
		}
		left2vhs[vp.first].emplace_back(vp.second);

	}


	std::unordered_map<CH, bool> visited_ch;
	for (auto& ep : visited_eh) {
		auto eh = ep.first;
		const auto& adjchs = hexmesh.NeighborCh(eh);
		for (auto& adjch : adjchs) {
			if (visited_ch.count(adjch)) continue;
			visited_ch[adjch] = true;

			bool is_involved = true;// 检查是否是需要被细分的体
			auto& cell = hexmesh.cells(adjch);
			const auto& c_vhs = cell.getVertexHandle();
			for (auto& c_ch : c_vhs) {
				if (!total_vhs.count(c_ch)) {
					is_involved = false;
					break;
				}
			}
			if (!is_involved) continue;

			const auto& c_fhs = cell.getFaceHandle();// 找到左边那个面
			for (auto& c_fh : c_fhs) {
				auto& face = hexmesh.faces(c_fh);
				const auto& f_vhs = face.getVertexHandle();
				bool is_left_face = true;
				for (auto& f_vh : f_vhs) {
					if (!left2right.count(f_vh)) {
						is_left_face = false;
						break;
					}
				}
				if (is_left_face) {
					for (int i = 1; i <= num_of_segments; ++i) {
						int j = i - 1;
						hexmesh.AddCell({ left2vhs[f_vhs[0]][j], left2vhs[f_vhs[1]][j] ,left2vhs[f_vhs[2]][j] ,left2vhs[f_vhs[3]][j],
							left2vhs[f_vhs[0]][i], left2vhs[f_vhs[1]][i] ,left2vhs[f_vhs[2]][i] ,left2vhs[f_vhs[3]][i] });
					}
					break;// 一个体只有一个左面, 只被细分一次
				}
			}


		}
	}

	for (auto& ep : visited_eh) {// 删去被细分的边
		if (hexmesh.isValid(ep.first)) {
			hexmesh.DeleteEdge(ep.first);
		}
	}

	//hexmesh.updateAllHandles();

}


void HexMesh_TopologicalOperators::get_dicing_vhs(MeshKernel::VolumeMesh& hexmesh, EH input_eh, std::unordered_map<VH, VH>& left2right) {

	std::unordered_set<VH> total_vhs;// 记录所有顶点 handle
	std::unordered_map<FH, bool> visited_fh;
	std::unordered_map<EH, bool> visited_eh;
	std::queue<EH> que;
	que.emplace(input_eh);
	VH vh_left = hexmesh.edges(input_eh).vh1(), vh_right = hexmesh.edges(input_eh).vh2();
	left2right[vh_left] = vh_right;
	visited_eh[input_eh] = true;
	total_vhs.insert(vh_left);
	total_vhs.insert(vh_right);

	while (!que.empty()) {// 使用 BFS, 把所有需要细分的边找到
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
					|| f_edge.vh1() == vh2 || f_edge.vh2() == vh2) continue;// 在同一个面上且不相邻的边

				vh_left = f_edge.vh1(), vh_right = f_edge.vh2();
				bool vh_left_is_ok = false;
				for (auto& adjvh : hexmesh.NeighborVh(vh_left)) {
					if (left2right.count(adjvh)) {
						vh_left_is_ok = true;
						break;
					}
				}
				if (!vh_left_is_ok) std::swap(vh_left, vh_right);// 调整方向, 确保细分时不会出现自交
				que.emplace(f_eh);
				left2right[vh_left] = vh_right;
				visited_eh[f_eh] = true;
				total_vhs.insert(vh_left);
				total_vhs.insert(vh_right);
			}
		}
	}

}

bool HexMesh_TopologicalOperators::is_dicing_ok(MeshKernel::VolumeMesh& hexmesh, EH input_eh) {


	std::unordered_map<VH, VH> left2right;
	get_dicing_vhs(hexmesh, input_eh, left2right);

	std::unordered_set<EH> total_ehs;
	for (auto& vp : left2right) {
		total_ehs.insert(hexmesh.getEdgeHandle(vp.first, vp.second));
	}
	for (auto& cp : hexmesh.allcells()) {
		int degree = 0;
		const auto& ehs = cp.second.getEdgeHandle();
		for (auto& eh : ehs) {
			if (total_ehs.count(eh)) {
				degree++;
			}
		}
		if (degree != 0 && degree != 4) {
			std::cout << "[HexMesh_TopologicalOperators]: Exist singular edges so we skip this dicing operation.\n";
			return false;
		}
	}

	return true;

}