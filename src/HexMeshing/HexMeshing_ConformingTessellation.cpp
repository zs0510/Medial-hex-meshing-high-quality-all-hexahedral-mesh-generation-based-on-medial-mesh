#include "HexMeshing_ConformingTessellation.h"

void HexMeshing_ConformingTessellation::conforming_split(VH algorithm_insert) {

	rebulid_bvh();
	std::cout << "[Conforming Tessellation]: Rebulid bvh success.\n";

	init_nodes_map(algorithm_insert);
	std::cout << "[Conforming Tessellation]: Init nodes map success.\n";

	// 选择需要判定的边, 以及每条边的最小段数
	calc_segnum();
	std::cout << "[Conforming Tessellation]: Calc the num of segments success.\n";

	calc_constraint();
	std::cout << "[Conforming Tessellation]: Calc constraints success.\n";

	opt_conforming_split();
	std::cout << "[Conforming Tessellation]: Optimize comforming split success.\n";

	hexmesh_origin = hexmesh;
	for (auto& ep : hexeh_to_minsegnum) {
		if (!hexmesh.isValid(ep.first) || ep.second < 2) continue;
		std::cout << "\tSubdivision Edge: EH = " << ep.first << ", the number of segments is " << ep.second << "\n";
		Subdivision_Edge::subdivision_edges(hexmesh, ep.first, ep.second);
	}

	// 为每段 分支-分支 曲线匹配源面与目标面
	assign_sweep_face();
	std::cout << "[Conforming Tessellation]: Assign sweep face success.\n";

	// 扫掠, 生成六面体网格
	sweep_curve_skel();
	std::cout << "[Conforming Tessellation]: Sweep curve skeleton success.\n";

	//// 避免特别尖锐的锥体
	//cone_detection();

}

void HexMeshing_ConformingTessellation::rebulid_bvh() {

	delete bvh_tree;
	bvh_tree = nullptr;
	bvh_tree = new BVH_Tree();
	bvh_tree->buildBVH_Tree(hexmesh);

}

void HexMeshing_ConformingTessellation::init_nodes_map(VH algorithm_insert) {

	// 标记所有纯曲线边
	for (auto& ep : sklmesh.alledges()) {
		if (sklmesh.NeighborFh(ep.first).empty()) {
			curve_skel_ehs.insert(ep.first);
		}
	}

	// 识别所有分支节点
	for (auto& vp : sklmesh.allvertices()) {

		auto& ehs = sklmesh.NeighborEh(vp.first);
		//int hanging_count = 0;// 记录悬挂边的数目
		vector<EH> hanging_ehs;
		for (auto& eh : ehs) {
			if (curve_skel_ehs.count(eh)) {
				hanging_ehs.emplace_back(eh);
			}
		}
		if (hanging_ehs.empty()) continue;// 无悬挂边, 不作为分支结点考虑
		int faces_count = sklmesh.NeighborFh(vp.first).size();

		if (faces_count != 0 || hanging_ehs.size() > 2 || is_corner(vp.first) || vp.first == algorithm_insert) {// 线面分支 或者 线线分支 或者 线线角点
			auto& v = sklmesh.vertices(vp.first);
			Vector3d src(v.x(), v.y(), v.z());
			for (auto& adjeh : hanging_ehs) {// 只需选择悬挂边作为射线
				auto& e = sklmesh.edges(adjeh);
				VH adjvh(e.vh1() + e.vh2() - vp.first);
				auto& adjv = sklmesh.vertices(adjvh);
				Vec dir = (adjv - v).normalized();// 射线由分支节点射向外面

				Ray ray(src, Vector3d(dir.x(), dir.y(), dir.z()));
				int rayid = rays.size();
				sklvh_to_rayid[vp.first].emplace_back(rayid);// 一个分支顶点可能会发射出多条射线
				sklvh_skleh_to_rayid[vp.first][adjeh] = rayid;
				rayid_to_sklvh_skleh.push_back(pair<VH, EH>(vp.first, adjeh));
				rays.emplace_back(ray);
				recast_ray(rayid, vp.first, adjeh);

			}

			branch_sklvh.insert(vp.first);

		} else if (hanging_ehs.size() == 2) {
			joint_sklvh.insert(vp.first);

		} else if (hanging_ehs.size() == 1) {
			end_sklvh.insert(vp.first);

		}

	}

	//rayid_to_quadfh.resize(rays.size());
	//rayid_to_intersection.resize(rays.size());
	std::cout << "[Conforming Tessellation]: branch nodes size = " << branch_sklvh.size() << ", ending nodes size = " << end_sklvh.size() << ".\n";
}

void HexMeshing_ConformingTessellation::calc_segnum() {

	// 提取所有需要判定的边
	int multi_intersect_face_count = 0;
	for (auto& fp : hexfh_to_rayid) {
		if (fp.second.empty()) continue;

		auto face = hexmesh.faces(fp.first);
		auto fvhs = face.getVertexHandle();
		EH eh1 = hexmesh.getEdgeHandle(fvhs[0], fvhs[1]);// 两条边拓扑应正交
		if (hexeh_to_minsegnum.count(hexmesh.getEdgeHandle(fvhs[2], fvhs[3]))) {// 尽量使用已确定为决策边的边(减小决策变量个数, 加速整数规划求解)
			eh1 = hexmesh.getEdgeHandle(fvhs[2], fvhs[3]);
		} else {
			hexeh_to_minsegnum[eh1] = 1;
		}
		EH eh2 = hexmesh.getEdgeHandle(fvhs[0], fvhs[3]);
		if (hexeh_to_minsegnum.count(hexmesh.getEdgeHandle(fvhs[1], fvhs[2]))) {
			eh2 = hexmesh.getEdgeHandle(fvhs[1], fvhs[2]);
		} else {
			hexeh_to_minsegnum[eh2] = 1;
		}
		hexfh_to_hexeh[fp.first] = pair<EH, EH>(eh1, eh2);

	}
	std::cout << "[Conforming Tessellation]: We find " << multi_intersect_face_count << " faces have multi intersections.\n";

	ray_constraints.reserve(rays.size());
	for (int rid = 0; rid < rays.size(); ++rid) {
		auto hexfh = rayid_to_hexfh[rid];
		auto hexehs = hexmesh.faces(hexfh).getEdgeHandle();
		for (auto& hexeh : hexehs) {
			if (hexeh == hexfh_to_hexeh[hexfh].first || hexeh == hexfh_to_hexeh[hexfh].second) {
				RayConstraint rc;
				rc.rid = rid;
				rc.hexfh = hexfh;
				rc.hexeh = hexeh;
				ray_constraints.emplace_back(rc);
			}
			
		}
	}
	std::cout << "[Conforming Tessellation]: We initialize " << ray_constraints.size() << " ray constraints.\n";

	// 确定面边段数 与 光线段数 的约束
	for (auto& hr : hexfh_to_rayid) {
		if (hr.second.empty()) continue;

		auto face = hexmesh.faces(hr.first);
		auto fvhs = face.getVertexHandle();
		EH eh1 = hexfh_to_hexeh[hr.first].first, eh2 = hexfh_to_hexeh[hr.first].second;
		// 计算边的方向
		Vec vec1 = (hexmesh.vertices(fvhs[1]) - hexmesh.vertices(fvhs[0])).normalized();
		Vec vec2 = (hexmesh.vertices(fvhs[3]) - hexmesh.vertices(fvhs[0])).normalized();
		// 计算交点间的线段在边上的最大投影长度
		double len1 = 0, len2 = 0;
		auto& points = hr.second;
		for (int i = 0; i < points.size(); ++i) {
			for (int j = i + 1; j < points.size(); ++j) {
				Vec diff = rayid_to_intersection[points[j]] - rayid_to_intersection[points[i]];
				len1 = max(len1, abs(diff.dot(vec1)));
				len2 = max(len2, abs(diff.dot(vec2)));
			}
		}
		if (points.size() > 1) {
			multi_intersect_face_count++;
		}
		// 取投影长度最大的边作为细分边
		EH eh_long = eh1, eh_short = eh2;
		if (len1 < len2) {
			swap(eh_long, eh_short);
		}

		hexeh_to_minsegnum[eh_long] = max(hexeh_to_minsegnum[eh_long], int(hr.second.size()));
		//vector<pair<VH, EH>> sklvh_skleh;
		//for (auto& rid : hr.second) {
		//	sklvh_skleh.emplace_back(rayid_to_sklvh_skleh[rid]);// 一条边是 求和
		//	hexeh_to_curve[eh_short].push_back({ rayid_to_sklvh_skleh[rid] });// 一条边是 全等
		//}
		//hexeh_to_curve[eh_long].emplace_back(sklvh_skleh);

		vector<int> rcids;
		for (auto& rid : hr.second) {
			int rcid1 = rid * 2, rcid2 = rid * 2 + 1;
			if (ray_constraints[rcid1].hexeh != eh_long) {
				swap(rcid1, rcid2);
			}
			if (ray_constraints[rcid1].hexeh != eh_long || ray_constraints[rcid2].hexeh != eh_short) {
				std::cerr << "\n[Conforming Tessellation]: Error: we get invalid ray constraints pair.\n";
				continue;
			}
			rcids.emplace_back(rcid1);
			hexeh_to_rcid[eh_short].push_back({ rcid2 });
			//hexeh_rayid_to_rcid[eh_short][rid] = rcid2;
			hexeh_rayid_to_rcid[eh_long][rid] = rcid1;
		}
		hexeh_to_rcid[eh_long].push_back(rcids);
		hexfh_to_subeh[hr.first] = eh_long;// 记录每个面被选做要细分的边
	}

	
}

void HexMeshing_ConformingTessellation::calc_constraint() {
 
	// 计算拓扑上相连的相等约束: 同一个面的对立边的段数应该保持一致
	for (auto& ep : hexeh_to_minsegnum) {
		auto eh = ep.first;
		if (visited.count(eh)) continue;
		// 计算应该与这条边有相同段数的边
		queue<EH> que;
		unordered_set<EH> visited_iter;// 记录本次迭代访问过的边
		que.emplace(eh);
		visited_iter.insert(eh);
		while (!que.empty()) {// 广度优先搜索
			auto cureh = que.front();
			que.pop();
			for (auto& fh : hexmesh.NeighborFh(cureh)) {
				auto& face = hexmesh.faces(fh);
				for (auto& adjeh : face.getEdgeHandle()) {
					if (!hexmesh.isConnected(adjeh, cureh) && !visited_iter.count(adjeh)) {
						que.emplace(adjeh);
						visited_iter.insert(adjeh);
					}
				}
			}
		}

		for (auto& cureh : visited_iter) {
			if (hexeh_to_minsegnum.count(cureh)) {
				visited.insert(cureh);
				constraint_equal_ehs[eh].insert(cureh);
			}
		}

	}

	// 通过中轴骨架曲线确定应该相等的光线约束
	constraint_equal_rcs.resize(ray_constraints.size());
	for (int i = 0; i < ray_constraints.size(); ++i) {
		constraint_equal_rcs[i] = i;
	}

	unordered_map<VH, unordered_set<EH>> visited_sklvh_skleh;
	for (auto& beg_vh : branch_sklvh) {
		for (auto& beg_eh : sklmesh.NeighborEh(beg_vh)) {
			if (!curve_skel_ehs.count(beg_eh) || visited_sklvh_skleh[beg_vh].count(beg_eh)) continue;
			bool is_branch_to_branch = true;
			auto edge = sklmesh.edges(beg_eh);
			vector<VH> curve_sklvhs;
			curve_sklvhs.emplace_back(beg_vh);
			curve_sklvhs.emplace_back(edge.vh1() == beg_vh ? edge.vh2() : edge.vh1());
			while (1) {
				auto vh = curve_sklvhs.back();
				if (end_sklvh.count(vh)) {
					is_branch_to_branch = false;
					break;
				} else if (branch_sklvh.count(vh)) {
					break;
				}
				for (auto& adjeh : sklmesh.NeighborEh(vh)) {
					edge = sklmesh.edges(adjeh);
					auto adjvh = VH(edge.vh1() + edge.vh2() - vh);
					if (adjvh == curve_sklvhs[curve_sklvhs.size() - 2]) continue;
					curve_sklvhs.emplace_back(adjvh);
					break;
				}
			}
			if (!branch_sklvh.count(curve_sklvhs.back()) && !end_sklvh.count(curve_sklvhs.back())) {
				std::cerr << "Error: We get a invalid curve beacuse its end vertex is not a branch or end node!\n";
				continue;
			}
			
			CurveLine cl;
			cl.beg_vh = beg_vh, cl.beg_eh = beg_eh;
			cl.vertices = curve_sklvhs;
			if (!is_branch_to_branch) curve_lines_b2e.emplace_back(cl);
			else {
				cl.end_vh = curve_sklvhs.back();
				cl.end_eh = sklmesh.getEdgeHandle(curve_sklvhs[curve_sklvhs.size() - 2], curve_sklvhs[curve_sklvhs.size() - 1]);
				curve_lines_b2b.emplace_back(cl);
				visited_sklvh_skleh[cl.end_vh].insert(cl.end_eh);
			} 

		}
	}
	std::cout << "\n[Conforming Tessellation]: We get " << curve_lines_b2b.size() << " curve lines between branching nodes.\n";
	std::cout << "[Conforming Tessellation]: We get " << curve_lines_b2e.size() << " curve lines between branching node and ending node.\n";

	//// 
	//for (auto& cl : curve_lines_b2e) {
	//	auto hexfh = rayid_to_hexfh[sklvh_skleh_to_rayid[cl.beg_vh][cl.beg_eh]];
	//	auto hexehs = hexmesh.faces(hexfh).getEdgeHandle();
	//	auto& ep = hexfh_to_hexeh[hexfh];
	//	for (auto& hexeh : hexehs) {
	//		if (hexeh == ep.first || hexeh == ep.second) {
	//			CurveConstraint cc;
	//			int ccid = curve_constraints.size();
	//			cc.is_b2b = false;
	//			hexeh_sklvh_skleh_to_ccid[hexeh][cl.beg_vh][cl.beg_eh] = ccid;
	//			curve_constraints.emplace_back(cc);
	//		}
	//	}
	//	
	//}

	for (auto& cl : curve_lines_b2b) {

		auto beg_fh = rayid_to_hexfh[sklvh_skleh_to_rayid[cl.beg_vh][cl.beg_eh]];
		auto beg_face = hexmesh.faces(beg_fh);
		auto beg_fvhs = beg_face.getVertexHandle();

		vector<Vex> plane_points[4];
		for (int i = 0; i < 4; ++i) {
			plane_points[i].emplace_back(hexmesh.vertices(beg_fvhs[i]));
		}
		for (int i = 1; i < cl.vertices.size(); ++i) {
			auto v = sklmesh.vertices(cl.vertices[i]);
			Vec dir = (v - sklmesh.vertices(cl.vertices[i - 1])).normalized();// 射线方向
			// 计算平面
			double a, b, c, d;// ax + by + cz + d = 0;
			if (i + 1 == cl.vertices.size()) {// 最后一个顶点
				a = dir.x(), b = dir.y(), c = dir.z();
				d = -(a * v.x() + b * v.y() + c * v.z());

			} else {
				Vec dir2 = (sklmesh.vertices(cl.vertices[i + 1]) - v).normalized();
				if (dir.dot(dir2) < 0) dir2 = dir2 * (-1);
				Vec normal = (dir + dir2).normalized();
				a = normal.x(), b = normal.y(), c = normal.z();
				d = -(a * v.x() + b * v.y() + c * v.z());

			}

			for (int j = 0; j < 4; ++j) {
				// 求出射线与平面的交点
				auto& pre_v = plane_points[j].back();
				double t = -(a * pre_v.x() + b * pre_v.y() + c * pre_v.z() + d)
					/ (a * dir.x() + b * dir.y() + c * dir.z());
				Vex intersection_v = pre_v + dir * t;
				plane_points[j].emplace_back(intersection_v);
			}
		}

		auto end_fh = rayid_to_hexfh[sklvh_skleh_to_rayid[cl.end_vh][cl.end_eh]];
		auto end_face = hexmesh.faces(end_fh);
		auto end_fvhs = end_face.getVertexHandle();
		vector<Vex> rear_points;
		for (int i = 0; i < 4; ++i) {
			rear_points.emplace_back(plane_points[i].back());
		}
		min_dist(rear_points, end_fvhs);
		//min_dir(rear_points, end_fvhs, (hexmesh.getFaceCenter(end_fh) - hexmesh.getFaceCenter(beg_fh)).normalized());
		// 记录配准信息
		cl.hexvhs_src = beg_fvhs;
		cl.hexvhs_end = end_fvhs;
		for (int i = 0; i < 4; ++i) {
			matched_vex.push_back(pair<Vex, Vex>(hexmesh.vertices(beg_fvhs[i]), hexmesh.vertices(end_fvhs[i])));
		}

		// 每个曲线段都对应两个决策变量: 两个方向
		vector<EH> beg_fehs = {
			hexmesh.getEdgeHandle(beg_fvhs[0], beg_fvhs[1]), hexmesh.getEdgeHandle(beg_fvhs[1], beg_fvhs[2]),
			hexmesh.getEdgeHandle(beg_fvhs[2], beg_fvhs[3]), hexmesh.getEdgeHandle(beg_fvhs[3], beg_fvhs[0])
		};
		vector<EH> end_fehs = {
			hexmesh.getEdgeHandle(end_fvhs[0], end_fvhs[1]), hexmesh.getEdgeHandle(end_fvhs[1], end_fvhs[2]),
			hexmesh.getEdgeHandle(end_fvhs[2], end_fvhs[3]), hexmesh.getEdgeHandle(end_fvhs[3], end_fvhs[0])
		};

		unordered_set<EH> decision_ehs;
		decision_ehs.insert(hexfh_to_hexeh[beg_fh].first);
		decision_ehs.insert(hexfh_to_hexeh[beg_fh].second);
		decision_ehs.insert(hexfh_to_hexeh[end_fh].first);
		decision_ehs.insert(hexfh_to_hexeh[end_fh].second);

		for (int i = 0; i < 4; ++i) {
			if (decision_ehs.count(beg_fehs[i])) {
				int j = i;
				if (!decision_ehs.count(end_fehs[j])) {
					j = (i + 2) % 4;
				}
				if (!decision_ehs.count(end_fehs[j])) {
					std::cerr << "Error: We can not find valid decision ending edge handle!!!\n";
					continue;
				}
				/*int ccid = curve_constraints.size();
				CurveConstraint cc;
				cc.is_b2b = true;
				hexeh_sklvh_skleh_to_ccid[beg_fehs[i]][cl.beg_vh][cl.beg_eh] = ccid;
				hexeh_sklvh_skleh_to_ccid[end_fehs[j]][cl.end_vh][cl.end_eh] = ccid;
				curve_constraints.emplace_back(cc);*/

				// 1. 首先找到对应的两条光线
				int rid1 = sklvh_skleh_to_rayid[cl.beg_vh][cl.beg_eh];
				int rid2 = sklvh_skleh_to_rayid[cl.end_vh][cl.end_eh];
				// 2. 确定匹配的变量
				int rcid1 = rid1 * 2;
				if (ray_constraints[rcid1].hexeh != beg_fehs[i]) rcid1++;
				int rcid2 = rid2 * 2;
				if (ray_constraints[rcid2].hexeh != end_fehs[j]) rcid2++;
				if (ray_constraints[rcid1].hexeh != beg_fehs[i] || ray_constraints[rcid2].hexeh != end_fehs[j]) {
					std::cerr << "Error: ray_constraints[rcid1].hexeh != beg_fehs[i] || ray_constraints[rcid2].hexeh != end_fehs[j] ...\n";
					continue;
				}
				constraint_equal_rcs[rcid1] = rcid2;
				hexeh_pairs.emplace_back(hexmesh.getEdgeMidpoint(beg_fehs[i]), hexmesh.getEdgeMidpoint(end_fehs[j]));
			}

		}
		hexfh_pairs.emplace_back(beg_fh, end_fh);

	}
	

}

void HexMeshing_ConformingTessellation::opt_conforming_split() {

	time_t beg_gurobi = clock();
	std::cout << "\n[Gurobi] begin...\n\n";

	int sz_ehs = hexeh_to_minsegnum.size(), sz_rcs = ray_constraints.size();
	unordered_map<EH, int> eh2vid;

	try {

		// Create an environment
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", "mip1_gurobi.log");
		env.set(GRB_DoubleParam::GRB_DoubleParam_TimeLimit, 300);// 设置最长优化时间为 300秒
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);

		// Create variables
		vector<GRBVar> vec_decision;// 每条边的段数
		vec_decision.reserve(sz_ehs);
		for (auto& ep : hexeh_to_minsegnum) {
			int vid = vec_decision.size();
			eh2vid[ep.first] = vid;
			GRBVar var = model.addVar(ep.second, 9999, NULL, GRB_INTEGER);// 决策变量
			vec_decision.emplace_back(var);
		}

		vector<GRBVar> vec_auxvar;// 每段曲线对应的六面体边的段数
		vec_auxvar.reserve(sz_rcs);
		for (int i = 0; i < sz_rcs; ++i) {
			int vid = vec_auxvar.size();
			GRBVar var = model.addVar(1, 9999, NULL, GRB_INTEGER);// 曲线段引入的辅助变量
			vec_auxvar.emplace_back(var);
		}
		

		GRBVar Z = model.addVar(1, INT_MAX, 0, GRB_INTEGER, "Z");

		// Set objective: 最小化边的段数和
		GRBLinExpr obj = Z;
		model.setObjective(obj, GRB_MINIMIZE);

		// 拓扑对立的 边 的段数 应相等
		for (auto& ep : constraint_equal_ehs) {
			auto& var1 = vec_decision[eh2vid[ep.first]];
			for (auto& eh2 : ep.second) {
				auto vid2 = eh2vid[eh2];
				auto& var2 = vec_decision[vid2];
				model.addConstr(var1 == var2);
			}
		}
		std::cout << "Set the constraints of equal ehs success.\n";

		// 相同中轴骨架曲线的边段数应该相等
		for (int rcid = 0; rcid < sz_rcs; ++rcid) {
			if (constraint_equal_rcs[rcid] != rcid) {
				model.addConstr(vec_auxvar[rcid] == vec_auxvar[constraint_equal_rcs[rcid]]);
			}
		}

		// 
		for (auto& ep : hexeh_to_rcid) {
			
			for (auto& rcids : ep.second) {
				GRBLinExpr le = 0;
				for (auto& rcid : rcids) {
					le += vec_auxvar[rcid];
				}
				model.addConstr(vec_decision[eh2vid[ep.first]] == le);
			}
		}
		
		/*for (auto& ep : hexeh_to_minsegnum) {
			auto hexeh = ep.first;
			
			for (auto& curve_line : hexeh_to_curve[hexeh]) {
				GRBLinExpr le = 0;
				bool is_valid = true;
				int cnt_le = 0;
				for (auto& skl_h : curve_line) {
					auto sklvh = skl_h.first;
					auto skleh = skl_h.second;
					if (!hexeh_sklvh_skleh_to_ccid.count(hexeh) || !hexeh_sklvh_skleh_to_ccid[hexeh].count(sklvh) || !hexeh_sklvh_skleh_to_ccid[hexeh][sklvh].count(skleh)) {
						std::cout << "Error: Map not count Hex handle or Skel handle!!! Hex_EH = " 
							<< hexeh << ", Skl_VH = " << sklvh << ", Skl_EH = " << skleh << "\n";
						is_valid = false;
						break;
					}
					int ccid = hexeh_sklvh_skleh_to_ccid[hexeh][sklvh][skleh];
					le += vec_auxvar[ccid];
					cnt_le++;
				}
				if (is_valid) {
					model.addConstr(vec_decision[eh2vid[hexeh]] == le);
					std::cout << "Hex_EH = " << hexeh << ", the num of le is " << cnt_le << ".\n";
				}
				
			}
			
		}
		std::cout << "Set the constraints of curve lines success.\n";*/

		GRBLinExpr edge_split_count = 0;// 细分段数应尽可能少
		for (int i = 0; i < sz_ehs; ++i) {
			edge_split_count += vec_decision[i];
		}
		model.addConstr(Z == edge_split_count);

		// Optimize model
		model.optimize();

		cout << "[Gurobi] Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;


		int result_z = 0;
		for (int i = 0; i < sz_ehs; ++i) {
			auto& var = vec_decision[i];
			if (var.get(GRB_DoubleAttr_X) > 0) {
				//is_selected[i] = true;
				result_z += var.get(GRB_DoubleAttr_X);
			}
		}

		for (auto& ep : hexeh_to_minsegnum) {
			auto& var = vec_decision[eh2vid[ep.first]];
			ep.second = var.get(GRB_DoubleAttr_X);
		}

		for (int rcid = 0; rcid < ray_constraints.size(); ++rcid) {
			ray_constraints[rcid].segment = vec_auxvar[rcid].get(GRB_DoubleAttr_X);
		}

		/*for (auto& auxvar : vec_auxvar) {
			int x = int(auxvar.get(GRB_DoubleAttr_X));
			if (x != 1) {
				std::cout << "Result: curve segment = " << x << ".\n";
			}
			 
		}*/
		/*for (int cid = 0; cid < constraint_count_b2b; ++cid) {
			auto& auxvar = vec_auxvar[cid];
			int x = int(auxvar.get(GRB_DoubleAttr_X));
			std::cout << "[Conforming Tessellation]: Segment of Branch to Branch is " << x << ".\n";
		}*/

		std::cout << "[Conforming Tessellation]: We get the minimum num of segments is " << result_z << " .\n";

	} catch (GRBException e) {
		cout << "[Gurobi] Error code = " << e.getErrorCode() << endl << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "[Gurobi] Exception during optimization" << endl;
	}
	time_t end_gurobi = clock();
	std::cout << "\n[Gurobi] end... cost " << int(end_gurobi - beg_gurobi) << " ms. \n\n";

}

void HexMeshing_ConformingTessellation::assign_sweep_face() {

	assign_face_dense_to_origin();
	
	assign_face_to_ray();

}

void HexMeshing_ConformingTessellation::sweep_curve_skel() {

	// 从 分支节点 出发, 向 终端节点 扫掠生成六面体
	for (auto& curve : curve_lines_b2e) {

		const auto& sklvhs = curve.vertices;// 路径上经过的骨架点
		int rid = sklvh_skleh_to_rayid[curve.beg_vh][curve.beg_eh];
		auto& faces_src = rayid_to_sub_fhs[rid];
		auto faces_path = faces_src;
		for (int i = 1; i < sklvhs.size(); ++i) {
			auto sklv = sklmesh.vertices(sklvhs[i]);
			// 计算射向
			auto dir = (sklv - sklmesh.vertices(sklvhs[i - 1])).normalized();
			// 计算切平面
			double a, b, c, d;
			if (i + 1 == sklvhs.size()) {
				a = dir.x(), b = dir.y(), c = dir.z();
				
			} else {
				Vec dir2 = (sklmesh.vertices(sklvhs[i + 1]) - sklv).normalized();
				if (dir.dot(dir2) < 0) dir2 = dir2 * (-1);
				Vec normal = (dir + dir2).normalized();
				a = normal.x(), b = normal.y(), c = normal.z();
				
			}
			d = -(a * sklv.x() + b * sklv.y() + c * sklv.z());
			double radius = sklmesh.radius[sklvhs[i]], r_max = 0;
			// 扫掠整个一层
			unordered_map<VH, VH> srcvh_to_endvh;
			for (auto& fhs : faces_path) {
				for (auto& fh : fhs) {
					auto src_vhs = hexmesh.faces(fh).getVertexHandle();
					auto end_vhs = src_vhs;
					for (int j = 0; j < 4; ++j) {
						if (!srcvh_to_endvh.count(src_vhs[j])) {
							auto srcv = hexmesh.vertices(src_vhs[j]);
							double t = -(a * srcv.x() + b * srcv.y() + c * srcv.z() + d)
								/ (a * dir.x() + b * dir.y() + c * dir.z());
							Vex end_v = srcv + dir * t;
							srcvh_to_endvh[src_vhs[j]] = hexmesh.AddVertex(end_v);
							r_max = max(r_max, (end_v - sklv).norm());
						}
						end_vhs[j] = srcvh_to_endvh[src_vhs[j]];
					}
					auto hexch = hexmesh.AddCell({ src_vhs[0], src_vhs[1], src_vhs[2], src_vhs[3], end_vhs[0], end_vhs[1], end_vhs[2], end_vhs[3] });
					auto& hexc = hexmesh.cells(hexch);
					// 记录下一次扫掠的路径面
					for (auto& hexfh_new : hexc.getFaceHandle()) {
						if (hexfh_new != fh && !hexmesh.isConnected(fh, hexfh_new)) {
							fh = hexfh_new;
							break;
						}
					}
				}
			}
			// 扫掠完后做半径的调整
			double r_scale = radius / r_max;
			for (auto& vp : srcvh_to_endvh) {
				auto& hexv = hexmesh.vertices(vp.second);
				auto diff = hexv - sklv;
				auto dir = diff.normalized();
				hexv = sklv + dir * diff.norm() * r_scale;
			}

		}

	}


	// 从 分支节点 出发, 向 分支节点 扫掠生成六面体
	for (auto& curve : curve_lines_b2b) {

		const auto& sklvhs = curve.vertices;// 路径上经过的骨架点
		int rid_src = sklvh_skleh_to_rayid[curve.beg_vh][curve.beg_eh];
		int rid_end = sklvh_skleh_to_rayid[curve.end_vh][curve.end_eh];
		auto faces_src = rayid_to_sub_fhs[rid_src];
		auto faces_end = rayid_to_sub_fhs[rid_end];
		auto faces_path = faces_src;
		unordered_map<FH, FH> src_to_rear;// 记录起始面到尾面的映射
		for (auto& fhs : faces_src) {
			for (auto& fh : fhs) {
				src_to_rear[fh] = fh;
			}
		}
		// 扫掠到分支节点前
		for (int i = 1; i + 1 < sklvhs.size(); ++i) {
			auto sklv = sklmesh.vertices(sklvhs[i]);
			// 计算射向
			auto dir = (sklv - sklmesh.vertices(sklvhs[i - 1])).normalized();
			// 计算切平面
			double a, b, c, d;
			Vec dir2 = (sklmesh.vertices(sklvhs[i + 1]) - sklv).normalized();
			if (dir.dot(dir2) < 0) dir2 = dir2 * (-1);
			Vec normal = (dir + dir2).normalized();
			a = normal.x(), b = normal.y(), c = normal.z();
			d = -(a * sklv.x() + b * sklv.y() + c * sklv.z());
			double radius = sklmesh.radius[sklvhs[i]], r_max = 0;
			// 扫掠整个一层
			unordered_map<VH, VH> srcvh_to_endvh;
			unordered_map<FH, FH> srcfh_to_endfh;
			for (auto& fhs : faces_path) {
				for (auto& fh : fhs) {
					auto src_vhs = hexmesh.faces(fh).getVertexHandle();
					auto end_vhs = src_vhs;
					for (int j = 0; j < 4; ++j) {
						if (!srcvh_to_endvh.count(src_vhs[j])) {
							auto srcv = hexmesh.vertices(src_vhs[j]);
							double t = -(a * srcv.x() + b * srcv.y() + c * srcv.z() + d)
								/ (a * dir.x() + b * dir.y() + c * dir.z());
							Vex end_v = srcv + dir * t;
							srcvh_to_endvh[src_vhs[j]] = hexmesh.AddVertex(end_v);
							r_max = max(r_max, (end_v - sklv).norm());
						}
						end_vhs[j] = srcvh_to_endvh[src_vhs[j]];
					}
					auto hexch = hexmesh.AddCell({ src_vhs[0], src_vhs[1], src_vhs[2], src_vhs[3], end_vhs[0], end_vhs[1], end_vhs[2], end_vhs[3] });
					auto& hexc = hexmesh.cells(hexch);
					// 记录下一次扫掠的路径面
					for (auto& hexfh_new : hexc.getFaceHandle()) {
						if (hexfh_new != fh && !hexmesh.isConnected(fh, hexfh_new)) {
							srcfh_to_endfh[fh] = hexfh_new;
							fh = hexfh_new;
							break;
						}
					}
				}
			}
			// 扫掠完后做半径的调整
			double r_scale = radius / r_max;
			for (auto& vp : srcvh_to_endvh) {
				auto& hexv = hexmesh.vertices(vp.second);
				auto diff = hexv - sklv;
				auto dir = diff.normalized();
				hexv = sklv + dir * diff.norm() * r_scale;
			}
			// 扫掠完后更新源面顶点与路径顶点的映射
			for (auto& vp : src_to_rear) {
				if (!srcfh_to_endfh.count(vp.second)) {
					std::cerr << "\nError: srcvh_to_endvh not count vp.seocnd!!!\n";
				}
				vp.second = srcfh_to_endfh[vp.second];

			}

		}

		// 将源面与目标面匹配
		// 调整面片顺序
		//auto hexvhs_src = curve.hexvhs_src, hexvhs_end = curve.hexvhs_end;
		unordered_map<VH, VH> hexvh_matched;
		for (int i = 0; i < 4; ++i) {
			hexvh_matched[curve.hexvhs_src[i]] = curve.hexvhs_end[i];
		}
		auto hexvhs_src_ordered = hexfh_origin_to_sub_order[rayid_to_hexfh[rid_src]];
		auto hexvhs_end_ordered = hexfh_origin_to_sub_order[rayid_to_hexfh[rid_end]];

		for (int i = 0; i < 4; ++i) {
			if (hexvhs_end_ordered[i] == hexvh_matched[hexvhs_src_ordered[0]]) {
				int j = (i + 1) % 4;
				if (hexvhs_end_ordered[j] != hexvh_matched[hexvhs_src_ordered[1]]) {// 需要反转
					reverse(hexvhs_end_ordered.begin(), hexvhs_end_ordered.end());
					// 行互换即可
					int left_r = 0, right_r = faces_end.size() - 1;
					while (left_r < right_r) {
						vector<FH> tmp = faces_end[left_r];
						faces_end[left_r] = faces_end[right_r];
						faces_end[right_r] = tmp;
						left_r++;
						right_r--;
					}
				}
				break;
			}
		}
		
		for (int i = 0; i < 4; ++i) {
			if (hexvhs_end_ordered[i] == hexvh_matched[hexvhs_src_ordered[0]]) {
				int j = (i + 1) % 4;
				if (hexvhs_end_ordered[j] != hexvh_matched[hexvhs_src_ordered[1]]) {
					std::cout << "\033[31m" << "Error: Hex VH is not matching!!!\n" << "\033[37m";
					break;
				}
				// 需要把 i 的位置移到 0
				if (i == 0) break;
				std::cout << "\033[31m" << "i = " << i << "\n" << "\033[37m";
				int r_n = faces_end.size(), c_n = faces_end[0].size();

				if (i == 1) {// 向左旋转 90 度
					vector<vector<FH>> faces_end_new(c_n, vector<FH>(r_n));
					for (int r = 0; r < r_n; ++r) {
						for (int c = 0; c < c_n; ++c) {
							faces_end_new[c_n - 1 - c][r] = faces_end[r][c];
						}
					}
					faces_end = faces_end_new;

				} else if (i == 2) {
					// 1. 行互换
					int left_r = 0, right_r = faces_end.size() - 1;
					while (left_r < right_r) {
						vector<FH> tmp = faces_end[left_r];
						faces_end[left_r] = faces_end[right_r];
						faces_end[right_r] = tmp;
						left_r++;
						right_r--;
					}
					// 2. 列反转
					for (auto& fhs : faces_end) {
						reverse(fhs.begin(), fhs.end());
					}

				} else if (i == 3) {// 向右旋转 90 度
					vector<vector<FH>> faces_end_new(c_n, vector<FH>(r_n));
					for (int r = 0; r < r_n; ++r) {
						for (int c = 0; c < c_n; ++c) {
							faces_end_new[c][r_n - 1 - r] = faces_end[r][c];
						}
					}
					faces_end = faces_end_new;

				}

			}

		}

		if (faces_src.size() == faces_end.size() && faces_src[0].size() == faces_end[0].size()) {
			int rn = faces_src.size(), cn = faces_src[0].size();
			for (int r = 0; r < rn; ++r) {
				for (int c = 0; c < cn; ++c) {
					auto fh1 = src_to_rear[faces_src[r][c]];
					auto fh2 = faces_end[r][c];
					VolumeEditer::add_cell_disconnected(hexmesh, fh1, fh2);
					//HexAddCell::add_cell_2fh(hexmesh, fh1, fh2);
					std::cout << "Add Cell success!!!\n";
				}
			}

		} else {
			std::cout << "\033[31m" << "faces_src.size() = " << faces_src.size() << ", faces_src[0].size() = " << faces_src[0].size()
				<< "; faces_end.size() = " << faces_end.size() << ", faces_end[0].size() = " << faces_end[0].size() << ".\n"
				<< "\tError: Hex VH is not matching!!!\n" << "\033[37m";

		}

	}


}

void HexMeshing_ConformingTessellation::assign_face_dense_to_origin() {

	// 为每个输入六面体面匹配新的面
	for (auto& hr : hexfh_to_rayid) {
		auto hexfh_origin = hr.first;
		auto rids = hr.second;
		auto face_origin = hexmesh_origin.faces(hexfh_origin);
		auto ehs_origin = face_origin.getEdgeHandle();
		auto vhs_ordered = face_origin.getVertexHandle();
		// 根据选择的细分边调整顺序
		auto subeh = hexfh_to_subeh[hr.first];
		auto sube = hexmesh_origin.edges(subeh);
		for (int i = 0; i < 4; ++i) {
			if (vhs_ordered[i] == sube.vh1()) {
				if (vhs_ordered[(i + 3) % 4] == sube.vh2()) {
					reverse(vhs_ordered.begin(), vhs_ordered.end());
				}
				break;
			}
		}
		for (int i = 0; i < 4; ++i) {
			if (vhs_ordered[i] == sube.vh1()) {
				auto vhs_origin_new = vhs_ordered;
				for (int j = 0; j < 4; ++j) {
					vhs_origin_new[j] = vhs_ordered[(i + j) % 4];
				}
				vhs_ordered = vhs_origin_new;
			}
		}
		// 检查顺序是否合理
		std::cout << "\nSub Edge: vh1 = " << sube.vh1() << ", vh2 = " << sube.vh2() <<
			"; vhs: " << vhs_ordered[0] << "->" << vhs_ordered[1] << "->" << vhs_ordered[2] << "->" << vhs_ordered[3] << ".\n";

		// 1. 首先使用最短路径, 获取所有的顶点 handle
		vector<vector<VH>> path_dense_vhs;// 0->1, 1->2, 3->2, 0->3
		for (int i = 0; i < 4; ++i) {
			auto vh_src = vhs_ordered[i], vh_end = vhs_ordered[(i + 1) % 4];
			if (i == 2 || i == 3) {// 更换顺序, 便于下一步组织面片
				swap(vh_src, vh_end);
			}
			unordered_map<VH, VH> vh_prev;// 记录前缀
			queue<VH> que;
			que.emplace(vh_src);
			vh_prev[vh_src] = VH(-1);
			while (!que.empty()) {
				auto vh_cur = que.front();
				que.pop();
				if (vh_cur == vh_end) break;
				for (auto& adjvh : hexmesh.NeighborVh(vh_cur)) {
					if (vh_prev.count(adjvh)) continue;
					que.emplace(adjvh);
					vh_prev[adjvh] = vh_cur;
				}
			}
			if (!vh_prev.count(vh_end)) {
				std::cerr << "Error: vh_src is not connected with vh_end!!!\n";
				continue;
			}

			VH vh_cur = vh_end;
			vector<VH> path;
			while (vh_prev.count(vh_cur)) {
				path.emplace_back(vh_cur);
				vh_cur = vh_prev[vh_cur];
			}
			reverse(path.begin(), path.end());
			path_dense_vhs.emplace_back(path);


		}

		int rows = path_dense_vhs[0].size(), cols = path_dense_vhs[1].size();
		vector<VH> path_last_top = path_dense_vhs[0];
		vector<vector<FH>> fhs_dense;
		
		for (int c = 1; c < cols; ++c) {
			int r = 1;
			vector<VH> path_cur_btm = { path_dense_vhs[3][c] };
			vector<FH> fhs_row;
			while (r < rows) {
				VH vh_left_top = path_last_top[r - 1];
				VH vh_right_top = path_last_top[r];
				VH vh_left_btm = path_cur_btm.back();
				FH fh_common = getFaceHandle(hexmesh, vh_left_top, vh_right_top, vh_left_btm);
				if (fh_common == -1) return;

				auto fvhs_tmp = hexmesh.faces(fh_common).getVertexHandle();
				VH vh_right_btm(fvhs_tmp[0] + fvhs_tmp[1] + fvhs_tmp[2] + fvhs_tmp[3] - vh_left_top - vh_right_top - vh_left_btm);
				fhs_row.emplace_back(fh_common);
				path_cur_btm.emplace_back(vh_right_btm);
				r++;
			}

			path_last_top = path_cur_btm;
			fhs_dense.emplace_back(fhs_row);

		}

		hexfh_origin_to_current[hexfh_origin] = fhs_dense;
		hexfh_origin_to_sub_order[hexfh_origin] = vhs_ordered;

		// 调试用: 绘制原先边对应的线
		for (auto& path_vhs : path_dense_vhs) {
			vector<Vex> path_vs;
			for (auto& vh : path_vhs) {
				auto& v = hexmesh.vertices(vh);
				path_vs.emplace_back(v);
			}
			path_dense_vex.push_back(path_vs);

		}

		// 合法性检验, 检查每个面对应的段数是否合法
		int fcnt_actually = fhs_dense.size() * fhs_dense[0].size();
		int fcnt_optimized = hexeh_to_minsegnum[hexfh_to_hexeh[hexfh_origin].first] * hexeh_to_minsegnum[hexfh_to_hexeh[hexfh_origin].second];
		if (fcnt_actually != fcnt_optimized) {
			std::cout << "\nError: face count between actually and optimized is not equal!!!\n";
		}

	}

}

void HexMeshing_ConformingTessellation::assign_face_to_ray() {

	// 为每条光线(中轴曲线)分配相应的扫掠面
	rayid_to_sub_fhs.resize(rays.size());

	for (auto& hr : hexfh_to_rayid) {

		auto hexfh_origin = hr.first;
		auto rids = hr.second;
		
		const auto subeh = hexfh_to_subeh[hexfh_origin];
		auto fhs_dense = hexfh_origin_to_current[hr.first];
		auto fhs_order = hexfh_origin_to_sub_order[hr.first];

		vector<pair<int, double>> rid_len;// 根据各交点在该方向上的投影决定分配顺序
		auto sube = hexmesh_origin.edges(subeh);
		auto sub_dir = (hexmesh_origin.vertices(sube.vh2()) - hexmesh_origin.vertices(sube.vh1())).normalized();
		auto sub_src = hexmesh_origin.vertices(sube.vh1());
		for (auto& rid : rids) {
			auto vex = rayid_to_intersection[rid];
			auto vec = vex - sub_src;
			double len_project = vec.dot(sub_dir);
			rid_len.push_back(pair<int, double>(rid, len_project));
		}
		// 按投影长度升序排序
		sort(rid_len.begin(), rid_len.end(), [&](pair<int, double>& p1, pair<int, double>& p2) {
			return p1.second < p2.second;
			});
		// 检查合法性
		vector<int> num_of_segment;
		for (auto& rp : rid_len) {
			int rid = rp.first, rcid = rid * 2;
			if (ray_constraints[rcid].hexeh != subeh) {
				rcid++;
				if (ray_constraints[rcid].hexeh != subeh) {
					std::cerr << "\nError: Ray_Constraint's eh is not a sub eh!!!\n";
					continue;
				}
			}
			
			num_of_segment.emplace_back(ray_constraints[rcid].segment);
		}

		int sum_of_ns = accumulate(num_of_segment.begin(), num_of_segment.end(), 0);

		if (sum_of_ns != fhs_dense[0].size()) {
			std::cerr << "\nError: Sum of segment is not equal with face dimension!!!\n";
			std::cout << "Sum_of_ns = " << sum_of_ns << ", fhs_dense.sz = " << fhs_dense.size() << ", fhs_dense[0].sz = " << fhs_dense[0].size() << "\n";
			continue;
		} else {
			std::cout << "OK!!! Sum of segment is equal with face dimension!!!\n";
			std::cout << "Sum_of_ns = " << sum_of_ns << ", fhs_dense.sz = " << fhs_dense.size() << ", fhs_dense[0].sz = " << fhs_dense[0].size() << "\n";
		}
		int col_idx = 0;
		for (int i = 0; i < num_of_segment.size(); ++i) {
			int rid = rid_len[i].first, scnt = num_of_segment[i];
			vector<vector<FH>> ray_fhs;
			for (int j = 0; j < fhs_dense.size(); ++j) {
				vector<FH> ray_fhs_r;
				for (int k = 0; k < scnt; ++k) {
					ray_fhs_r.push_back(fhs_dense[j][k + col_idx]);
				}
				ray_fhs.push_back(ray_fhs_r);
			}
			col_idx += scnt;
			rayid_to_sub_fhs[rid] = ray_fhs;
		}

	}

}

void HexMeshing_ConformingTessellation::min_dist(vector<Vex> rear_vex, vector<VH>& end_vhs) {

	vector<vector<VH>> modes = {
		{ end_vhs[0], end_vhs[1], end_vhs[2], end_vhs[3] }, { end_vhs[0], end_vhs[3], end_vhs[2], end_vhs[1] }, 
		{ end_vhs[1], end_vhs[2], end_vhs[3], end_vhs[0] }, { end_vhs[1], end_vhs[0], end_vhs[3], end_vhs[2] },
		{ end_vhs[2], end_vhs[3], end_vhs[0], end_vhs[1] }, { end_vhs[2], end_vhs[1], end_vhs[0], end_vhs[3] },
		{ end_vhs[3], end_vhs[0], end_vhs[1], end_vhs[2] }, { end_vhs[3], end_vhs[2], end_vhs[1], end_vhs[0] }
	};
	double sum_norm2_min = 99999999.0;
	int mode = -1;
	for (int i = 0; i < modes.size(); ++i) {
		double sum_norm2 = 0;
		for (int j = 0; j < 4; ++j) {
			sum_norm2 += (rear_vex[j] - hexmesh.vertices(modes[i][j])).norm2();
		}
		if (sum_norm2 < sum_norm2_min) {
			sum_norm2_min = sum_norm2;
			mode = i;
		}
	}

	end_vhs = modes[mode];

}

void HexMeshing_ConformingTessellation::min_dir(vector<Vex> rear_vex, vector<VH>& end_vhs, Vec dir) {

	vector<vector<VH>> modes = {
		{ end_vhs[0], end_vhs[1], end_vhs[2], end_vhs[3] }, { end_vhs[0], end_vhs[3], end_vhs[2], end_vhs[1] },
		{ end_vhs[1], end_vhs[2], end_vhs[3], end_vhs[0] }, { end_vhs[1], end_vhs[0], end_vhs[3], end_vhs[2] },
		{ end_vhs[2], end_vhs[3], end_vhs[0], end_vhs[1] }, { end_vhs[2], end_vhs[1], end_vhs[0], end_vhs[3] },
		{ end_vhs[3], end_vhs[0], end_vhs[1], end_vhs[2] }, { end_vhs[3], end_vhs[2], end_vhs[1], end_vhs[0] }
	};
	double sum_abs_dot_max = 0;
	int mode = -1;
	for (int i = 0; i < modes.size(); ++i) {
		double sum_abs_dot = 0;
		for (int j = 0; j < 4; ++j) {
			Vec vec = (hexmesh.vertices(modes[i][j]) - rear_vex[j]).normalized();
			sum_abs_dot = abs(dir.dot(vec));
		}
		if (sum_abs_dot > sum_abs_dot_max) {
			sum_abs_dot_max = sum_abs_dot;
			mode = i;
		}
	}
	end_vhs = modes[mode];

}

bool HexMeshing_ConformingTessellation::is_corner(VH vh) {

	auto& adjehs = sklmesh.NeighborEh(vh);
	auto& v = sklmesh.vertices(vh);
	vector<Vec> dirs;
	for (auto& adjeh : adjehs) {
		auto& e = sklmesh.edges(adjeh);
		VH adjvh(e.vh1() + e.vh2() - vh);
		auto& adjv = sklmesh.vertices(adjvh);
		Vec vec = (adjv - v).normalized();
		dirs.emplace_back(vec);
	}

	int sz = dirs.size();
	for (int i = 0; i < sz; ++i) {
		for (int j = i + 1; j < sz; ++j) {
			if (dirs[i].dot(dirs[j]) > std::cos(120.0 / 180.0 * M_PI)) return true;
		}
	}

	return false;

}

void HexMeshing_ConformingTessellation::recast_rays() {

	hexfh_to_rayid.clear();

	for (int i = 0; i < rays.size(); ++i) {
		recast_ray(i);
	}

}

void HexMeshing_ConformingTessellation::recast_ray(int rayid, int sklvh, int skleh) {

	if (rayid >= rays.size()) {
		std::cerr << "Error: Ray ID is invalid!\n";
		return;
	}

	const auto& ray = rays[rayid];

	auto intersection = bvh_tree->getIntersection(ray);
	if (intersection.happened == false) {
		std::cerr << "Ray no intersect with hex mesh.\n";
		return;
	}
	FH fh(intersection.index);
	rayid_to_hexfh.emplace_back(fh);
	rayid_to_intersection.push_back(Vex(intersection.pos[0], intersection.pos[1], intersection.pos[2]));
	hexfh_to_rayid[fh].emplace_back(rayid);
	/*if (vh != -1 && eh != -1) {
		hexfh_to_sklvh_skleh[fh].push_back(pair<VH, EH>(VH(vh), EH(eh)));
	}*/

}

bool HexMeshing_ConformingTessellation::dfs_curve_vhs_between_branch(std::vector<VH>& vhs) {
	//std::cerr << "We are dfs curve vhs.\n";
	if (branch_sklvh.count(vhs.back())) {// 只统计分支节点与分支节点之间的曲线路径
		return true;
	} else if (end_sklvh.count(vhs.back())) {
		return false;
	}
	VH pre = vhs[vhs.size() - 2], cur = vhs.back();// 当前结点只会是 joint node
	for (auto& eh : sklmesh.NeighborEh(cur)) {
		if (curve_skel_ehs.count(eh)) {
			auto& e = sklmesh.edges(eh);
			auto next = VH(e.vh1() + e.vh2() - cur);
			if (next == pre) continue;
			vhs.push_back(next);
			return dfs_curve_vhs_between_branch(vhs);
		}
	}
	return false;
}

FH HexMeshing_ConformingTessellation::getFaceHandle(MeshKernel::VolumeMesh& mesh, EH eh1, EH eh2) {

	FH res(-1);
	if (!mesh.isValid(eh1) || !mesh.isValid(eh2)) {
		std::cout << "Error: Func_Name = getFaceHandle, exist invalid edge handle!!!\n";
		return res;
	}

	auto& fhs1 = mesh.NeighborFh(eh1);
	auto& fhs2 = mesh.NeighborFh(eh2);

	for (auto& fh : fhs1) {
		if (fhs2.count(fh)) {
			if (res != -1) {
				std::cerr << "Error : Func_Name = getFaceHandle, exist more than one common face!!!\n";
			}
			res = fh;
		}
	}

	if (res == -1) {
		std::cerr << "Error : Func_Name = getFaceHandle, not exist common face!!!\n";
	}

	return res;

}

FH HexMeshing_ConformingTessellation::getFaceHandle(MeshKernel::VolumeMesh& mesh, VH vh1, VH vh2, VH vh3) {

	FH res(-1);
	if (!mesh.isValid(vh1) || !mesh.isValid(vh2) || !mesh.isValid(vh3)) {
		std::cout << "Error: Func_Name = getFaceHandle, exist invalid vertex handle!!!\n";
		return res;
	}

	auto& fhs1 = mesh.NeighborFh(vh1);
	auto& fhs2 = mesh.NeighborFh(vh2);
	auto& fhs3 = mesh.NeighborFh(vh3);

	for (auto& fh : fhs1) {
		if (fhs2.count(fh) && fhs3.count(fh)) {
			if (res != -1) {
				std::cerr << "Error : Func_Name = getFaceHandle, exist more than one common face!!!\n";
			}
			res = fh;
		}
	}

	if (res == -1) {
		std::cerr << "Error : Func_Name = getFaceHandle, not exist common face!!!\n";
		std::cerr << "\t vh1 = " << vh1 << ", vh2 = " << vh2 << ", vh3 = " << vh3 << ".\n";
		debug_should_be_com_faces.emplace_back(vh1);
		debug_should_be_com_faces.emplace_back(vh2);
		debug_should_be_com_faces.emplace_back(vh3);
	}

	return res;

}

void HexMeshing_ConformingTessellation::cone_detection() {

	std::unordered_map<FH, double> faces_area;
	for (auto& fp : hexmesh.allfaces()) {
		faces_area[fp.first] = hexmesh.getQuadArea(fp.first);
	}

	std::unordered_map<EH, int> edges_need_sub;// 需要细分的边以及相应的细分层数
	for (auto& cp : hexmesh.allcells()) {
		const auto& adjchs = hexmesh.NeighborCh(cp.first);
		if (adjchs.size() > 2) continue;
		const auto& _c_fhs = cp.second.getFaceHandle();
		vector<FH> c_fhs(_c_fhs.begin(), _c_fhs.end());
		for (int i = 0; i < 6; ++i) {
			for (int j = i + 1; j < 6; ++j) {
				bool is_connected_via_v = false;
				const auto& f_vhs0 = hexmesh.faces(c_fhs[i]).getVertexHandle();
				const auto& f_vhs1 = hexmesh.faces(c_fhs[j]).getVertexHandle();
				for (auto& f_vh0 : f_vhs0) {
					if (std::find(f_vhs1.begin(), f_vhs1.end(), f_vh0) != f_vhs1.end()) {
						is_connected_via_v = true;
						break;
					}
				}
				if (is_connected_via_v) continue;// 不是相对立的面
				// 判断是否是两个暴露在外面的面, 若是, 则不细分
				auto f_chs0 = hexmesh.NeighborCh(c_fhs[i]);
				const auto& f_chs1 = hexmesh.NeighborCh(c_fhs[j]);
				f_chs0.insert(f_chs1.begin(), f_chs1.end());
				if (f_chs0.size() < 2) continue;
				double ratio = faces_area[c_fhs[i]] / faces_area[c_fhs[j]];
				if (ratio < 1.f) ratio = 1.0 / ratio;
				if (ratio >= 4.0) {
					int num = (ratio / 4.0) + 1;
					EH eh(-1);
					for (auto& f_vh0 : f_vhs0) {
						for (auto& f_vh1 : f_vhs1) {
							if (hexmesh.isConnected(f_vh0, f_vh1)) {
								eh = hexmesh.getEdgeHandle(f_vh0, f_vh1);
								break;
							}
						}
					}
					edges_need_sub[eh] = std::min(num, 3);
				}
				break;
			}
		}
	}

	std::cout << "Detect " << edges_need_sub.size() << " cones!!!" << std::endl;

	for (auto& ep : edges_need_sub) {
		if (hexmesh.isValid(ep.first)) {
			Subdivision_Edge::subdivision_edges(hexmesh, ep.first, ep.second);
		}
	}

	hexmesh.updateAllHandles();

}