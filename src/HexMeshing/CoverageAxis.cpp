#include "CoverageAxis.h"

CoverageAxis::CoverageAxis(MeshKernel::SurfaceMesh& _mesh) : mesh(_mesh) {

	offset_radius = 0.02;
	init();

}

void CoverageAxis::init() {

	// 1. 初始化 AABB & BVH 树
	std::vector<Vector3f> vertices;
	for (auto& fp : mesh.allfaces()) {
		auto vhs = fp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& v = mesh.vertices(vh);
			vertices.push_back(Vector3f(v.x(), v.y(), v.z()));
		}
	}
	ab_tree = new AABB_Tree(vertices);

	bvh_tree = new BVH_Tree();
	bvh_tree->buildBVH_Tree(mesh);

	// 2. 初始化半径偏移
	mesh.initBBox();
	double diagonal_length = (mesh.BBoxMax - mesh.BBoxMin).norm();
	offset_radius *= diagonal_length;
	//offset_radius *= 1.2;

	// 3. 初始化采样点

	//int it = 0;
	//for (auto& fp : mesh.allfaces()) {// 均匀采样
	//	if (it % 3 == 0) {
	//		Vex center_face = mesh.getFaceCenter(fp.first);
	//		SurfacePoint* sp = new SurfacePoint();
	//		sp->pos = Vector3d(center_face.x(), center_face.y(), center_face.z());
	//		sp->id = surface_points.size();
	//		surface_points.emplace_back(sp);
	//	}
	//	it++;
	//}

	// 随机采样
	int sample_points_count = std::min(3500, int(mesh.fsize() * 1.25));
	map<double, FH> area2fh;
	double area_sum = 0.0;
	for (auto& fp : mesh.allfaces()) {
		const auto& vhs = fp.second.getVertexHandle();
		auto& v0 = mesh.vertices(vhs[0]);
		auto& v1 = mesh.vertices(vhs[1]);
		auto& v2 = mesh.vertices(vhs[2]);
		Vec vec01 = v1 - v0;
		Vec vec02 = v2 - v0;
		double area = (vec01.cross(vec02)).norm();
		area_sum += area;
		area2fh[area_sum] = fp.first;
	}

	while (surface_points.size() < sample_points_count) {
		double area = (rand() % RAND_MAX) * 1.0 / RAND_MAX * area_sum;
		auto it = area2fh.lower_bound(area);
		FH fh = (*it).second;
		if (!mesh.isValid(fh)) continue;
		double wa = 1.0 * (rand() % 10000007);
		double wb = 1.0 * (rand() % 10000007);
		double wc = 1.0 * (rand() % 10000007);
		double ws = wa + wb + wc;
		wa /= ws;
		wb /= ws;
		wc /= ws;
		const auto& vhs = mesh.faces(fh).getVertexHandle();
		Vex point = mesh.vertices(vhs[0]) * wa + mesh.vertices(vhs[1]) * wb + mesh.vertices(vhs[2]) * wc;
		SurfacePoint* sp = new SurfacePoint();
		sp->pos = Vector3d(point.x(), point.y(), point.z());
		sp->id = surface_points.size();
		surface_points.emplace_back(sp);
	}

	std::cout << "[Coverage Axis]: We get " << surface_points.size() << " surface points.\n";

}

void CoverageAxis::compute_coverage_axis() {

	inner_points_generation();

	point_selection_based_on_set_coverage_MILP();

	//point_selection_based_on_set_coverage_heap();

	connection_establishment();
	medial_mesh.writeMaFile("C:/My Files/Graphics/model_data/!test_data/!CoverageAxis_NoPostprocessing.ca.ma");
	postprocessing();
	medial_mesh.writeMaFile("C:/My Files/Graphics/model_data/!test_data/!CoverageAxis_Final.ca.ma");

	is_computed = true;

}

void CoverageAxis::inner_points_generation() {

	// 方式一: 无面信息, 运算较快, 但最终结果可能存在空洞
	vector<Vector3d> sample_points;
	for (auto* sp : surface_points) {
		sample_points.emplace_back(sp->pos);
	}
	VoronoiDiagram::get_medial_mesh(sample_points, medial_mesh);

	//// 方式二: 有面信息
	//vector<Eigen::Vector3d> sample_points;
	//for (auto* sp : surface_points) {
	//	sample_points.emplace_back(Eigen::Vector3d(sp->pos.x(), sp->pos.y(), sp->pos.z()));
	//}
	//VoronoiMesh::get_voronoi_mesh(sample_points, medial_mesh);
	//medial_mesh.updateAllHandles();

	Vector3d dir(1, 1, 1);

	int vcnt = medial_mesh.vsize();
	vector<bool> is_selected(vcnt, true);
#pragma omp parallel for
	for (int i = 0; i < vcnt; ++i) {// 并行版本
		VH vh(i);
		auto& v = medial_mesh.vertices(vh);
		Vector3d pos(v.x(), v.y(), v.z());
		Ray ray(pos, dir);
		auto intersection = bvh_tree->getIntersection(ray);
		if (intersection.intersection_count % 2 == 0) {
			is_selected[i] = false;
		}
	}
	//std::cout << "[Coverage Axis]: mark outer vertices success.\n";

	// 方式一: 直接删除顶点, 可能会引入空洞
	for (int i = 0; i < vcnt; ++i) {
		if (!is_selected[i]) {
			VH vh(i);
			medial_mesh.DeleteVertex(vh);
		}
	}
	medial_mesh.updateAllHandles();

	//// 方式二: 使用边坍缩, 可能会存在大量面折叠
	//connection_establishment(is_selected);

	for (auto& vp : medial_mesh.allvertices()) {
		InnerPoint* point = new InnerPoint();
		point->pos = Vector3d(vp.second.x(), vp.second.y(), vp.second.z());
		point->id = inner_points.size();
		inner_points.emplace_back(point);
		ip_to_vh[point] = vp.first;
		vh_to_ip[vp.first] = point;
	}

	// 初始化内部点的半径
	int inner_pcnt = inner_points.size();
#pragma omp parallel for
	for (int i = 0; i < inner_pcnt; ++i) {
		auto* inner_point = inner_points[i];
		Vector3f pos(inner_point->pos.x(), inner_point->pos.y(), inner_point->pos.z()), np;
		double radius = ab_tree->findNearstPoint(pos, np);
		inner_point->radius = radius + offset_radius;// 膨胀后的半径
		medial_mesh.radius[ip_to_vh[inner_point]] = radius;// 未膨胀的半径
	}

	std::cout << "[Coverage Axis]: We get " << medial_mesh.vsize() << " inner points.\n";
}

void CoverageAxis::point_selection_based_on_set_coverage_heap() {

#pragma omp parallel for
	for (int i = 0; i < inner_points.size(); ++i) {// 并行执行
		auto* inner_point = inner_points[i];
		inner_point->coverage_node = new CoverageNode();
		auto* node = inner_point->coverage_node;
		node->centroid = inner_point;
		for (auto* surface_point : surface_points) {
			double dist = (surface_point->pos - inner_point->pos).Length();
			if (dist <= inner_point->radius) {
				node->coverage_points.insert(surface_point);
			}
		}
	}
	for (auto* inner_point : inner_points) {
		heap.push(inner_point->coverage_node);
		for (auto* surface_point : inner_point->coverage_node->coverage_points) {
			surface_point->coverage_nodes.emplace_back(inner_point->coverage_node);
		}
	}

	std::cout << "Initialize heap success. heap size = " << heap.size() << "\n";

	//unordered_set<SurfacePoint*> coveraged_surface_points;
	int iter = 0;
	while (!heap.empty() && surface_points_coveraged.size() < surface_points.size()) {

		auto* coverage_node = heap.top();
		heap.pop();

		if (coverage_node->is_valid == false) {
			unordered_set<SurfacePoint*> srf_ps = coverage_node->coverage_points;
			for (auto* sp : srf_ps) {
				if (surface_points_coveraged.count(sp)) {
					coverage_node->coverage_points.erase(sp);
				}
			}
			if (coverage_node->coverage_points.size() > 0) {
				coverage_node->is_valid = true;
				heap.push(coverage_node);
			}
			continue;
		}

		//if (coverage_node->coverage_points.size() <= 3) break;

		auto* inner_point = coverage_node->centroid;

		if (inner_points_selected.count(inner_point)) {
			std::cerr << "Excounter same inner point\n";
			break;
		}

		inner_points_selected.insert(inner_point);
		vhs_selected.insert(ip_to_vh[inner_point]);

		//std::cout << "Iter " << iter++ << ": inner point id = " << inner_point->id << ", coveraged points size increase " << coverage_node->coverage_points.size() << std::endl;

		for (auto* surface_point : coverage_node->coverage_points) {
			if (surface_points_coveraged.count(surface_point)) continue;
			for (auto* c_node : surface_point->coverage_nodes) {
				c_node->is_valid = false;// 旧的无效的
			}
			surface_point->coverage_nodes.clear();
			surface_points_coveraged.insert(surface_point);
		}

		delete coverage_node;
		coverage_node = nullptr;


	}

	std::cout << "We selecte " << inner_points_selected.size() << " points coverage "
		<< surface_points_coveraged.size() << " sample points.\n";

}

void CoverageAxis::point_selection_based_on_set_coverage_MILP() {

	int m = surface_points.size(), n = inner_points.size();

	vector<vector<int>> is_coveraged(m, vector<int>(n, 0));// 0 表示未被覆盖, 1 表示被覆盖

#pragma omp parallel for
	for (int i = 0; i < inner_points.size(); ++i) {// 并行执行
		auto* inner_point = inner_points[i];
		for (auto* surface_point : surface_points) {
			double dist = (surface_point->pos - inner_point->pos).Length();
			if (dist <= inner_point->radius) {
				is_coveraged[surface_point->id][inner_point->id] = 1;
			}
		}
	}

	for (int i = 0; i < m; ++i) {
		int j = 0;
		while (j < n) {
			if (is_coveraged[i][j] == 1) break;
			j++;
		}
		if (j == n) {
			std::cerr << "Error: Exist surface points are not coveraged! So model is infeasible.\n";
			return;
		}
	}

	/*std::ofstream off("C:/My Files/Graphics/model_data/!test_data/coverage_matrix.txt", std::ios::out);
	off << m << " " << n << std::endl;
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			off << is_coveraged[i][j];
			if (j == n - 1) off << std::endl;
			else off << " ";
		}
	}
	off.close();*/
	time_t beg_gurobi = clock();
	std::cout << "\n[Gurobi] begin...\n\n";

	try {

		// Create an environment
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", "mip1_gurobi.log");
		env.set(GRB_DoubleParam::GRB_DoubleParam_TimeLimit, 200);// 设置最长优化时间为 300秒
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);

		// Create variables
		vector<GRBVar> vec_decision;
		vec_decision.reserve(n);
		for (int i = 0; i < n; ++i) {
			GRBVar var = model.addVar(0.0, 1.0, NULL, GRB_BINARY);// 决策变量 0-1
			vec_decision.emplace_back(var);
		}

		GRBVar Z = model.addVar(1, n, 0, GRB_INTEGER, "Z");

		// Set objective: 最小化中轴球的个数
		GRBLinExpr obj = Z;
		model.setObjective(obj, GRB_MINIMIZE);

		for (int i = 0; i < m; ++i) {// 每个采样点都需要被覆盖
			GRBLinExpr lhs = 0;
			for (int j = 0; j < n; ++j) {
				lhs += is_coveraged[i][j] * vec_decision[j];
			}
			model.addConstr(lhs >= 1);
		}

		GRBLinExpr inner_point_cnt = 0;// 中轴点的个数应尽可能少
		for (int i = 0; i < n; ++i) {
			inner_point_cnt += vec_decision[i];
		}
		model.addConstr(Z == inner_point_cnt);

		// Optimize model
		model.optimize();

		cout << "[Gurobi] Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

		for (int i = 0; i < n; ++i) {
			auto& var = vec_decision[i];
			if (var.get(GRB_DoubleAttr_X) > 0) {
				inner_points_selected.insert(inner_points[i]);
				vhs_selected.insert(ip_to_vh[inner_points[i]]);
			}
		}

		std::cout << "[Coverage Axis]: We selecte " << inner_points_selected.size() << " points.\n";

	} catch (GRBException e) {
		cout << "[Gurobi] Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "[Gurobi] Exception during optimization" << endl;
	}
	time_t end_gurobi = clock();
	std::cout << "\n[Gurobi] end... cost " << int(end_gurobi - beg_gurobi) << " ms. \n\n";

}

void CoverageAxis::connection_establishment() {

	// 1. 收缩未被选择的边, 只要这条边的两个顶点不是全被选择的, 那就要被收缩
	while (1) {

		unordered_set<EH> ehs_candidated;
		for (auto& ep : medial_mesh.alledges()) {
			//if (vhs_selected.count(ep.second.vh1()) && vhs_selected.count(ep.second.vh2())) continue;
			int degree = 0;
			if (vhs_selected.count(ep.second.vh1())) degree++;
			if (vhs_selected.count(ep.second.vh2())) degree++;
			if (degree == 1) {
				ehs_candidated.insert(ep.first);
			}
		}

		if (ehs_candidated.empty()) {
			std::cout << "Edges collapse completed.\n";
			break;
		}

		for (auto& eh : ehs_candidated) {

			if (medial_mesh.isValid(eh)) {
				auto& edge = medial_mesh.edges(eh);
				VH vh_keep = vhs_selected.count(edge.vh1()) ? edge.vh1() : edge.vh2();
				VH vh_remove(edge.vh1() + edge.vh2() - vh_keep);
				for (auto& adjvh : medial_mesh.NeighborVh(vh_remove)) {
					if (adjvh == vh_keep) continue;
					medial_mesh.AddEdge(adjvh, vh_keep);
				}
				/*for (auto& adjeh : medial_mesh.NeighborEh(vh_remove)) {
					medial_mesh.DeleteEdge(adjeh);
				}*/
				medial_mesh.DeleteVertex(vh_remove);
			}

		}


	}

	medial_mesh.updateAllHandles();

	// 2. 填补空洞
	hole_filling();


	std::cout << "[Coverage Axis]: The medial mesh has " << medial_mesh.vsize() << " vertices, "
		<< medial_mesh.esize() << " edges and "
		<< medial_mesh.fsize() << " faces after collaspsing.\n";


}

void CoverageAxis::connection_establishment(vector<bool>& is_selected) {

	// 方式二: 有面信息, 需要维护面拓扑
	std::queue<EH> ehs_need_collapse;
	for (auto& ep : medial_mesh.alledges()) {
		if (!is_selected[ep.second.vh1()] || !is_selected[ep.second.vh2()]) {
			ehs_need_collapse.push(ep.first);
		}
	}
	std::cout << "[Coverage Axis]: mark non-selected edges success.\n";

	while (!ehs_need_collapse.empty()) {

		auto eh = ehs_need_collapse.front();
		ehs_need_collapse.pop();
		if (!medial_mesh.isValid(eh)) continue;
		auto& edge = medial_mesh.edges(eh);
		VH vh_keep = edge.vh1();
		VH vh_remove = edge.vh2();
		if (!is_selected[vh_keep]) {
			swap(vh_keep, vh_remove);
		}
		std::vector<std::vector<VH>> new_triangles;
		for (auto& adjeh : medial_mesh.NeighborEh(vh_remove)) {
			auto& adje = medial_mesh.edges(adjeh);
			VH vh_other(adje.vh1() + adje.vh2() - vh_remove);
			if (vh_other == vh_keep || medial_mesh.isConnected(vh_other, vh_keep)) continue;
			EH eh_new = medial_mesh.AddEdge(vh_keep, vh_other);
			if (!is_selected[vh_keep] || !is_selected[vh_other]) {
				ehs_need_collapse.push(eh_new);
			}
		}
		for (auto& fh : medial_mesh.NeighborFh(vh_remove)) {
			auto vex = (medial_mesh.faces(fh)).getVertexHandle();
			if (std::find(vex.begin(), vex.end(), vh_keep) != vex.end()) continue;
			for (int i = 0; i < vex.size(); ++i) {
				if (vex[i] == vh_remove) {
					vex[i] = vh_keep;
					new_triangles.push_back(vex);
					break;
				}
			}
		}
		medial_mesh.DeleteVertex(vh_remove);

		// 添加三角形
		for (auto& tri : new_triangles) {
			medial_mesh.AddFace(tri);
		}

	}
	std::cout << "[Coverage Axis]: collapse non-selected edges success.\n";

	medial_mesh.updateAllHandles();

}

void CoverageAxis::get_inner_points_selected(vector<Vector3d>& inner_points_pos) {

	if (!is_computed) {
		compute_coverage_axis();
	}

	inner_points_pos.clear();
	inner_points_pos.reserve(inner_points_selected.size());

	for (auto* inner_point : inner_points_selected) {

		Vector3d pos = inner_point->pos;
		inner_points_pos.emplace_back(pos);

	}

}

void CoverageAxis::get_medial_mesh(SkeletalMesh& _medial_mesh) {

	if (!is_computed) {
		compute_coverage_axis();
	}

	_medial_mesh = medial_mesh;

}

bool CoverageAxis::find_boundary_edges_loop(EH eh, vector<VH>& vhs, vector<EH>& ehs) {
	if (vhs.size() > 4) return false;// 暂时只考虑四边及以内的空洞
	// 应该采取回溯的思想
	VH vh_curr = vhs.back();
	EH eh_next(-1);
	unordered_set<EH> eh_in_same_face;// 在同一个面上的边, 不应该被考虑
	for (auto& fh : medial_mesh.NeighborFh(eh)) {
		auto& face = medial_mesh.faces(fh);
		const auto& ehs = face.getEdgeHandle();
		for (auto& eh_same : ehs) {
			eh_in_same_face.insert(eh_same);
		}
	}
	for (auto& adjeh : medial_mesh.NeighborEh(eh)) {
		if (eh_in_same_face.count(adjeh)) continue;
		if (medial_mesh.NeighborFh(adjeh).size() < 2) {
			auto& adje = medial_mesh.edges(adjeh);
			if (adje.vh1() == vh_curr || adje.vh2() == vh_curr) {// 与这条边邻接, 且不是上一条边
				eh_next = adjeh;
				auto& edge_next = medial_mesh.edges(eh_next);
				VH vh_next(edge_next.vh1() + edge_next.vh2() - vh_curr);
				if (vhs.front() == vh_next) return true;// 已形成循环, 可以开始添加面了
				vhs.push_back(vh_next);
				ehs.push_back(eh_next);
				if (find_boundary_edges_loop(eh_next, vhs, ehs)) {
					return true;
				} else {
					// 回溯
					vhs.pop_back();
					ehs.pop_back();
				}
			}
		}
	}

	return false;

}

void CoverageAxis::hole_filling() {

	for (auto& ep : medial_mesh.alledges()) {// 填补三边形空洞
		if (medial_mesh.NeighborFh(ep.first).size() >= 2) continue;
		auto& edge = medial_mesh.edges(ep.first);
		VH vh1 = edge.vh1(), vh2 = edge.vh2();
		const auto& adjvhs1 = medial_mesh.NeighborVh(vh1);
		const auto& adjvhs2 = medial_mesh.NeighborVh(vh2);
		unordered_set<VH> vhs_common;
		for (auto& adjvh : adjvhs1) {
			if (adjvhs2.count(adjvh)) {
				vhs_common.insert(adjvh);
			}
		}
		for (auto& adjvh : vhs_common) {
			medial_mesh.AddFace({ vh1, vh2, adjvh });
		}
	}

	for (auto& ep : medial_mesh.alledges()) {// 填补不在边界且不与空洞邻接的空洞
		if (medial_mesh.NeighborFh(ep.first).size() >= 2) continue;
		vector<VH> vhs = { ep.second.vh1(), ep.second.vh2() };
		vector<EH> ehs = { ep.first };
		bool is_bdy_loop = find_boundary_edges_loop(ep.first, vhs, ehs);// 目前只考虑四边空洞
		if (is_bdy_loop) {
			for (int i = 2; i < vhs.size(); ++i) {// 将这圈循环的空洞三角化
				vector<VH> vhs_tri = { vhs[0], vhs[i - 1], vhs[i] };
				medial_mesh.AddFace(vhs_tri);
			}
		}
	}

}

bool CoverageAxis::neighbor_with_hanging_edge(VH vh) {

	for (auto& eh : medial_mesh.NeighborEh(vh)) {
		if (medial_mesh.NeighborFh(eh).size() < 2) {
			return true;
		}
	}
	return false;
}

bool CoverageAxis::is_all_hanging(vector<EH> ehs) {

	for (auto& eh : ehs) {
		if (medial_mesh.NeighborFh(eh).size() >= 2) {
			return false;
		}
	}
	return true;
}

void CoverageAxis::postprocessing() {

	bool topo_change = false;
	Vector3d dir(1, 1, 1);
	auto edges_origin = medial_mesh.alledges();
	for (auto& ep : edges_origin) {// 删除暴露在外面的边
		if (medial_mesh.NeighborFh(ep.first).size() > 1) continue;// 非边界
		auto& v1 = medial_mesh.vertices(ep.second.vh1());
		auto& v2 = medial_mesh.vertices(ep.second.vh2());
		bool outer_surface = false;

		for (double r = 0.1; r < 0.99; r += 0.1) {// 在边上采样
			Vex v = v1 * r + v2 * (1.0 - r);
			Ray ray(Vector3d(v.x(), v.y(), v.z()), dir);
			auto isec = bvh_tree->getIntersection(ray);
			if (isec.intersection_count % 2 == 0) {
				outer_surface = true;
				break;
			}
		}

		if (outer_surface) {
			medial_mesh.DeleteEdge(ep.first);
			topo_change = true;
		}
	}

	/*auto faces_origin = medial_mesh.allfaces();
	vector<vector<double>> weight_triangle = { { 0.2, 0.2, 0.6 }, { 0.2, 0.6, 0.2 }, { 0.6, 0.2, 0.2 }, { 0.4, 0.4, 0.2 },{ 0.4, 0.2, 0.4 }, { 0.2, 0.4, 0.4 }, { 0.333, 0.333, 0.333 } };
	for (auto& fp : faces_origin) {
		bool outer_surface = false;
		const auto& vhs = fp.second.getVertexHandle();
		vector<Vex> vs;
		for (auto& vh : vhs) {
			vs.emplace_back(medial_mesh.vertices(vh));
		}
		for (auto& wt : weight_triangle) {
			Vex v = vs[0] * wt[0] + vs[1] * wt[1] + vs[2] * wt[2];
			Ray ray(Vector3d(v.x(), v.y(), v.z()), dir);
			auto isec = bvh_tree->getIntersection(ray);
			if (isec.intersection_count % 2 == 0) {
				outer_surface = true;
				break;
			}
		}
		if (outer_surface) {
			medial_mesh.DeleteFace(fp.first);
			topo_change = true;
		}
	}*/

	if (topo_change) medial_mesh.updateAllHandles();

}


CoverageAxis_Vorn::CoverageAxis_Vorn(MeshKernel::SurfaceMesh& _surface_mesh, SkeletalMesh& _medial_mesh)
	: mesh(_surface_mesh), medial_mesh(_medial_mesh) {

	// 1. 初始化 AABB & BVH 树
	std::vector<Vector3f> vertices;
	for (auto& fp : mesh.allfaces()) {
		auto vhs = fp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& v = mesh.vertices(vh);
			vertices.push_back(Vector3f(v.x(), v.y(), v.z()));
		}
	}
	ab_tree = new AABB_Tree(vertices);

	bvh_tree = new BVH_Tree();
	bvh_tree->buildBVH_Tree(mesh);

	// 2. 初始化半径偏移
	mesh.initBBox();
	double diagonal_length = (mesh.BBoxMax - mesh.BBoxMin).norm();
	offset_radius *= diagonal_length;

	// 随机采样
	int sample_points_count = std::min(1500, int(mesh.fsize() * 1.25));
	map<double, FH> area2fh;
	double area_sum = 0.0;
	for (auto& fp : mesh.allfaces()) {
		const auto& vhs = fp.second.getVertexHandle();
		auto& v0 = mesh.vertices(vhs[0]);
		auto& v1 = mesh.vertices(vhs[1]);
		auto& v2 = mesh.vertices(vhs[2]);
		Vec vec01 = v1 - v0;
		Vec vec02 = v2 - v0;
		double area = (vec01.cross(vec02)).norm();
		area_sum += area;
		area2fh[area_sum] = fp.first;
	}

	while (surface_points.size() < sample_points_count) {
		double area = (rand() % RAND_MAX) * 1.0 / RAND_MAX * area_sum;
		auto it = area2fh.lower_bound(area);
		FH fh = (*it).second;
		if (!mesh.isValid(fh)) continue;
		double wa = 1.0 * (rand() % 10000007);
		double wb = 1.0 * (rand() % 10000007);
		double wc = 1.0 * (rand() % 10000007);
		double ws = wa + wb + wc;
		wa /= ws;
		wb /= ws;
		wc /= ws;
		const auto& vhs = mesh.faces(fh).getVertexHandle();
		Vex point = mesh.vertices(vhs[0]) * wa + mesh.vertices(vhs[1]) * wb + mesh.vertices(vhs[2]) * wc;
		SurfacePoint* sp = new SurfacePoint();
		sp->pos = Vector3d(point.x(), point.y(), point.z());
		sp->id = surface_points.size();
		surface_points.emplace_back(sp);
	}

	std::cout << "We get " << surface_points.size() << " surface points.\n";

}

void CoverageAxis_Vorn::compute_coverage_axis() {

	inner_points_generation();

	point_selection_based_on_set_coverage_MILP();
	medial_mesh.writeMaFile("C:/My Files/Graphics/model_data/!test_data/!CoverageAxis_NoPostprocessing.ca.ma");
	postprocessing();
	medial_mesh.writeMaFile("C:/My Files/Graphics/model_data/!test_data/!CoverageAxis_Final.ca.ma");
}

void CoverageAxis_Vorn::inner_points_generation() {

	int vcnt = medial_mesh.vsize();
	vector<bool> is_selected(vcnt, true);
	Vector3d dir(1, 1, 1);
#pragma omp parallel for
	for (int i = 0; i < vcnt; ++i) {// 并行版本
		VH vh(i);
		auto& v = medial_mesh.vertices(vh);
		Vector3d pos(v.x(), v.y(), v.z());
		Ray ray(pos, dir);
		auto intersection = bvh_tree->getIntersection(ray);
		if (intersection.intersection_count % 2 == 0) {
			is_selected[i] = false;
		}
	}

	connection_establishment(is_selected);

	for (int vid = 0; vid < medial_mesh.vsize(); vid++) {
		VH vh(vid);
		auto& v = medial_mesh.vertices(vh);
		InnerPoint* point = new InnerPoint();
		point->pos = Vector3d(v.x(), v.y(), v.z());
		point->id = vid;
		inner_points.emplace_back(point);
	}

	// 初始化内部点的半径
	int inner_pcnt = inner_points.size();
#pragma omp parallel for
	for (int i = 0; i < inner_pcnt; ++i) {
		auto* inner_point = inner_points[i];
		Vector3f pos(inner_point->pos.x(), inner_point->pos.y(), inner_point->pos.z()), np;
		double radius = ab_tree->findNearstPoint(pos, np);
		inner_point->radius = radius + offset_radius;// 膨胀后的半径
		medial_mesh.radius[VH(i)] = radius;// 未膨胀的半径
	}

	std::cout << "[Coverage Axis]: We get " << medial_mesh.vsize() << " inner points.\n";

}

void CoverageAxis_Vorn::point_selection_based_on_set_coverage_MILP() {

	int m = surface_points.size(), n = inner_points.size();
	vector<bool> is_selected(n, false);
	vector<vector<int>> is_coveraged(m, vector<int>(n, 0));// 0 表示未被覆盖, 1 表示被覆盖

#pragma omp parallel for
	for (int i = 0; i < inner_points.size(); ++i) {// 并行执行
		auto* inner_point = inner_points[i];
		for (auto* surface_point : surface_points) {
			double dist = (surface_point->pos - inner_point->pos).Length();
			if (dist <= inner_point->radius) {
				is_coveraged[surface_point->id][inner_point->id] = 1;
			}
		}
	}

	for (int i = 0; i < m; ++i) {
		int j = 0;
		while (j < n) {
			if (is_coveraged[i][j] == 1) break;
			j++;
		}
		if (j == n) {
			std::cerr << "Error: Exist surface points are not coveraged! So model is infeasible.\n";
			return;
		}
	}

	/*std::ofstream off("C:/My Files/Graphics/model_data/!test_data/coverage_matrix.txt", std::ios::out);
	off << m << " " << n << std::endl;
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			off << is_coveraged[i][j];
			if (j == n - 1) off << std::endl;
			else off << " ";
		}
	}
	off.close();*/
	time_t beg_gurobi = clock();
	std::cout << "\n[Gurobi] begin...\n\n";

	try {

		// Create an environment
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", "mip1_gurobi.log");
		env.set(GRB_DoubleParam::GRB_DoubleParam_TimeLimit, 100);// 设置最长优化时间为 300秒
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);

		// Create variables
		vector<GRBVar> vec_decision;
		vec_decision.reserve(n);
		for (int i = 0; i < n; ++i) {
			GRBVar var = model.addVar(0.0, 1.0, NULL, GRB_BINARY);// 决策变量 0-1
			vec_decision.emplace_back(var);
		}

		GRBVar Z = model.addVar(1, n, 0, GRB_INTEGER, "Z");

		// Set objective: 最小化中轴球的个数
		GRBLinExpr obj = Z;
		model.setObjective(obj, GRB_MINIMIZE);

		for (int i = 0; i < m; ++i) {// 每个采样点都需要被覆盖
			GRBLinExpr lhs = 0;
			for (int j = 0; j < n; ++j) {
				lhs += is_coveraged[i][j] * vec_decision[j];
			}
			model.addConstr(lhs >= 1);
		}

		GRBLinExpr inner_point_cnt = 0;// 中轴点的个数应尽可能少
		for (int i = 0; i < n; ++i) {
			inner_point_cnt += vec_decision[i];
		}
		model.addConstr(Z == inner_point_cnt);

		// Optimize model
		model.optimize();

		cout << "[Gurobi] Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

		for (int i = 0; i < n; ++i) {
			auto& var = vec_decision[i];
			if (var.get(GRB_DoubleAttr_X) > 0) {
				is_selected[i] = true;
			}
		}

		connection_establishment(is_selected);

		std::cout << "[Coverage Axis]: We selecte " << medial_mesh.vsize() << " points.\n";

	} catch (GRBException e) {
		cout << "[Gurobi] Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "[Gurobi] Exception during optimization" << endl;
	}
	time_t end_gurobi = clock();
	std::cout << "\n[Gurobi] end... cost " << int(end_gurobi - beg_gurobi) << " ms. \n\n";

}

void CoverageAxis_Vorn::connection_establishment(vector<bool>& is_selected) {

	// 方式二: 有面信息, 需要维护面拓扑
	std::queue<EH> ehs_need_collapse;
	for (auto& ep : medial_mesh.alledges()) {
		if (!is_selected[ep.second.vh1()] || !is_selected[ep.second.vh2()]) {
			ehs_need_collapse.push(ep.first);
		}
	}
	std::cout << "[Coverage Axis]: mark non-selected edges success.\n";

	while (!ehs_need_collapse.empty()) {

		auto eh = ehs_need_collapse.front();
		ehs_need_collapse.pop();
		if (!medial_mesh.isValid(eh)) continue;
		auto& edge = medial_mesh.edges(eh);
		VH vh_keep = edge.vh1();
		VH vh_remove = edge.vh2();
		if (!is_selected[vh_keep]) {
			swap(vh_keep, vh_remove);
		}
		std::vector<std::vector<VH>> new_triangles;
		for (auto& adjeh : medial_mesh.NeighborEh(vh_remove)) {
			auto& adje = medial_mesh.edges(adjeh);
			VH vh_other(adje.vh1() + adje.vh2() - vh_remove);
			if (vh_other == vh_keep || medial_mesh.isConnected(vh_other, vh_keep)) continue;
			EH eh_new = medial_mesh.AddEdge(vh_keep, vh_other);
			if (!is_selected[vh_keep] || !is_selected[vh_other]) {
				ehs_need_collapse.push(eh_new);
			}
		}
		for (auto& fh : medial_mesh.NeighborFh(vh_remove)) {
			auto vex = (medial_mesh.faces(fh)).getVertexHandle();
			if (std::find(vex.begin(), vex.end(), vh_keep) != vex.end()) continue;
			for (int i = 0; i < vex.size(); ++i) {
				if (vex[i] == vh_remove) {
					vex[i] = vh_keep;
					new_triangles.push_back(vex);
					break;
				}
			}
		}
		medial_mesh.DeleteVertex(vh_remove);

		// 添加三角形
		for (auto& tri : new_triangles) {
			medial_mesh.AddFace(tri);
		}

	}
	std::cout << "[Coverage Axis]: collapse non-selected edges success.\n";

	medial_mesh.updateAllHandles();

}

void CoverageAxis_Vorn::postprocessing() {

	bool topo_change = false;
	auto edges_origin = medial_mesh.alledges();
	Vector3d dir(1, 1, 1);
	for (auto& ep : edges_origin) {// 删除暴露在外面的边
		//if (medial_mesh.NeighborFh(ep.first).size() > 1) continue;// 非边界
		auto& v1 = medial_mesh.vertices(ep.second.vh1());
		auto& v2 = medial_mesh.vertices(ep.second.vh2());
		bool outer_surface = false;

		for (double r = 0.1; r < 0.99; r += 0.1) {// 在边上采样
			Vex v = v1 * r + v2 * (1.0 - r);
			Ray ray(Vector3d(v.x(), v.y(), v.z()), dir);
			auto isec = bvh_tree->getIntersection(ray);
			if (isec.intersection_count % 2 == 0) {
				outer_surface = true;
				break;
			}
		}

		if (outer_surface) {
			medial_mesh.DeleteEdge(ep.first);
			topo_change = true;
		}
	}

	if (topo_change) medial_mesh.updateAllHandles();

}
