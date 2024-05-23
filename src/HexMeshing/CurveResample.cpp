#include "CurveResample.h"

void CurveResample::resample() {

	
	init_nodes_map();

	std::unordered_map<VH, std::unordered_map<EH, bool>> visited_bvh;

	std::vector<std::pair<VH, VH>> new_edges;

	for (auto& vp : branch_vhs) {
		VH bvh(vp.first);
		const auto& adjehs = mesh.NeighborEh(bvh);
		for (auto& adjeh : adjehs) {
			if (visited_bvh[bvh].count(adjeh) || !curve_ehs.count(adjeh)) continue;
			
			EH eh_cur = adjeh;
			auto e_cur = mesh.edges(eh_cur);
			VH vh_cur(e_cur.vh1() + e_cur.vh2() - bvh);
			std::vector<VH> line_vhs = { bvh, vh_cur };
			std::vector<EH> line_ehs = { eh_cur };
			while (joint_vhs.count(vh_cur)) {

				EH eh_next(-1);
				const auto& adjehs_t = mesh.NeighborEh(vh_cur);
				for (auto& adjeh_t : adjehs_t) {
					if (!curve_ehs.count(adjeh_t) || adjeh_t == eh_cur) continue;
					eh_next = adjeh_t;
					break;
				}

				if (eh_next == -1) {
					std::cout << "[Curve Resample]: Invalid lint eh!!!\n";
					return;
				}
				eh_cur = eh_next;
				e_cur = mesh.edges(eh_cur);
				vh_cur = VH(e_cur.vh1() + e_cur.vh2() - vh_cur);

				line_vhs.emplace_back(vh_cur);
				line_ehs.emplace_back(eh_cur);

			}
			if (branch_vhs.count(vh_cur)) {
				visited_bvh[vh_cur][eh_cur] = true;
			}
			int line_vsz = line_vhs.size();
			std::vector<bool> used(line_vsz, false);
			used.front() = used.back() = true;
			recusive_resample(1, line_vsz, line_vhs, used);
			VH pre_used = line_vhs.front();
			for (int i = 1; i < line_vhs.size(); ++i) {
				if (used[i]) {
					new_edges.emplace_back(pre_used, line_vhs[i]);
					pre_used = line_vhs[i];
				}
			}
		}
	}
	
	for (auto& eh : curve_ehs) {
		mesh.DeleteEdge(eh.first);
	}

	for (auto& vp : new_edges) {
		mesh.AddEdge(vp.first, vp.second);
	}

	mesh.updateAllHandles();

	SkeltalMesh_Operation::branch_merge(mesh);

}

void CurveResample::curve_subdivision(EH eh, int iter) {

	int n_seg = std::pow(2, iter);
	int added_vcnt =  n_seg - 1;
	if (added_vcnt <= 0) return;
	auto& edge = mesh.edges(eh);
	auto vh1 = edge.vh1();
	auto vh2 = edge.vh2();
	auto& v1 = mesh.vertices(vh1);
	auto& v2 = mesh.vertices(vh2);

	std::vector<VH> added_vhs;
	for (int i = 1; i <= added_vcnt; ++i) {
		double weight = i * 1.0 / n_seg;
		Vex pos = v1 * (1 - weight) + v2 * weight;
		double radius = mesh.radius[vh1] * (1 - weight) + mesh.radius[vh2] * weight;
		VH newvh = mesh.AddVertex(pos);
		mesh.radius[newvh] = radius;
		added_vhs.push_back(newvh);
	}

	mesh.DeleteEdge(eh);
	mesh.AddEdge(vh1, added_vhs.front());
	for (int i = 1; i < added_vhs.size(); ++i) {
		mesh.AddEdge(added_vhs[i - 1], added_vhs[i]);
	}
	mesh.AddEdge(vh2, added_vhs.back());

}

void CurveResample::recusive_resample(int left, int right, std::vector<VH>& vhs, std::vector<bool>& used) {
	/*
	We start with a dense sampling of the curve-skeleton
	(left) and we re-sample it by iteratively splitting half-way each of
	its portions. The splitting process ends when the maximal sphere
	centered at the new sample intersects both the spheres centered at
	the two end-points of the current segment (middle). The resulting
	coarse sampling of the skeleton (right) determines the connectivity
	of the final hexahedral mesh.
	*/
	if (left > right) return;
	int mid = (left + right) / 2;
	if (!used[mid]) {
		int pre = mid - 1, rear = mid + 1;
		while (pre >= 0 && used[pre] == false) pre--;
		while (rear < vhs.size() && used[rear] == false) rear++;
		assert(pre >= 0 && rear < vhs.size());
		bool pre_intersect, rear_intersect;
		double pre_dis = (mesh.vertices(vhs[mid]) - mesh.vertices(vhs[pre])).norm();
		double rear_dis = (mesh.vertices(vhs[mid]) - mesh.vertices(vhs[rear])).norm();
		pre_intersect = (pre_dis < mesh.radius[vhs[mid]] + mesh.radius[vhs[pre]]);
		rear_intersect = (rear_dis < mesh.radius[vhs[mid]] + mesh.radius[vhs[rear]]);
		if (!pre_intersect && !rear_intersect) {
			used[mid] = true;
		}
	}
	recusive_resample(left, mid - 1, vhs, used);
	recusive_resample(mid + 1, right, vhs, used);
}

void CurveResample::init_nodes_map() {

	for (auto& ep : mesh.alledges()) {
		const auto& adjfhs = mesh.NeighborFh(ep.first);
		if (adjfhs.empty()) {
			curve_ehs[ep.first] = true;
		}
	}



	for (auto& vp : mesh.allvertices()) {

		VH vh(vp.first);

		const auto& adjfhs = mesh.NeighborFh(vh);
		const auto& adjehs = mesh.NeighborEh(vh);

		bool neighbor_with_curve = false;
		for (auto& adjeh : adjehs) {
			if (curve_ehs.count(adjeh)) {
				neighbor_with_curve = true;
				break;
			}
		}

		if (!neighbor_with_curve) continue;// 纯纯的三角面片顶点, 无需处理

		if (adjfhs.empty()) {

			if (adjehs.size() > 2) {
				branch_vhs[vh] = true;// 线线分支
			} else if (adjehs.size() == 2) {
				joint_vhs[vh] = true;// 线上结点
			} else if (adjehs.size() == 1) {
				end_vhs[vh] = true;// 线端结点
			}

		} else {
			branch_vhs[vh] = true;// 线面分支
		}


	}

	std::cout << "[Curve Resample]: branch nodes size = " << branch_vhs.size() << ", "
		<< "joint nodes size = " << joint_vhs.size() << ", "
		<< "end nodes size = " << end_vhs.size() << std::endl;

}

//{
//	auto oldehs = mesh.alledges();
//	for (auto& ep : oldehs) {
//		if (mesh.isValid(ep.first) && mesh.NeighborFh(ep.first).empty()) {
//			curve_subdivision(ep.first, 3);
//		}
//	}
//
//
//	std::unordered_map<VH, bool> branch_vhs;
//	std::unordered_map<VH, bool> joint_vhs;
//	std::unordered_map<VH, bool> end_vhs;
//
//	for (auto& vp : mesh.allvertices()) {
//
//		VH vh(vp.first);
//		const auto& adjfhs = mesh.NeighborFh(vh);
//		const auto& adjehs = mesh.NeighborEh(vh);
//		if (adjfhs.empty()) {
//			if (adjehs.size() == 2) {// 暂时只考虑出现在一条曲线上的关节节点
//				joint_vhs[vh] = true;
//			} else if (adjehs.size() == 1) {
//				end_vhs[vh] = true;
//			} else if (adjehs.size() > 2) {
//				branch_vhs[vh] = true;
//			}
//		} else {
//			bool is_on_curve = false;
//			for (auto& eh : adjehs) {
//				if (mesh.NeighborFh(eh).empty()) {
//					is_on_curve = true;
//					break;
//				}
//			}
//			if (is_on_curve) branch_vhs[vh] = true;
//		}
//
//	}
//	std::cout << "Branch vhs.size = " << branch_vhs.size() << std::endl;
//	std::cout << "Joint vhs.size = " << joint_vhs.size() << std::endl;
//	std::cout << "End vhs.size = " << end_vhs.size() << std::endl;
//
//	std::unordered_map<VH, bool> visited;
//	std::vector<std::pair<VH, VH>> new_edges;
//	std::unordered_map<EH, bool> need_to_erased;
//	std::unordered_map<VH, std::unordered_map<EH, bool>> visited_bvh;
//	for (auto& vp : branch_vhs) {// 从分支结点出发
//		VH bvh = vp.first;
//		const auto& adjehs = mesh.NeighborEh(bvh);
//		for (auto& eh : adjehs) {// 从分支结点的边出发
//			if (!mesh.isValid(eh) || visited_bvh[bvh].count(eh)) continue;// 无效或已被访问过
//			auto& e = mesh.edges(eh);
//			VH vh_cur(e.vh1() + e.vh2() - bvh);
//			if (!joint_vhs.count(vh_cur)) continue;// 若第一个邻接结点为分支结点或终端结点则不简化
//			std::vector<VH> curve_vhs;// 这条曲线上所有的点
//			std::vector<EH> curve_ehs;// 这条曲线上所有的边
//			curve_vhs.push_back(bvh);
//			std::unordered_map<EH, bool> visited_ehs;
//			EH eh_cur = eh;
//			while (1) {
//				curve_ehs.push_back(eh_cur);
//				curve_vhs.push_back(vh_cur);
//				if (branch_vhs.count(vh_cur)) {
//					visited_bvh[vh_cur][eh_cur] = true;
//					break;
//				}
//				EH eh_next(-1);
//				visited_ehs[eh_cur] = true;
//				for (auto& eh_n : mesh.NeighborEh(eh_cur)) {
//					if (visited_ehs.count(eh_n) || !mesh.NeighborFh(eh_n).empty()) continue;
//					auto& e_n = mesh.edges(eh_n);
//					VH vh_next(e_n.vh1() + e_n.vh2() - vh_cur);
//					if (!joint_vhs.count(vh_next) && !branch_vhs.count(vh_next) && !end_vhs[vh_next]) continue;
//					eh_next = eh_n;
//					vh_cur = vh_next;
//					break;
//				}
//				if (eh_next == -1) break;
//				eh_cur = eh_next;
//			}
//			std::cout << "vhs.size = " << curve_vhs.size() << ", ehs.size = " << curve_ehs.size() << std::endl;
//			visited[curve_vhs.back()] = true;
//			// 简化
//			for (auto& oldeh : curve_ehs) {
//				need_to_erased[oldeh] = true;
//			}
//			std::vector<bool> used(curve_vhs.size(), false);
//			used.front() = used.back() = true;// 首尾两点保证被使用
//			//used[curve_vhs.size() / 2.0] = true;
//			recusive_resample(1, curve_vhs.size() - 2, curve_vhs, used);
//			VH pre_used = curve_vhs.front();
//			for (int i = 1; i < curve_vhs.size(); ++i) {
//				if (used[i]) {
//					new_edges.emplace_back(pre_used, curve_vhs[i]);
//					pre_used = curve_vhs[i];
//				}
//			}
//		}
//	}
//
//	// 删除旧的边
//	for (auto& oldeh : need_to_erased) {
//		if (mesh.isValid(oldeh.first)) {
//			mesh.DeleteEdge(oldeh.first);
//		}
//	}
//
//	// 添加新的边
//	for (auto& vp : new_edges) {
//		mesh.AddEdge(vp.first, vp.second);
//	}
//
//	mesh.updateAllHandles();
//}