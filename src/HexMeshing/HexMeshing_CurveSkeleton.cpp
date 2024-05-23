#include "HexMeshing_CurveSkeleton.h"
#include "Volume_Algorithm/VolumeEditer.h"

using namespace std;

void HexMeshing_CurveSkeleton::hexmeshing() {

	init_nodes_map();
	std::cout << "Init nodes mapping success.\n";
	set_joint_plane();
	set_end_plane();
	std::cout << "Init joint and end nodes' planes success.\n";
	set_branch_face();
	std::cout << "Init branching nodes' faces success.\n";
	lanuch_from_branching_node();
	//std::cout << "Init branching nodes' faces success.\n";

	/*set_end_cube();
	set_joint_plane();*/

	cone_detection();

	

}

void HexMeshing_CurveSkeleton::init_nodes_map() {

	for (auto& ep : sklmesh.alledges()) {
		if (sklmesh.NeighborFh(ep.first).empty()) {
			curve_skel_ehs[ep.first] = true;// 曲线边的度为0
		}
	}

	for (auto& vp : sklmesh.allvertices()) {

		auto& ehs = sklmesh.NeighborEh(vp.first);
		int hanging_count = 0;// 记录悬挂边的数目
		for (auto& eh : ehs) {
			if (curve_skel_ehs.count(eh)) {
				hanging_count++;
			}
		}
		if (hanging_count == 0) continue;// 面节点
		int faces_count = sklmesh.NeighborFh(vp.first).size();

		if (faces_count != 0) {
			branching_nodes[vp.first] = true;
		} else {

			if (hanging_count == 1) {
				end_nodes[vp.first] = true;// 骨架线的端点
			} else if (hanging_count > 2 || is_corner(vp.first)) {
				branching_nodes[vp.first] = true;// 骨架线的分支节点
			} else if (hanging_count == 2) {
				joint_nodes[vp.first] = true;
			}

		}

	}

	std::cout << "Branching nodes size = " << branching_nodes.size() << std::endl;
	std::cout << "Joint nodes size = " << joint_nodes.size() << std::endl;
	std::cout << "End nodes size = " << end_nodes.size() << std::endl;

}

void HexMeshing_CurveSkeleton::set_branch_face() {

	BVH_Tree* bvh_tree = new BVH_Tree();
	bvh_tree->buildBVH_Tree(hexmesh);

	for (auto& branching_node : branching_nodes) {
		VH branching_vh = branching_node.first;
		Vex branching_vex = sklmesh.vertices(branching_vh);
		for (auto& eh : sklmesh.NeighborEh(branching_vh)) {
			if (curve_skel_ehs.count(eh)) {
				//std::cout << "Branching vh = " << branching_vh << ", eh = " << eh << "\n";
				auto& e = sklmesh.edges(eh);
				auto& adjv = sklmesh.vertices(VH(e.vh1() + e.vh2() - branching_vh));
				Vec dir = adjv - branching_vex;
				Ray ray(Vector3d(branching_vex.x(), branching_vex.y(), branching_vex.z())
					, Vector3d(dir.x(), dir.y(), dir.z()));
				auto intersection = bvh_tree->getIntersection(ray);
				if (intersection.happened == false) {
					Ray ray_invert(Vector3d(adjv.x(), adjv.y(), adjv.z())
						, -Vector3d(dir.x(), dir.y(), dir.z()));
					intersection = bvh_tree->getIntersection(ray_invert);
					if (intersection.happened == false) {
						std::cerr << "[HexMeshing Curve]: Error!!! Branching node no intersection with hexmesh.\n\n";
						continue;
					}
				}
				FH intersect_fh = FH(intersection.index);
				branchingnode_to_fh[branching_vh][eh] = intersect_fh;
				//// 将面的中心移至分支结点与面的交点
				//Vex center = hexmesh.getFaceCenter(intersect_fh);
				//auto& face = hexmesh.faces(intersect_fh);
				//Vec move = Vex(intersection.pos[0], intersection.pos[1], intersection.pos[2]) - center;
				//for (auto& hexvh : face.getVertexHandle()) {
				//	auto& hexv = hexmesh.vertices(hexvh);
				//	hexv += move;
				//}

			}
		}

	}

}

void HexMeshing_CurveSkeleton::lanuch_from_branching_node() {
	std::cerr << "We are launching from branching node. Current cells size = " << hexmesh.csize() << "\n";
	
	unordered_map<VH, unordered_map<EH, bool>> visited_b;

	for (auto& branching_node : branching_nodes) {
		VH branching_vh = branching_node.first;
		Vex branching_vex = sklmesh.vertices(branching_vh);

		for (auto& eh : sklmesh.NeighborEh(branching_vh)) {
			if (curve_skel_ehs.count(eh)) {
				if (visited_b[branching_vh].count(eh)) continue;

				//std::cout << "Branching vh = " << branching_vh << ", eh = " << eh << "\n";
				auto& e = sklmesh.edges(eh);
				auto adjvh = VH(e.vh1() + e.vh2() - branching_vh);
				
				std::vector<VH> vhs = { branching_vh, adjvh };

				if (dfs_curve_vhs(vhs)) {
					extrude_curve_vhs(vhs);
					visited_b[branching_vh][eh] = true;
					if (branching_nodes.count(vhs.back())) {
						visited_b[vhs.back()][sklmesh.getEdgeHandle(vhs[vhs.size() - 2], vhs.back())] = true;// 将这个分支结点的这个方向标记为已遍历
					}
				} else {
					std::cerr << "Invalid vhs.\n";
				}
			}
		}

	}
	std::cerr << "We are launched from branching node. Current cells size = " << hexmesh.csize() << "\n";

}

void HexMeshing_CurveSkeleton::extrude_curve_vhs(std::vector<VH>& vhs) {

	//std::cout << "Curve vhs size = " << vhs.size() << std::endl;
	int sz = vhs.size();

	if (sz == 2 && branching_nodes.count(vhs.back())) {// 两个分支结点直接相连
		EH skleh = sklmesh.getEdgeHandle(vhs.front(), vhs.back());
		FH fh0 = branchingnode_to_fh[vhs.front()][skleh];
		FH fh1 = branchingnode_to_fh[vhs.back()][skleh];
		HexAddCell::add_cell_2fh(hexmesh, fh0, fh1);
		return;
	}

	for (int i = 1; i < sz - 1; ++i) {
		VH prevh = vhs[i - 1], curvh = vhs[i];
		Vec dir = sklmesh.vertices(curvh) - sklmesh.vertices(prevh);
		vector<pair<VH, VH>> pair_vh;
		if (!branchingnode_to_fh.count(prevh) && !jointnode_to_fh.count(prevh)) continue;
		FH pre_fh = (branchingnode_to_fh.count(prevh)) ? 
			branchingnode_to_fh[prevh][sklmesh.getEdgeHandle(prevh, curvh)] : jointnode_to_fh[prevh];
		auto pre_fvhs = hexmesh.faces(pre_fh).getVertexHandle();
		for (auto& fvh : pre_fvhs) {
			Vex src = hexmesh.vertices(fvh);
			if (!jointnode_to_jp.count(curvh)) {
				std::cerr << "Invalid curvh!\n";
				return;
			}
			auto& plane = jointnode_to_jp[curvh];
			double t = -(plane.a * src.x() + plane.b * src.y() + plane.c * src.z() + plane.d) 
				/ ( plane.a * dir.x() + plane.b * dir.y() + plane.c * dir.z());
			Vex intersect_v = src + dir * t;
			Vec dir = (intersect_v - sklmesh.vertices(curvh)).normalized();
			Vex new_v = sklmesh.vertices(curvh) + dir * sklmesh.radius[curvh] * radius_scale;
			pair_vh.emplace_back( fvh, hexmesh.AddVertex(new_v));
		}
		jointnode_to_fh[curvh] = hexmesh.AddFace({ pair_vh[0].second, pair_vh[1].second, pair_vh[2].second, pair_vh[3].second });
		hexmesh.AddCell({ pair_vh[0].first, pair_vh[1].first, pair_vh[2].first, pair_vh[3].first, 
			pair_vh[0].second, pair_vh[1].second, pair_vh[2].second, pair_vh[3].second });
		
	}

	if (branching_nodes.count(vhs.back())) {
		VH b_vh = vhs.back(), j_vh = vhs[vhs.size() - 2];
		EH b_eh = sklmesh.getEdgeHandle(vhs.back(), vhs[vhs.size() - 2]);
		FH b_fh = branchingnode_to_fh[b_vh][b_eh];
		FH j_fh = jointnode_to_fh.count(j_vh) ? jointnode_to_fh[j_vh] : branchingnode_to_fh[j_vh][b_eh];
		if (!hexmesh.isValid(b_fh) || !hexmesh.isValid(j_fh)) {
			return;
		}
		VolumeEditer::add_cell(hexmesh, b_fh, j_fh);
		/*vector<pair<VH, VH>> pair_vh;
		unordered_set<VH> uniquevhs;
		Vec dir = (hexmesh.getFaceCenter(b_fh) - hexmesh.getFaceCenter(j_fh)).normalized();
		for (auto& bfvh : hexmesh.faces(b_fh).getVertexHandle()) {
			auto& bfv = hexmesh.vertices(bfvh);
			double cosine_max = -2;
			VH s_vh(-1);
			for (auto& jfvh : hexmesh.faces(j_fh).getVertexHandle()) {
				auto& jfv = hexmesh.vertices(jfvh);
				Vec vec = (bfv - jfv).normalized();
				double cosine = vec.dot(dir);
				if (cosine > cosine_max) {
					cosine_max = cosine;
					s_vh = jfvh;
				}
			}
			uniquevhs.insert(bfvh);
			uniquevhs.insert(s_vh);
			pair_vh.emplace_back(bfvh, s_vh);
		}
		if (uniquevhs.size() == 8) {
			hexmesh.AddCell({ pair_vh[0].first, pair_vh[1].first, pair_vh[2].first, pair_vh[3].first,
			pair_vh[0].second, pair_vh[1].second, pair_vh[2].second, pair_vh[3].second });
		}*/

	} else if (end_nodes.count(vhs.back())) {

		VH endvh = vhs.back(), prevh = vhs[sz - 2];
		auto& endv = sklmesh.vertices(endvh);
		auto& prev = sklmesh.vertices(prevh);
		FH prefh = jointnode_to_fh.count(prevh) ? jointnode_to_fh[prevh] 
			: branchingnode_to_fh[prevh][sklmesh.getEdgeHandle(endvh, prevh)];
		Vec dir = (endv - prev).normalized();
		vector<pair<VH, VH>> pair_vh;
		if (!hexmesh.isValid(prefh)) {
			std::cerr << "Line 218: Prefh is invalid!\n";
			return;
		}
		for (auto& hexvh : hexmesh.faces(prefh).getVertexHandle()) {
			auto& hexv = hexmesh.vertices(hexvh);
			auto& plane = endnode_to_jp[endvh];
			double t = -(plane.a * hexv.x() + plane.b * hexv.y() + plane.c * hexv.z() + plane.d)
				/ (plane.a * dir.x() + plane.b * dir.y() + plane.c * dir.z());
			Vex new_v = hexv + dir * t;
			pair_vh.emplace_back(hexvh, hexmesh.AddVertex(new_v));
		}
		hexmesh.AddCell({ pair_vh[0].first, pair_vh[1].first, pair_vh[2].first, pair_vh[3].first,
			pair_vh[0].second, pair_vh[1].second, pair_vh[2].second, pair_vh[3].second });

	}


}

bool HexMeshing_CurveSkeleton::dfs_curve_vhs(std::vector<VH>& vhs) {
	//std::cerr << "We are dfs curve vhs.\n";
	if (branching_nodes.count(vhs.back()) || end_nodes.count(vhs.back())) {
		return true;
	}
	VH pre = vhs[vhs.size() - 2], cur = vhs.back();// 当前结点只会是 joint node
	for (auto& eh : sklmesh.NeighborEh(cur)) {
		if (curve_skel_ehs.count(eh)) {
			auto& e = sklmesh.edges(eh);
			auto next = VH(e.vh1() + e.vh2() - cur);
			if (next == pre) continue;
			vhs.push_back(next);
			return dfs_curve_vhs(vhs);
		}
	}
	return false;
}

void HexMeshing_CurveSkeleton::set_joint_plane() {
	
	for (auto& joint_node : joint_nodes) {

		VH joint_vh = joint_node.first;// 骨架的曲线连接点
		Vex joint_vex = sklmesh.vertices(joint_vh);
		std::vector<Vec> dirs;
		for (auto& eh : sklmesh.NeighborEh(joint_vh)) {
			if (curve_skel_ehs.count(eh)) {
				auto& e = sklmesh.edges(eh);
				auto& adjv = sklmesh.vertices(VH(e.vh1() + e.vh2() - joint_vh));
				dirs.push_back((adjv - joint_vex).normalized());
				if (dirs.size() == 2) {
					if (dirs[1].dot(dirs[0]) < 0) {
						dirs[1] *= -1;
					}
					break;
				}
			}
		}
		if (dirs.size() < 2) {
			std::cerr << "[HexMeshing Curve]: Not a valid joint vertex.\n";
			continue;
		}
		Vec N = (dirs[0] + dirs[1]).normalized();
		double d = -( joint_vex.x() * N.x() + joint_vex.y() * N.y() + joint_vex.z() * N.z());
		jointnode_to_jp[joint_vh] = JointPlane(N.x(), N.y(), N.z(), d);

	}


}

void HexMeshing_CurveSkeleton::set_end_plane() {

	for (auto& end_node : end_nodes) {

		VH end_vh = end_node.first;// 骨架的曲线连接点
		Vex end_vex = sklmesh.vertices(end_vh);
		Vec dir;
		for (auto& eh : sklmesh.NeighborEh(end_vh)) {
			if (curve_skel_ehs.count(eh)) {
				auto& e = sklmesh.edges(eh);
				auto& adjv = sklmesh.vertices(VH(e.vh1() + e.vh2() - end_vh));
				dir = (end_vex - adjv).normalized();
				break;
			}
		}
		
		Vec N = dir.normalized();
		end_vex += N * sklmesh.radius[end_vh];// 往外移, 逼近原曲面
		double d = -(end_vex.x() * N.x() + end_vex.y() * N.y() + end_vex.z() * N.z());
		endnode_to_jp[end_vh] = JointPlane(N.x(), N.y(), N.z(), d);

	}

}

void HexMeshing_CurveSkeleton::cone_detection() {

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

bool HexMeshing_CurveSkeleton::is_corner(VH vh) {

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
			if (dirs[i].dot(dirs[j]) > 0) return true;
		}
	}

	return false;

}