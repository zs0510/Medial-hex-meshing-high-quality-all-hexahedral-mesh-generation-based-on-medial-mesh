#include "BoundarySmoother.h"

void BoundarySmoother::smoothing_2D(MeshKernel::SurfaceMesh& mesh, int iters, double move_speed) {

	std::unordered_map<VH, Vex> bd_vhs;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			bd_vhs[vp.first] = vp.second;
		}
	}


	std::cout << "Iter Max = " << iters << ", move speed = " << move_speed << std::endl;
	std::cout << "Boundary vertices size = " << bd_vhs.size() << std::endl;
	double cos_PI_4 = 0.7071;

	for (int it = 0; it < iters; ++it) {
		std::cout << "Iter " << it << " is running..."<< std::endl;
		std::unordered_map<VH, Vex> ref_bd_vhs = bd_vhs;

		bool update_flag = false;
		for (auto& vp : ref_bd_vhs) {
			std::vector<Vec> bd_vec;
			std::vector<Vex> bd_vex;
			auto& v = vp.second;
			const auto& adjehs = mesh.NeighborEh(vp.first);
			for (auto& adjeh : adjehs) {
				if (mesh.isOnBoundary(adjeh)) {
					auto& adje = mesh.edges(adjeh);
					VH adjvh(adje.vh1() + adje.vh2() - vp.first);
					auto& adjv = ref_bd_vhs[adjvh];
					bd_vec.push_back((adjv - v).normalized());
					bd_vex.push_back(adjv);
				}
			}
			
			if (bd_vec.size() != 2) {
				std::cerr << "[Error]: Exist vertex has less than 2 neighbor vertices. vsize = " << bd_vec.size() << std::endl;
				continue;
			}
			

			/*Vec np;
			double dis;
			distance_to_line(v, bd_vex[0], bd_vex[1], dis, np);
			update_flag = true;
			bd_vhs[vp.first] = v * (1 - move_speed) + np * move_speed;*/

			/*double cos_theta = bd_vec[0].dot(bd_vec[1]);
			if (cos_theta > 0.25) {
				
			} else {
				std::cout << "Need not to smooth " << vp.first << std::endl;
			}*/
			Vec np;
			double dis;
			distance_to_line(v, bd_vex[0], bd_vex[1], dis, np);
			update_flag = true;
			bd_vhs[vp.first] = v * (1 - move_speed) + np * move_speed;
			
		}

		if (!update_flag) break;

	}

	for (auto& vp : bd_vhs) {
		auto& v = mesh.vertices(vp.first);
		v = vp.second;
	}

}

void BoundarySmoother::smoothing_2D(SkeletalMesh& sklmesh, std::vector<int>& non_mainfold_vhs, int iters, double move_speed) {

	std::unordered_map<VH, Vex> vertices_bdy;

	std::unordered_map<EH, bool> curve_skel_ehs;
	std::unordered_map<EH, bool> face_skel_ehs;

	for (auto& ep : sklmesh.alledges()) {
		auto adjfhs = sklmesh.NeighborFh(ep.first);
		if (adjfhs.empty()) {
			curve_skel_ehs[ep.first] = true;// 曲线边的度为0
		} else if (adjfhs.size() == 1) {
			face_skel_ehs[ep.first] = true;
		}
	}

	for (auto& vp : sklmesh.allvertices()) {

		auto& ehs = sklmesh.NeighborEh(vp.first);
		//int hanging_count = 0;// 记录悬挂边的数目
		int hanging_ecnt = 0, boundary_ecnt = 0;
		for (auto& eh : ehs) {
			if (curve_skel_ehs.count(eh)) {
				hanging_ecnt++;
			} else if (face_skel_ehs.count(eh)) {
				boundary_ecnt++;
			}
		}
		
		if (boundary_ecnt > 1) {
			vertices_bdy[vp.first] = vp.second;
		}

	}

	std::unordered_map<VH, std::vector<VH>> adj_bdy_vhs;

	for (auto& vp : vertices_bdy) {
		/*const auto& adjvhs = sklmesh.NeighborVh(vp.first);
		for (auto& adjvh : adjvhs) {
			if (vertices_bdy.count(adjvh)) {
				adj_bdy_vhs[vp.first].emplace_back(adjvh);
			}
		}*/
		const auto& adjehs = sklmesh.NeighborEh(vp.first);
		for (auto& adjeh : adjehs) {
			if (face_skel_ehs.count(adjeh)) {
				auto& e = sklmesh.edges(adjeh);
				VH adjvh(e.vh1() + e.vh2() - vp.first);
				adj_bdy_vhs[vp.first].emplace_back(adjvh);
			}
		}

		if (adj_bdy_vhs[vp.first].size() != 2) {
			std::cerr << "[Error]: Exist vertex has less than 2 neighbor vertices. vsize = " 
				<< adj_bdy_vhs[vp.first].size() << std::endl;
			non_mainfold_vhs.emplace_back(vp.first);
		}

	}

	if (!non_mainfold_vhs.empty()) return;

	auto tmp_sklmesh = sklmesh;
	Project_KD_Tree_SkeletalMesh kdtree(tmp_sklmesh);

	for (int it = 0; it < iters; ++it) {
		std::unordered_map<VH, Vex> ref_bd_vhs = vertices_bdy;
		for (auto& vp : ref_bd_vhs) {
			VH vh = vp.first;
			auto& v = vp.second;
			auto& adjv0 = ref_bd_vhs[adj_bdy_vhs[vh][0]];
			auto& adjv1 = ref_bd_vhs[adj_bdy_vhs[vh][1]];
			Vec np;
			double dis;
			distance_to_line(v, adjv0, adjv1, dis, np);
			vertices_bdy[vp.first] = v * (1 - move_speed) + np * move_speed;
			vertices_bdy[vp.first] = kdtree.get_nearest(vertices_bdy[vp.first]);
		}

	}

	
	for (auto& vp : vertices_bdy) {
		auto& v = sklmesh.vertices(vp.first);
		v = vp.second;
		sklmesh.radius[vp.first] = kdtree.get_radius(v);// 更新半径
	}

}

bool BoundarySmoother::distance_to_line(const Vex& p, const Vex& v0, const Vex& v1, double& dist, Vex& np) {
	/*auto& p = mesh.vertices(vh);
	auto& v0 = mesh.vertices(vh0);
	auto& v1 = mesh.vertices(vh1);*/
	Vec v0v1(v1 - v0), pv0(v0 - p), pv1(v1 - p);
	double area = fabs(v0v1.cross(pv0).norm());
	if (v0v1.norm() > 1e-12) {
		dist = area / v0v1.norm();
		double t = (pv0.dot(pv0) - pv0.dot(pv1)) / (pv0.dot(pv0) + pv1.dot(pv1) - 2 * pv0.dot(pv1));
		np = v0 * (1 - t) + v1 * t;
		return true;
	}
	return false;
}

void BoundarySmoother::smoothing_branch(SkeletalMesh& sklmesh) {

	unordered_set<EH> curve_ehs;
	for (auto& ep : sklmesh.alledges()) {
		if (sklmesh.NeighborFh(ep.first).empty()) {
			curve_ehs.insert(ep.first);
		}
	}
	sklmesh.genAllFacesNormal();
	int branch_cnt = 0, smooth_cnt = 0;
	for (auto& vp : sklmesh.allvertices()) {
		auto& v = sklmesh.vertices(vp.first);
		auto& ehs = sklmesh.NeighborEh(vp.first);
		int hanging_count = 0;// 记录悬挂边的数目
		Vex dir(0, 0, 0);
		for (auto& eh : ehs) {
			if (curve_ehs.count(eh)) {
				auto& e = sklmesh.edges(eh);
				VH vh_t(e.vh1() + e.vh2() - vp.first);
				auto& v_t = sklmesh.vertices(vh_t);
				dir += (v_t - v).normalized();
				hanging_count++;
			}
		}
		if (hanging_count == 0) continue;// 面节点
		dir.normalize();

		const auto& adjfhs = sklmesh.NeighborFh(vp.first);
		if (adjfhs.empty()) continue;
		branch_cnt++;
		bool is_parallel_with_normal = false;// 与法向量是否足够平行
		for (auto& fh : adjfhs) {
			auto& face = sklmesh.faces(fh);
			Vex normal = face.getNormal();
			if (abs(normal.norm2() - 1.0) > 1e-6f) {
				std::cerr << "Error: Normal is not unit!!!" << std::endl;
				continue;
			}
			if (abs(dir.dot(normal)) > 0.95) {
				is_parallel_with_normal = true;
				break;
			}
		}
		if (is_parallel_with_normal) continue;// 与法向量接近平行的不需要做光滑

		
		const auto& adjvhs = sklmesh.NeighborVh(vp.first);
		bool is_perp_with_neighbor = false;// 与邻域面边界边是否足够垂直
		vector<Vex> neighbors_bdy;
		for (auto& adjvh : adjvhs) {
			for (auto& adjeh : sklmesh.NeighborEh(adjvh)) {
				if (sklmesh.NeighborFh(adjeh).size() == 1) {
					auto& adjv = sklmesh.vertices(adjvh);
					neighbors_bdy.emplace_back(adjv);
					Vex vec = (adjv - v).normalized();
					if (abs(dir.dot(vec)) < 0.05) {
						is_perp_with_neighbor = true;
					}
					break;
				}
			}
		}
		if (is_perp_with_neighbor) continue;// 和面边界边太垂直的边不需要做光滑

		if (neighbors_bdy.size() == 2) {
			smooth_cnt++;
			Vex v_new = (neighbors_bdy.front() + neighbors_bdy.back()) * 0.5;
			v =  v_new * 0.95 + v * 0.05;
		}

	}

	std::cout << "We find " << branch_cnt << " branch vertices and smooth " << smooth_cnt << " vertices.\n";

}