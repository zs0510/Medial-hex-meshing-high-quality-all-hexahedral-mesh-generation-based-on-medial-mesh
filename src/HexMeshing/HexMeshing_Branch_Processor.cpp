#include "HexMeshing_Branch_Processor.h"

void HexMeshing_Branch_Processor::optimize(VH& algorithm_insert, bool only_place__not_split) {


	/*MeshKernel::IO io_input;
	io_input.iGameWriteVolumeFile(hexmesh, "C:/My Files/Graphics/model_data/!test_data/hexmesh_input.mesh");*/

	init_nodes_map();

	optimize_curve_branch();// 在线线分支结点上放置六面体

	if (hexmesh.csize() == 0) {
		// 没有分支节点, 需要选择一个分支节点放置六面体
		double most_orthogonal = 1.0;
		VH vh_selected(-1);
		vector<pair<VH, Vec>> vecs_selected;
		for (auto& vp : sklmesh.allvertices()) {
			VH sklvh = vp.first;
			auto& sklv = sklmesh.vertices(sklvh);
			auto& ehs = sklmesh.NeighborEh(sklvh);
			vector<pair<VH, Vec>> vecs_tmp;
			for (auto& eh : ehs) {
				if (curve_skel_ehs.count(eh)) {
					auto& edge = sklmesh.edges(eh);
					VH adjvh(edge.vh1() + edge.vh2() - sklvh);
					auto& adjv = sklmesh.vertices(adjvh);
					vecs_tmp.emplace_back(adjvh, (adjv - sklv).normalized());
				}
			}
			if (vecs_tmp.size() == 2) {
				double cos = vecs_tmp[0].second.dot(vecs_tmp[1].second);
				if (cos < most_orthogonal) {
					most_orthogonal = cos;
					vh_selected = sklvh;
					vecs_selected = vecs_tmp;
				}
			}
		}

		if (vh_selected != -1) {
			place_branch_hex(vh_selected, vecs_selected);
			algorithm_insert = vh_selected;
		}

	}

	if (only_place__not_split) return;

	optimize_face_branch();// 同时也要细分曲线上的分支结点

	/*MeshKernel::IO io_output;
	io_output.iGameWriteVolumeFile(hexmesh, "C:/My Files/Graphics/model_data/!test_data/hexmesh_output.mesh");*/


}

void HexMeshing_Branch_Processor::init_nodes_map() {

	for (auto& ep : sklmesh.alledges()) {
		if (sklmesh.NeighborFh(ep.first).empty()) {
			curve_skel_ehs[ep.first] = true;// 曲线边的度为0
		}
	}

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

		if (faces_count != 0 || hanging_ehs.size() > 2 || is_corner(vp.first)) {// 线面分支 或者 线线分支 或者 线线角点
			auto& v = sklmesh.vertices(vp.first);
			Vector3d src(v.x(), v.y(), v.z());
			for (auto& adjeh : hanging_ehs) {// 只需选择悬挂边作为射线
				auto& e = sklmesh.edges(adjeh);
				VH adjvh(e.vh1() + e.vh2() - vp.first);
				auto& adjv = sklmesh.vertices(adjvh);
				Vec dir = (adjv - v).normalized();// 射线由分支节点射向外面
				
				Ray ray(src, Vector3d(dir.x(), dir.y(), dir.z()));
				int rayid = rays.size();
				rays.emplace_back(ray);
				sklvh_to_ray[vp.first].emplace_back(rayid);// 一个分支顶点可能会发射出多条射线

			}
			//std::cout << "SKL_VH = " << vp.first << ", Ray size = " << sklvh_to_ray[vp.first].size() << std::endl;


			if (faces_count == 0) {// 线线分支结点
				branch_curvh[vp.first] = true;
			}
		}

		

	}

	rayid_to_hexfh.resize(rays.size());// 每个光线和与之相交的六面体面
	rayid_to_intersection.resize(rays.size());// 每个光线和与之相交的交点

	std::cout << "Total branch size = " << sklvh_to_ray.size() << std::endl;
	std::cout << "Pure curve branch size = " << branch_curvh.size() << std::endl;

}

void HexMeshing_Branch_Processor::optimize_face_branch() {

	
	while (1) {

		rebuild_bvh_tree();
		recast_rays();// 重新发射所有射线
		auto multi_intersect_faces = check_unique_intersect();// 获取有多个交点的面

		if (multi_intersect_faces.empty()) break;

		std::cout << "We find " << multi_intersect_faces.size() << " multi intersect faces.\n";

		for (auto& mul_fh : multi_intersect_faces) {
			if (!hexmesh.isValid(mul_fh)) continue;// 面无效, 可能已经被细分掉了
			auto& face = hexmesh.faces(mul_fh);
			// 每个面就两个细分方向
			const auto& f_vhs = face.getVertexHandle();
			vector<Vex> f_vs = { hexmesh.vertices(f_vhs[0]), hexmesh.vertices(f_vhs[1]), 
				hexmesh.vertices(f_vhs[2]), hexmesh.vertices(f_vhs[3])};
			Vec vec01 = f_vs[1] - f_vs[0], dir01 = vec01.normalized();// 方向 1
			Vec vec03 = f_vs[3] - f_vs[0], dir03 = vec03.normalized();// 方向 2
			Vec normal_face = (dir01.cross(dir03)).normalized();// 面的法向

			// 找出距离最近的两个交点
			double dist_min = 999999999.f;
			int rayid0 = -1, rayid1 = -1;
			auto& f_rays = hexfh_to_ray[mul_fh];
			if (f_rays.size() < 2) continue;
			for (int i = 0; i < f_rays.size(); ++i) {
				int rid0 = f_rays[i];
				for (int j = i + 1; j < f_rays.size(); ++j) {
					int rid1 = f_rays[j];
					assert(rid1 != rid0);
					double dist = (rayid_to_intersection[rid0] - rayid_to_intersection[rid1]).norm2();
					if (dist < dist_min) {
						dist_min = dist;
						rayid0 = rid0;
						rayid1 = rid1;
					}
				}
			}

			if (rayid0 == -1 || rayid1 == -1 || dist_min == static_cast<double>(0)) continue;

			Vec vec_ist = rayid_to_intersection[rayid1] - rayid_to_intersection[rayid0];
			Vec dir_ist = vec_ist.normalized();
			
			

			double len01 = 0, len03 = 0;
			//bool rev01 = false, rev03 = false;// 是否需要翻转
			if (dir_ist.dot(dir01) >= 0.0) {
				len01 = vec_ist.dot(dir01);// 投影距离
			} else {
				len01 = vec_ist.dot(dir01 * -1);
				//rev01 = true;
			}
			if (dir_ist.dot(dir03) >= 0.0) {
				len03 = vec_ist.dot(dir03);
			} else {
				len03 = vec_ist.dot(dir03 * -1);
				//rev03 = true;
			}

			// 选择边长更长的方向作为细分方向
			Vex midpoint_ist = (rayid_to_intersection[rayid1] + rayid_to_intersection[rayid0]) * 0.5;// 两个交点的中点
			Vex point_nearest;
			if (len01 >= len03) {// 以边 eh01 细分
				dist_point_line_segment(midpoint_ist, f_vs[0], f_vs[1], point_nearest);
				double dist2v0 = (point_nearest - f_vs[0]).norm();
				double dist2v1 = (point_nearest - f_vs[1]).norm();
				double weight_left = dist2v1 / (dist2v0 + dist2v1);
				EH eh01 = hexmesh.getEdgeHandle(f_vhs[0], f_vhs[1]);
				auto& e01 = hexmesh.edges(eh01);
				if (e01.vh1() != f_vhs[0]) weight_left = 1 - weight_left;
				Subdivision_Edge::subdivision_edges_weighted(hexmesh, eh01, weight_left);

			} else {// 以边 eh03 细分
				dist_point_line_segment(midpoint_ist, f_vs[0], f_vs[3], point_nearest);
				double dist2v0 = (point_nearest - f_vs[0]).norm();
				double dist2v3 = (point_nearest - f_vs[3]).norm();
				double weight_left = dist2v3 / (dist2v0 + dist2v3);
				EH eh03 = hexmesh.getEdgeHandle(f_vhs[0], f_vhs[3]);
				auto& e03 = hexmesh.edges(eh03);
				if (e03.vh1() != f_vhs[0]) weight_left = 1 - weight_left;
				Subdivision_Edge::subdivision_edges_weighted(hexmesh, eh03, weight_left);

			}

		}


		//// 测试用, 只细分一次
		//break;
	}

	//// 任务3: 处理 线面分支 节点
 //   //        两个在不同体且相邻的分支面, 中间应该插入一个面作为缓冲地带, 这样在后面的拟合可以取得更好的效果
	//while (1) {

	//	rebuild_bvh_tree();
	//	recast_rays();// 重新发射所有射线

	//	bool is_ok = true;
	//	for (auto it1 = hexfh_to_ray.begin(); it1 != hexfh_to_ray.end() && is_ok; ++it1) {
	//		auto it2 = it1;
	//		++it2;
	//		auto fh1 = it1->first;
	//		const auto& chs1 = hexmesh.NeighborCh(fh1);

	//		while (it2 != hexfh_to_ray.end() && is_ok) {
	//			auto fh2 = it2->first;
	//			if (hexmesh.isConnected(fh1, fh2)) {
	//				bool in_same_cell = false;
	//				const auto& chs2 = hexmesh.NeighborCh(fh2);
	//				for (auto& ch2 : chs2) {
	//					if (chs1.count(ch2)) {
	//						in_same_cell = true;
	//						break;
	//					}
	//				}

	//				if (in_same_cell) {// 在同一个体上, 不用做细分处理
	//					it2++;
	//					continue;
	//				}
	//				// 找到相邻的那条边
	//				const auto& ehs1 = hexmesh.faces(fh1).getEdgeHandle();
	//				const auto& ehs2 = hexmesh.faces(fh2).getEdgeHandle();
	//				EH com_eh(-1);
	//				for (auto& eh1 : ehs1) {
	//					for (auto& eh2 : ehs2) {
	//						if (eh1 == eh2) {
	//							com_eh = eh1;
	//							break;
	//						}
	//					}
	//					if (com_eh != -1) break;
	//				}
	//				if (com_eh != -1) {
	//					Vex midpoint_ist = (rayid_to_intersection[it1->second[0]] + hexmesh.getEdgeMidpoint(com_eh)) * 0.5;
	//					EH sub_eh(-1);
	//					for (auto& eh : ehs1) {
	//						if (eh != com_eh && hexmesh.isConnected(eh, com_eh)) {
	//							sub_eh = eh;
	//							break;
	//						}
	//					}
	//					if (sub_eh != -1) {
	//						Vex point_nearest;
	//						auto& sub_e = hexmesh.edges(sub_eh);
	//						dist_point_line_segment(midpoint_ist, hexmesh.vertices(sub_e.vh1()), hexmesh.vertices(sub_e.vh2()), point_nearest);
	//						double dist2v0 = (point_nearest - hexmesh.vertices(sub_e.vh1())).norm();
	//						double dist2v1 = (point_nearest - hexmesh.vertices(sub_e.vh2())).norm();
	//						double weight_left = dist2v1 / (dist2v0 + dist2v1);
	//						Subdivision_Edge::subdivision_edges_weighted(hexmesh, sub_eh, weight_left);
	//						is_ok = false;
	//					}
	//				}

	//			}
	//			it2++;

	//		}
	//		
	//	}

	//	if (is_ok) break;

	//}

	std::cout << "Optimize face branch success.\n" << std::endl;

}



void HexMeshing_Branch_Processor::optimize_curve_branch() {

	for (auto& cvh : branch_curvh) {
		VH sklvh = cvh.first;
		auto& sklv = sklmesh.vertices(sklvh);
		vector<pair<VH, Vec>> adj_vecs;
		const auto& adjvhs = sklmesh.NeighborVh(sklvh);

		for (auto& adjvh : adjvhs) {
			auto& adjv = sklmesh.vertices(adjvh);
			adj_vecs.emplace_back(adjvh, (adjv - sklv).normalized());
		}

		//double abs_cosine_min = 1;// 越接近垂直, abs_cosine 越接近 0
		//for (int i = 0; i < adj_vecs.size(); ++i) {
		//	for (int j = i + 1; j < adj_vecs.size(); ++j) {
		//		double abs_cosine = std::abs(adj_vecs[i].second.dot(adj_vecs[j].second));
		//		if (abs_cosine < abs_cosine_min) {
		//			abs_cosine_min = abs_cosine;
		//		}
		//	}
		//}

		place_branch_hex(sklvh, adj_vecs);
	}

}

vector<FH> HexMeshing_Branch_Processor::check_unique_intersect() {

	vector<FH> not_unique_faces;

	for (auto& fp : hexfh_to_ray) {
		// 需要剔除那些被细分掉的无效面
		if (hexmesh.isValid(fp.first) && fp.second.size() > 1) {
			not_unique_faces.emplace_back(fp.first);
		}
	}

	return not_unique_faces;

}

void HexMeshing_Branch_Processor::rebuild_bvh_tree() {

	delete bvh_tree;
	bvh_tree = nullptr;
	bvh_tree = new BVH_Tree();
	bvh_tree->buildBVH_Tree(hexmesh);

}

void HexMeshing_Branch_Processor::recast_rays() {

	hexfh_to_ray.clear();// 清除旧数据

	for (int i = 0; i < rays.size(); ++i) {
		recast_ray(i);
		//std::cout << "Ray " << i << " --> face " << rayid_to_hexfh[i] << std::endl;
	}

}

void HexMeshing_Branch_Processor::recast_ray(int rayid) {
	if (rayid >= rays.size()) {
		std::cerr << "Ray ID is invalid!\n";
		return;
	}

	const auto& ray = rays[rayid];

	auto intersection = bvh_tree->getIntersection(ray);
	if (intersection.happened == false) {
		std::cerr << "Ray no intersect with hex mesh.\n";
		return;
	}
	FH fh(intersection.index);
	rayid_to_hexfh[rayid] = fh;
	rayid_to_intersection[rayid] = Vex(intersection.pos[0], intersection.pos[1], intersection.pos[2]);
	hexfh_to_ray[fh].emplace_back(rayid);

}

double HexMeshing_Branch_Processor::dist_point_line_segment(Vex& v, Vex& v0, Vex& v1, Vex& nearest_vertex) {
	Vec vec1 = v - v0;
	Vec vec2 = v1 - v0;
	double t = vec2 * vec2;
	Vex min_v = v0;
	if (t > 1E-6F) {
		t = vec1 * vec2 / t;
		if (t > 1.f) {
			min_v = v1;
			vec1 = v - min_v;
		} else if (t > 0.f) {
			min_v = v0 + vec2 * t;
			vec1 = v - min_v;
		}
	}
	nearest_vertex = min_v;
	return vec1.norm();
}

void HexMeshing_Branch_Processor::place_branch_hex(VH vh, vector<pair<VH, Vec>>& adj_vecs) {

	/*if (adj_vecs.size() < 3) {
		std::cerr << "It is not a curve branch!!!\n";
		return;
	}*/

	

	vector<Eigen::Vector3d> input_pca;
	for (auto& adj_vec : adj_vecs) {
		auto& vec = adj_vec.second;
		input_pca.emplace_back(Eigen::Vector3d(vec.x(), vec.y(), vec.z()));
	}
	vector<Eigen::Vector3d> output_pca;
	Math_PCA::get_principal_components(input_pca, output_pca);

	/*for (auto& vec : input_pca) {
		std::cout << "PCA Input: " << vec << ", and norm = " << vec.norm() << std::endl;
	}
	for (auto& vec :output_pca) {
		std::cout << "PCA Output: " << vec << ", and norm = " << vec.norm() << std::endl;
	}*/

	/*if (output_pca.size() < 2) {
		std::cerr << "PCA result dimension is too low!!!\n";
	}*/
	Vec normal0_hex(output_pca[0].x(), output_pca[0].y(), output_pca[0].z());
	normal0_hex.normalize();
	Vec normal1_hex, normal2_hex;
	
	//std::cout << "Norm: Normal 0 = " << normal0_hex.norm() << ", Normal 1 = " << normal1_hex.norm() << std::endl;
	//std::cout << "Cosine between Normal 0 and Normal 1 = " << normal0_hex.dot(normal1_hex) << std::endl;
	//if (normal0_hex.norm() < 0.999) std::cerr << "Normal 0 is invalid!\n";
	//if (normal1_hex.norm() < 0.999) std::cerr << "Normal 0 is invalid!\n";

	if (output_pca.size() == 2) {
		normal1_hex = Vec(output_pca[1].x(), output_pca[1].y(), output_pca[1].z());
		normal1_hex.normalize();
	} else {

		vector<double> rads(adj_vecs.size(), 0);// 越大说明方向越接近
		for (int i = 0; i < adj_vecs.size(); ++i) {
			Vec& vec0 = adj_vecs[i].second;
			for (int j = i + 1; j < adj_vecs.size(); ++j) {
				Vec& vec1 = adj_vecs[j].second;
				double cosine = vec0.dot(vec1);
				if (cosine < 0.f) cosine = vec0.dot(vec1 * -1);
				rads[i] += cosine;
				rads[j] += cosine;
			}
		}

		int rad_max_idx = 0;// 最适合作为六面体法向的一条射线
		for (int i = 1; i < adj_vecs.size(); ++i) {
			if (rads[i] > rads[rad_max_idx]) {
				rad_max_idx = i;
			}
		}
		normal1_hex = adj_vecs[rad_max_idx].second;
		normal2_hex = (normal0_hex.cross(normal1_hex)).normalized();
		normal1_hex = (normal0_hex.cross(normal2_hex)).normalized();

	}

	/*if (output_pca.size() == 3) {
		normal2_hex = Vec(output_pca[2].x(), output_pca[2].y(), output_pca[2].z());
		normal2_hex.normalize();
	} else {
		normal2_hex = (normal0_hex.cross(normal1_hex)).normalized();
	}*/

	// 改用暴力求解
	Eigen::Vector3d vec_pca(normal0_hex.x(), normal0_hex.y(), normal0_hex.z());
	Eigen::Vector3d vec_base(normal1_hex.x(), normal1_hex.y(), normal1_hex.z());
	//auto vecs_opt = optimize_orientation(vec_pca, vec_base, input_pca);
	auto vecs_opt = optimize_orientation_zmy(vec_pca, input_pca);
	normal0_hex = Vex(vecs_opt[0].x(), vecs_opt[0].y(), vecs_opt[0].z());
	normal1_hex = Vex(vecs_opt[1].x(), vecs_opt[1].y(), vecs_opt[1].z());
	normal2_hex = Vex(vecs_opt[2].x(), vecs_opt[2].y(), vecs_opt[2].z());

	vector<Vec> normals{ normal0_hex, normal1_hex, normal2_hex };
	//// 测试, 轴对齐box
	//normals = { Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1) };
	vector<vector<int>> dirs{ {-1, -1, -1}, { 1, -1, -1 }, {1, -1, 1}, {-1, -1, 1},
		{-1, 1, -1}, { 1, 1, -1 }, {1, 1, 1}, {-1, 1, 1} };
	vector<VH> hex_vhs;
	vector<Vex> hex_vs;
	auto& sklv = sklmesh.vertices(vh);
	for (auto& dir : dirs) {
		Vec move_dir = (normals[0] * dir[0] + normals[1] * dir[1] + normals[2] * dir[2]).normalized();
		Vex hexv = sklv + move_dir * sklmesh.radius[vh] * radius_scale / 1.414;
		VH hexvh = hexmesh.AddVertex(hexv);
		hex_vhs.push_back(hexvh);
		hex_vs.push_back(hexv);
	}

	CH hexch = hexmesh.AddCell(hex_vhs);

}

vector<Eigen::Vector3d> HexMeshing_Branch_Processor::optimize_orientation(Eigen::Vector3d vec_pca, Eigen::Vector3d vec_base, 
	vector<Eigen::Vector3d>& dirs) {

	vec_pca.normalize();
	vec_base.normalize();
	// 暴力遍历旋转180°, 找出局部最优解
	double angle2rad = PI / 180.0;
	
	Eigen::Vector3d u = vec_pca, v, w;
	vector<Eigen::Vector3d> result;
	double target_z_min = (double)std::numeric_limits<double>::max();
	for (int it = 0; it < 15; ++it) {
		for (int it_r = 0; it_r < 1; ++it_r) {
			// 固定 u, 优化 v
			double nx = u.x(), ny = u.y(), nz = u.z();
			for (int i = 0; i < 90; ++i) {
				double rad = angle2rad * i;
				Eigen::Matrix3d mat_rotate;
				mat_rotate << nx * nx * (1 - cos(rad)) + cos(rad), nx* ny* (1 - cos(rad)) + nz * sin(rad), nx* nz* (1 - cos(rad)) - ny * sin(rad),
					nx* ny* (1 - cos(rad)) - nz * sin(rad), ny* ny* (1 - cos(rad)) + cos(rad), ny* nz* (1 - cos(rad)) + nx * sin(rad),
					nx* nz* (1 - cos(rad)) + ny * sin(rad), ny* nz* (1 - cos(rad)) - nx * sin(rad), nz* nz* (1 - cos(rad)) + cos(rad);
				v = mat_rotate * vec_base;// vec_base 是被优化的向量
				w = (u.cross(v)).normalized();
				//std::cout << "v.norm() = " << v.norm() << std::endl;
				double target_z = 0;
				for (auto& dir : dirs) {
					target_z += (fabs(dir.dot(u)) + fabs(dir.dot(v)) + fabs(dir.dot(w)));
					//target_z += std::fmin(fabs(dir.dot(u)), std::fmin(fabs(dir.dot(v)), fabs(dir.dot(w))));// by zs
				}
				if (target_z < target_z_min) {
					target_z_min = target_z;
					result = { u, v, w };
				}
			}
			vec_base = w;
			u = v;
		}
	}
	
	std::cout << "[Hex Branch]: Dirs size = " << dirs.size() << ", target_z = " << target_z_min << std::endl;

	return result;

}

vector<Eigen::Vector3d> HexMeshing_Branch_Processor::optimize_orientation_zmy(Eigen::Vector3d vec_pca, vector<Eigen::Vector3d>& dirs) {

	Eigen::Vector3d local_z(vec_pca.x(), vec_pca.y(), vec_pca.z()), local_x;
	local_z.normalize();
	// find one local_x, local_y;
	if (local_z.x() == 0)
		local_x = Eigen::Vector3d(1, 0, 0);
	else if (local_z.y() == 0)
		local_x = Eigen::Vector3d(0, 1, 0);
	else if (local_z.z() == 0)
		local_x = Eigen::Vector3d(0, 0, 1);
	else
		local_x = Eigen::Vector3d(-local_z.y(), local_z.x(), 0).normalized();
	// lamda
	auto calValue = [&dirs](Eigen::Vector3d _local_z, Eigen::Vector3d _local_x) {
		Eigen::Vector3d _local_y = (_local_x.cross(_local_z)).normalized();
		double ans = 0;
		for (const auto& toward : dirs)
			ans += fabs(_local_x.dot(toward)) + fabs(_local_y.dot(toward)) + fabs(_local_z.dot(toward));
		return ans;
	};
	double all_min = DBL_MAX;
	Eigen::Vector3d all_local_x, all_local_z;
	int iter = 20;
	while (iter--) {
		// 固定 local_z, 优化 local_x, 即 local_x 绕着 local_z 旋转 
		glm::vec3 RotationAxis(local_z.x(), local_z.y(), local_z.z());
		Eigen::Vector3d ans_local_x = local_x;// 记录本次迭代的最佳向量
		double min = calValue(local_z, local_x);// 记录本次迭代的最小点乘值
		for (int angle = 1; angle < 90; angle++) {
			glm::quat q = inAngleAxis(RotationAxis, 1.0);
			rotateByQuat(q, local_x);
			double now = calValue(local_z, local_x);
			if (min > now) {
				min = now;
				ans_local_x = local_x;
			}
		}
		//local_x = local_z;  // 二维转换
		local_x = (ans_local_x.cross(local_z)).normalized();
		local_z = ans_local_x;
		// iter min
		if (all_min > min) {
			all_min = min;
			all_local_x = local_x;
			all_local_z = local_z;
		}
	}
	auto all_local_y = (all_local_x.cross(all_local_z)).normalized();
	return { all_local_x, all_local_y, all_local_z };
}

bool HexMeshing_Branch_Processor::is_corner(VH vh) {

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

glm::quat HexMeshing_Branch_Processor::inAngleAxis(glm::vec3 RotationAxis, double RotationAngle) {
	RotationAngle = RotationAngle * M_PI / double(180.0);
	RotationAxis = normalize(RotationAxis);
	glm::quat t;
	t.x = RotationAxis.x * sin(RotationAngle / 2);
	t.y = RotationAxis.y * sin(RotationAngle / 2);
	t.z = RotationAxis.z * sin(RotationAngle / 2);
	t.w = cos(RotationAngle / 2);
	return t;
}

void HexMeshing_Branch_Processor::rotateByQuat(const glm::quat& q, Eigen::Vector3d& in) {
	glm::mat4 model = glm::mat4(1.0f);
	model = glm::mat4_cast(q) * model; // 旋转模型矩阵
	glm::vec4 p0(in.x(), in.y(), in.z(), 0);
	glm::vec4 out = model * p0;
	in = Eigen::Vector3d(out.x, out.y, out.z);
}