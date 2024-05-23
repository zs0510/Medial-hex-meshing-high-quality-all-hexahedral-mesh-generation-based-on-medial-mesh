#include "QuadMesher.h"

bool QuadMesher::meshing(MeshKernel::SurfaceMesh& quadmesh) {

	if (!preprocessing()) return false;

	calculate_cross_field();
	
	generate_initial_streamlines();
	simplify_streamlines();

	return true;

}

bool QuadMesher::preprocessing() {

	if (!trimesh.isTriangleMesh()) {
		std::cerr << "QuadMesher: Input mesh isn't triangle mesh!" << std::endl;
		return false;
	}

	for (auto& vp : trimesh.allvertices()) {
		if (std::abs(vp.second.z()) > DBL_EPSILON) {// 非二维
			std::cerr << "QuadMesher: Input mesh isn't a 2D plane." << std::endl;
			return false;
		}
	}

	return true;

}

void QuadMesher::quad_meshing() {



}

void QuadMesher::simplify_streamlines() {

	int before_sl_count = streamlines.size();
	std::cout << "We are going to merge!\n";
	merge_neighbor_edgehandle();

	// 记录奇异点流出与流进的流线 id
	std::unordered_map<int, std::vector<int>> out_streamline;
	std::unordered_map<int, std::vector<int>> in_streamline;
	// 1. 参数化每条流线
	for (int i = 0; i < streamlines.size(); ++i) {
		auto& sl = streamlines[i];
		sl.init_length_ratio();
		out_streamline[sl.beg_sid].push_back(i);
		if (!sl.end_on_boundary) in_streamline[sl.end_sid].push_back(i);
	}

	// 2. 通过拓扑和几何来简化流线
	std::vector<bool> need_to_erase(streamlines.size(), false);// 记录旧流线是否可直接作为新流线使用
	std::vector<Streamline> streamlines_new;

	// 2.1 简化端点相同的流线
	for (int i = 0; i < streamlines.size(); ++i) {
		auto& sl0 = streamlines[i];
		if (need_to_erase[i] || sl0.end_on_boundary || sl0.beg_sid == sl0.end_sid) continue;
		int beg_sid = sl0.beg_sid, end_sid = sl0.end_sid;
		
		std::vector<std::pair<double, int>> candiate_slid;

		candiate_slid.push_back( { sl0.length.back(), i } );

		for (int j = i + 1; j < streamlines.size(); ++j) {
			auto& sl1 = streamlines[j];
			if (need_to_erase[j] || sl1.end_on_boundary) continue;
			if (sl1.beg_sid == end_sid && sl1.end_sid == beg_sid) {// 找到两条意义相同, 方向相反的流线
				candiate_slid.push_back( { sl1.length.back(), j } );
			} else if (sl1.beg_sid == beg_sid && sl1.end_sid == end_sid) {// 找到两条意义相同, 方向相同的流线
				candiate_slid.push_back({ sl1.length.back(), j });
			}
		}

		if (candiate_slid.size() > 1) {
			std::sort(candiate_slid.begin(), candiate_slid.end(), [&](std::pair<double, int>& s0, std::pair<double, int>& s1) {
				return s0.first > s1.first;
				});
			for (auto& slid : candiate_slid) need_to_erase[slid.second] = true;

			int slid0 = candiate_slid[0].second;
			int slid1 = candiate_slid[1].second;
			auto& sl_0 = streamlines[slid0];
			auto& sl_1 = streamlines[slid1];
			if (sl_0.end_sid == sl_1.beg_sid) {
				auto& points = sl_1.points;
				std::reverse(points.begin(), points.end());
				sl_1.beg_sid = sl_0.beg_sid, sl_1.end_sid = sl_0.end_sid;
				sl_1.init_length_ratio();
			}

			Streamline sl_new;
			get_streamline_itrp(sl_0, sl_1, sl_new);
			streamlines_new.push_back(sl_new);
		}

	}

	for (int i = 0; i < streamlines.size(); ++i) {// 擦去冗余的流线
		if (!need_to_erase[i]) {
			streamlines_new.push_back(streamlines[i]);
		}
	}
	
	streamlines = streamlines_new;
	int after_sl_count = streamlines.size();
	std::cout << "Streamlines simplification: Before = " << before_sl_count << ", after = " << after_sl_count << "\n";

	// 将边界组织成流线
	

}

void QuadMesher::generate_initial_streamlines() {

	init_singular_neighbor_handles();// 将奇异点邻域的流线延伸至奇异点

	for (auto& singular_point : singular_points) {

		if (singular_point.mode == 0) generate_initial_streamlines_v(singular_point);
		else if (singular_point.mode == 1) generate_initial_streamlines_e(singular_point);
		else if (singular_point.mode == 2) generate_initial_streamlines_f(singular_point);

	}

	std::vector<Streamline> new_slines;
	for (int i = 0; i < streamlines.size(); ++i) {
		if (streamlines[i].end_on_fault) bad_streamlines.push_back(streamlines[i]);
		else new_slines.push_back(streamlines[i]);
	}
	streamlines = new_slines;

	std::cout << "Generate streamlines success.\n";
	for (auto& sl : streamlines) {
		std::cout << "Points size = " << sl.points.size() << "\n";
	}

}

bool QuadMesher::get_first_point_streamline(Vec source, EH eh, FH fh, std::vector<Vec>& march_dir, std::vector<double>& weight_vh1) {// 源点, 目标边, 所在面
	auto& edge = trimesh.edges(eh);
	const auto& fvhs = trimesh.faces(fh).getVertexHandle();
	VH vh1 = edge.vh1();
	VH vh2 = edge.vh2();
	auto& v1 = trimesh.vertices(vh1);
	auto& v2 = trimesh.vertices(vh2);
	
	int fvid1 = -1, fvid2 = -1;
	for (int i = 0; i < 3; ++i) {
		if (fvhs[i] == vh1) fvid1 = i;
		else if (fvhs[i] == vh2) fvid2 = i;
	}
	if (fvid1 == -1 || fvid2 == -1) return false;// 边面不一致
	std::vector<Vec> cross_edge_dirs(2);// 保存跟最垂直于边的标架分量
	std::vector<std::vector<Vec>> cross_dirs = {// 获得边端点的标架
		get_cross_dirs(crossfield[fh][fvid1]) , get_cross_dirs(crossfield[fh][fvid2])
	};

	Vec dire = (v2 - v1).normalize();
	Vec dirz(0, 0, 1);
	Vec dir_t_e = (dire.cross(dirz)).normalized();// 垂直于边的方向

	for (int i = 0; i < 2; ++i) {
		double cos_max = -2;
		for (auto& dirv : cross_dirs[i]) {
			double _cos = dir_t_e.dot(dirv);
			if (_cos > cos_max) {
				cos_max = _cos;
				cross_edge_dirs[i] = dirv;
			}
		}
	}

	for (double w = 0.f; w <= 1.f; w += 0.001) {
		Vec point = v1 * w + v2 * (1 - w);
		Vec sp = (point - source).normalize();
		Vec dirp = cross_edge_dirs[0] * w + cross_edge_dirs[1] * (1 - w);// 标架分量
		dirp.normalize();
		std::vector<Vec> cross_p = get_cross_dirs(dirp);// 四个标架
		for (auto& dir : cross_p) {
			double cosine = dir.dot(sp);
			if (cosine > 0.99) {
				weight_vh1.push_back(w);
				march_dir.push_back(dir);
				w += 0.3;
				break;
			}
		}
	}
	
	return !weight_vh1.empty();
}

std::vector<Vec> QuadMesher::get_cross_dirs(const Vec& point, EH eh, FH fh) {// 返回 point 在边上插值得到的标架场 
	auto& edge = trimesh.edges(eh);
	auto& face = trimesh.faces(fh);
	auto vh1 = edge.vh1(), vh2 = edge.vh2();
	auto& v1 = trimesh.vertices(vh1);
	auto& v2 = trimesh.vertices(vh2);
	const auto& fvhs = face.getVertexHandle();
	int fvid1 = -1, fvid2 = -1;
	for (int i = 0; i < 3; ++i) {
		if (fvhs[i] == vh1) fvid1 = i;
		else if (fvhs[i] == vh2) fvid2 = i;
	}
	if (fvid1 == -1 || fvid2 == -1) {
		std::cerr << "QuadMesher: EH is not belong to FH.\n";
		return {};
	}
	std::vector<Vec> cross_edge_dirs(2);
	std::vector<std::vector<Vec>> cross_dirs = {
		get_cross_dirs(crossfield[fh][fvid1]) , get_cross_dirs(crossfield[fh][fvid2])
	};
	for (int i = 0; i < 2; ++i) {// 选出两个跟边方向最接近的标架分量
		Vec dire = (v2 - v1).normalize();
		double cos_max = -2;
		for (auto& dirv : cross_dirs[i]) {
			double _cos = dire.dot(dirv);
			if (_cos > cos_max) {
				cos_max = _cos;
				cross_edge_dirs[i] = dirv;
			}
		}
	}
	// 理论上使用横坐标之差即可
	double dis1 = (point - v1).norm();
	double dis2 = (point - v2).norm();
	Vec cross_dirp = cross_edge_dirs[0] * dis2 + cross_edge_dirs[1] * dis1;
	cross_dirp.normalize();
	return get_cross_dirs(cross_dirp);
}

std::vector<Vec> QuadMesher::get_cross_dirs(VH vh, FH fh) {

	auto& face = trimesh.faces(fh);
	const auto& vhs = face.getVertexHandle();
	int local_id = -1;
	if (vhs[0] == vh) local_id = 0;
	else if (vhs[1] == vh) local_id = 1;
	else if (vhs[2] == vh) local_id = 2;
	
	if (local_id == -1) {
		std::cout << "Error: vh is not belong to fh" << std::endl;
		return {};
	}

	return get_cross_dirs(crossfield[fh][local_id]);
}

std::vector<Vec> QuadMesher::get_cross_dirs(const Vec& vec) {
	std::vector<Vec> res(4);
	Vec dirz(0, 0, 1);
	res[0] = vec;
	res[1] = dirz.cross(vec);
	res[1].normalize();
	res[2] = res[0] * -1;
	res[3] = res[1] * -1;
	return res;
}

void QuadMesher::calculate_cross_field() {

	/* VH is equal to ID */

	std::vector<Vec3> points;
	std::vector<ID3> triangles;
	std::vector<ID2> feature_edges;
	std::vector<ID3> singularities;// 奇异点的 id
	std::vector<std::array<double, 9>> global_triangle_dir;

	// 求解标架场的基本数据
	points.reserve(trimesh.vsize());
	triangles.reserve(trimesh.fsize());

	for (int i = 0; i < trimesh.vsize(); ++i) {
		VH vh(i);
		auto& v = trimesh.vertices(vh);
		points.push_back(Vec3{ v.x(), v.y(), v.z() });
	}

	for (int i = 0; i < trimesh.fsize(); ++i) {
		FH fh(i);
		auto& f = trimesh.faces(fh);
		const auto& vhs = f.getVertexHandle();
		triangles.push_back(ID3{ (ID)vhs[0], (ID)vhs[1], (ID)vhs[2] });
	}

	// 求解标架场的约束数据
	for (auto& eh : feature_ehs) {
		auto& e = trimesh.edges(eh);
		feature_edges.push_back(ID2{ (ID)e.vh1(), (ID)e.vh2() });
	}

	std::cout << "Calculate cross field - Input Data: #V = " << points.size() << ", #T = " << triangles.size() << "; constrained #E = " << feature_edges.size() << std::endl;

	MeshMath::calc_cross_field(points, triangles, feature_edges, singularities, global_triangle_dir);

	int sv_cnt = 0, se_cnt = 0, sf_cnt = 0;
	SingularPoint point;
	// 通过求解标架场得到的奇异点
	for (auto& singular_id : singularities) {// 所有的奇异点
		if (singular_id[0] != NO_ID) {// 是顶点奇异点
			auto& v = trimesh.vertices(VH(singular_id[0]));
			point.pos = v;
			point.mode = 0;
			point.handle = singular_id[0];
			point.id = singular_points.size();
			vh_to_spid[VH(singular_id[0])] = singular_points.size();
			sv_cnt++;
		} else if (singular_id[1] != NO_ID) {// 是边奇异点
			int fidx = singular_id[1] / 3;// 将 ID_SU 转换为 ID_IGAME
			int eidx = singular_id[1] % 3;
			FH fh(fidx);
			EH eh = trimesh.faces(fh).eh(eidx);
			auto& edge = trimesh.edges(eh);
			auto& v1 = trimesh.vertices(edge.vh1());
			auto& v2 = trimesh.vertices(edge.vh2());
			point.pos = (v1 + v2) * 0.5;
			point.mode = 1;
			point.handle = eh;
			point.id = singular_points.size();
			eh_to_spid[eh] = singular_points.size();
			se_cnt++;
		} else if (singular_id[2] != NO_ID) {// 是面奇异点
			FH fh(singular_id[2]);
			const auto& vhs = trimesh.faces(fh).getVertexHandle();
			auto& v0 = trimesh.vertices(vhs[0]);
			auto& v1 = trimesh.vertices(vhs[1]);
			auto& v2 = trimesh.vertices(vhs[2]);
			point.pos = (v0 + v1 + v2) / 3.f;// 取面片中心
			point.mode = 2;
			point.handle = fh;
			point.id = singular_points.size();
			fh_to_spid[fh] = singular_points.size();
			sf_cnt++;
		}
		singular_points.push_back(point);
	}

	for (int i = 0; i < trimesh.fsize(); ++i) {
		crossfield.push_back({
			Vec(global_triangle_dir[i][0], global_triangle_dir[i][1], global_triangle_dir[i][2]),
			Vec(global_triangle_dir[i][3], global_triangle_dir[i][4], global_triangle_dir[i][5]),
			Vec(global_triangle_dir[i][6], global_triangle_dir[i][7], global_triangle_dir[i][8])
			});
	}

	//// unit check
	//for (auto& face_cross : crossfield) {
	//	for (auto& vex_cross : face_cross) {
	//		if (std::abs(vex_cross.norm2() - 1.f) > 1e-6f) std::cerr << "Error: No normalize!\n";
	//	}
	//}

	std::cout << "QuadMesher: Calculate cross field success." << std::endl;
	std::cout << "\tSingular vertices size = " << sv_cnt << ".\n";
	std::cout << "\tSingular edges size = " << se_cnt << ".\n";
	std::cout << "\tSingular faces size = " << sf_cnt << ".\n";

	int cv_cnt = 0;
	for (auto& vp : trimesh.allvertices()) {
		if (trimesh.isOnBoundary(vp.first)) {
			EH eh0(-1), eh1(-1);// 两条边界边
			for (auto& adjeh : trimesh.NeighborEh(vp.first)) {
				if (trimesh.isOnBoundary(adjeh)) {
					if (eh0 == -1) eh0 = adjeh;
					else {
						eh1 = adjeh;
						break;
					}
				}
			}
			if (eh0 == -1 || eh1 == -1) {
				std::cerr << "Boundary vh no 2 boundary ehs.\n";
				continue;
			}
			auto& edge0 = trimesh.edges(eh0);
			auto& edge1 = trimesh.edges(eh1);
			auto& adjv0 = trimesh.vertices(VH(edge0.vh1() + edge0.vh2() - vp.first));
			auto& adjv1 = trimesh.vertices(VH(edge1.vh1() + edge1.vh2() - vp.first));
			auto& v = trimesh.vertices(vp.first);
			Vec vec0 = (adjv0 - v).normalize();
			Vec vec1 = (adjv1 - v).normalize();
			double cosine = vec0.dot(vec1);
			if (std::abs(cosine) < 0.5) {// cos85 == 0.087, 
				point.pos = v;
				point.mode = 0;// 属于顶点奇异点
				point.handle = vp.first;
				point.id = singular_points.size();
				vh_to_spid[vp.first] = singular_points.size();
				cv_cnt++;
				singular_points.push_back(point);
			}
		}
	}

	std::cout << "QuadMesher: Detect corner vertex success.\n\tCorner vertices size = " << cv_cnt << std::endl;
	

}

bool QuadMesher::get_intersection_line2(const Vex& src, const Vec& dir, EH eh, Vex& intersection) {

	double k0, b0, k1, b1;// 直线方程 y = kx + b
	bool line0_is_vertical = false, line1_is_vertical = false;
	double line0_x, line1_x;

	double numeric_error = 1E-6F;// 数值精度误差

	auto& edge = trimesh.edges(eh);
	auto& v1 = trimesh.vertices(edge.vh1());
	auto& v2 = trimesh.vertices(edge.vh2());
	// 求边的直线方程
	if (v2.x() != v1.x()) {
		k0 = (v2.y() - v1.y()) / (v2.x() - v1.x());
		b0 = v2.y() - k0 * v2.x();
	} else {
		line0_is_vertical = true;
		line0_x = v2.x();
	}
	// 求射线的直线方程
	if (dir.x() != 0.f) {
		k1 = dir.y() / dir.x();
		b1 = src.y() - k1 * src.x();
	} else {
		line1_is_vertical = true;
		line1_x = src.x();
	}

	if (!line0_is_vertical && !line1_is_vertical) {// 两条线都是斜的

		double diffk = k0 - k1;
		if (diffk == 0) return false;// 平行
		double x = (b1 - b0) / diffk;
		double y = k0 * x + b0;
		if (x < std::min(v1.x(), v2.x()) - numeric_error || 
			x > std::max(v1.x(), v2.x()) + numeric_error ||
			y < std::min(v1.y(), v2.y()) - numeric_error || 
			y > std::max(v1.y(), v2.y()) + numeric_error) {
			return false;// 交点在线段外
		}
		intersection = Vec(x, y, 0);

	} else if (line0_is_vertical && line1_is_vertical) {// 两条线都是竖直的

		if (std::abs(line0_x - line1_x) < numeric_error) {
			double dis1 = (v1 - src).norm();
			double dis2 = (v2 - src).norm();
			if (dis1 < dis2) intersection = v2;
			else intersection = v1;
		} else return false;

	} else if (line0_is_vertical) {

		double y = line0_x * k1 + b1;
		if (y < std::min(v1.y(), v2.y()) - numeric_error || 
			y > std::max(v1.y(), v2.y()) + numeric_error) return false;
		intersection = Vec(line0_x, y, 0);

	} else if (line1_is_vertical) {

		double y = line1_x * k0 + b0;
		if (y < std::min(v1.y(), v2.y()) - numeric_error ||
			y > std::max(v1.y(), v2.y()) + numeric_error) return false;
		intersection = Vec(line1_x, y, 0);

	}

	/*std::cout << "Intersection! source = (" << src.x() << ", " << src.y() << "), dir = (" << dir.x() << ", " << dir.y() << "), v1 = ("
		<< v1.x() << ", " << v1.y() << "), v2 = (" << v2.x() << ", " << v2.y() << ")" << std::endl;
	std::cout << "\tIntersection = (" << intersection.x() << ", " << intersection.y() << ")\n";
	std::cout << "\tLine0: " << intersection.y() - (intersection.x() * k0 + b0) << ", Line1: " << intersection.y() - (intersection.x() * k1 + b1) << "\n";*/
	return true;

}

Vec QuadMesher::get_most_similar_dir(const Vec& ref_dir, const std::vector<Vec>& dirs) {
	double cosine_max = -2;
	Vec res;
	for (auto& dir : dirs) {
		double cosine = ref_dir.dot(dir);
		if (cosine > cosine_max) {
			cosine_max = cosine;
			res = dir;
		}
	}
	if (cosine_max < 0.71) std::cout << "Warning: cosine is too small!\n";
	return res;
}

std::unordered_set<EH> QuadMesher::get_smallest_bounding_ehs(const EH& _eh) {
	std::unordered_set<EH> res;
	auto& edge = trimesh.edges(_eh);
	VH vh1 = edge.vh1(), vh2 = edge.vh2();
	std::unordered_set<FH> adjfhs = trimesh.NeighborFh(vh1);
	for (auto& adjfh : trimesh.NeighborFh(vh2)) adjfhs.insert(adjfh);
	for (auto& fh : adjfhs) {
		const auto& ehs = trimesh.faces(fh).getEdgeHandle();
		for (auto& eh : ehs) {
			//auto& e = trimesh.edges(eh);
			//if (e.vh1() == vh1 || e.vh2() == vh2 || e.vh1() == vh2 || e.vh2() == vh1) continue;// 不与 eh 邻接
			if (_eh == eh) continue;
			res.insert(eh);
		}
	}
	return res;
}

std::vector<std::pair<EH, FH>> QuadMesher::get_samllest_bounding_ehs_fhs(const EH& _eh) {
	std::vector<std::pair<EH, FH>> res;
	auto& edge = trimesh.edges(_eh);
	VH vh1 = edge.vh1(), vh2 = edge.vh2();
	std::unordered_set<FH> adjfhs = trimesh.NeighborFh(vh1);
	for (auto& adjfh : trimesh.NeighborFh(vh2)) adjfhs.insert(adjfh);
	for (auto& fh : adjfhs) {
		auto& face = trimesh.faces(fh);
		const auto& ehs = face.getEdgeHandle();
		for (auto& eh : ehs) {
			auto& e = trimesh.edges(eh);
			if (e.vh1() == vh1 || e.vh2() == vh1 || e.vh1() == vh2 || e.vh2() == vh1) continue;
			res.push_back(std::pair<EH, FH>(eh, fh));
			break;
		}
	}
	/*for (auto& fh : trimesh.NeighborFh(_eh)) {
		for (auto& eh : trimesh.faces(fh).getEdgeHandle()) {
			if (eh != _eh) res.push_back(std::pair<EH, FH>(eh, fh));
		}
	}*/
	return res;
}

void QuadMesher::get_streamline_itrp(const Streamline& sl0, const Streamline& sl1, Streamline& sl_new) {
	sl_new.beg_sid = sl0.beg_sid;
	sl_new.end_on_boundary = sl0.end_on_boundary && sl1.end_on_boundary;
	sl_new.end_sid = sl0.end_sid;
	int num = sl0.points.size() + sl1.points.size();
	if (num == 0) return;
	auto& points = sl_new.points;
	points.push_back((sl0.points[0] + sl1.points[0]) * 0.5f);
	double step = 1.f / num;
	int pid0 = 1, pid1 = 1;
	for (double r = step; r <= 1.f; r += step) {
		while (pid0 < sl0.ratio.size() && sl0.ratio[pid0] < r) pid0++;// 找到大于等于当前 r 的点
		while (pid1 < sl1.ratio.size() && sl1.ratio[pid1] < r) pid1++;
		if (pid0 == sl0.ratio.size() || pid1 == sl1.ratio.size()) break;
		double rdiff_00 = r - sl0.ratio[pid0 - 1];
		double rdiff_01 = sl0.ratio[pid0] - r;
		if (rdiff_00 < 0 || rdiff_01 < 0) {
			std::cerr << "QuadMesher: Wrong ratio!" << std::endl;
			break;
		}
		if (rdiff_00 + rdiff_01 == 0) {
			std::cerr << "QuadMesher: Devide zero!" << std::endl;
			break;
		} 
		Vex pos0 = (sl0.points[pid0 - 1] * rdiff_01 + sl0.points[pid0] * rdiff_00) / (rdiff_00 + rdiff_01);
		double rdiff_10 = r - sl1.ratio[pid1 - 1];
		double rdiff_11 = sl1.ratio[pid1] - r;
		if (rdiff_10 < 0 || rdiff_11 < 0) {
			std::cerr << "QuadMesher: Wrong num!" << std::endl;
			break;
		}
		if (rdiff_10 + rdiff_11 == 0) {
			std::cerr << "QuadMesher: Devide zero!" << std::endl;
			break;
		}
		Vex pos1 = (sl1.points[pid1 - 1] * rdiff_11 + sl1.points[pid1] * rdiff_10) / (rdiff_10 + rdiff_11);
		points.push_back((pos0 + pos1) * 0.5f);
	}
	if (sl_new.end_on_boundary) {// 将端点延伸至边界
		int n = sl_new.points.size();
		Vec dir_end = sl_new.points[n-1] - sl_new.points[n-2];
		dir_end.normalize();
		Vex intersection;
		EH eh0(sl0.end_edgehandle), eh1(sl1.end_edgehandle);
		std::set<EH> bdehs;
		get_among_ehs(eh0, eh1, bdehs);
		bool intersection_flag = false;
		for (auto& bdeh : bdehs) {
			if (get_intersection_line2(sl_new.points.back(), dir_end, bdeh, intersection)) {
				intersection_flag = true;
				sl_new.end_edgehandle = bdeh;
				break;
			}
		}
		if (!intersection_flag) {
			std::cerr << "QuadMesher: Boundary ehs not intersection with the ray.\n";
			return;
		}
		sl_new.points.push_back(intersection);
	} else {// 确保端点在奇异点上
		Vex endpoint = (sl0.points.back() + sl1.points.back()) * 0.5;
		if ((sl_new.points.back() - endpoint).norm() > 1E-6F) {
			sl_new.points.push_back(endpoint);
		}
	}
	sl_new.init_length_ratio();

}

void QuadMesher::get_among_ehs(const EH& beg, const EH& end, std::set<EH>& res) {
	if (beg == end) {
		res = { beg };
		return;
	}
	if (!trimesh.isOnBoundary(beg) || !trimesh.isOnBoundary(end)) {
		std::cerr << "QuadMesher: There exist eh not on the boundary.\n";
		return;
	}
	std::set<EH> ehs0, ehs1;
	ehs0.insert(beg), ehs1.insert(beg);
	EH eh0(-1), eh1(-1);
	for (auto& eh : trimesh.NeighborEh(beg)) {
		if (trimesh.isOnBoundary(eh)) {
			if (eh0 == -1) eh0 = eh;
			else if (eh1 == -1) eh1 = eh;
			else {
				std::cerr << "QuadMesher: It is not a 2D-manifold mesh.\n";
				return;
			}
		}
	}
	while (1) {
		ehs0.insert(eh0);
		if (eh0 == end) {
			res = ehs0;
			break;
		}

		ehs1.insert(eh1);
		if (eh1 == end) {
			res = ehs1;
			break;
		}

		int degree = 0;
		EH eh_tmp(-1);

		for (auto& eh : trimesh.NeighborEh(eh0)) {
			if (trimesh.isOnBoundary(eh)) {
				if (!ehs0.count(eh)) {
					eh_tmp = eh;
					degree++;
				}
			}
		}
		if (degree != 1 || eh_tmp == -1) {
			std::cerr << "QuadMesher: It is not a 2D-manifold mesh.\n";
			return;
		}
		eh0 = eh_tmp;

		degree = 0;
		eh_tmp = EH(-1);
		for (auto& eh : trimesh.NeighborEh(eh1)) {
			if (trimesh.isOnBoundary(eh)) {
				if (!ehs1.count(eh)) {
					eh_tmp = eh;
					degree++;
				}
			}
		}
		if (degree != 1 || eh_tmp == -1) {
			std::cerr << "QuadMesher: It is not a 2D-manifold mesh.\n";
			return;
		}
		eh1 = eh_tmp;

	}
}


void QuadMesher::generate_initial_streamlines_v(const SingularPoint& point) {

	VH vh(point.handle);

	std::vector<FH> fhs;
	std::vector<EH> ehs;// 最小多边形的边
	for (auto& adjfh : trimesh.NeighborFh(vh)) {
		const auto& f_ehs = (trimesh.faces(adjfh)).getEdgeHandle();
		for (auto& eh : f_ehs) {
			auto& edge = trimesh.edges(eh);
			if (edge.vh1() != vh && edge.vh2() != vh) {
				ehs.emplace_back(eh);
				break;
			}
		}
		fhs.emplace_back(adjfh);
	}
	assert(fhs.size() == ehs.size());
	auto& v = trimesh.vertices(vh);
	std::vector<MarchNode> nodes;
	for (int i = 0; i < fhs.size(); ++i) {
		std::vector<double> weight_vh1;
		std::vector<Vec> march_dir;
		if (get_first_point_streamline(v, ehs[i], fhs[i], march_dir, weight_vh1)) {
			auto& edge = trimesh.edges(ehs[i]);
			MarchNode tmp_node;
			for (int j = 0; j < march_dir.size(); ++j) {
				tmp_node.dir = march_dir[j];
				tmp_node.pre_fh = fhs[i];
				tmp_node.eh = ehs[i];
				tmp_node.pos = (trimesh.vertices(edge.vh1())) * weight_vh1[j] + (trimesh.vertices(edge.vh2())) * (1 - weight_vh1[j]);
				nodes.push_back(tmp_node);
			}
			std::cout << "get first points size = " << march_dir.size() << std::endl;
		}
	}

	std::vector<Vec> bdy_vec;
	for (auto& eh : trimesh.NeighborEh(vh)) {
		if (trimesh.isOnBoundary(eh)) {
			auto& edge = trimesh.edges(eh);
			auto& adjvh = edge.vh1() == vh ? edge.vh2() : edge.vh1();
			auto& adjv = trimesh.vertices(adjvh);
			Vec out_vec = (adjv - v).normalized();
			bdy_vec.push_back(out_vec);
		}
	}
	

	for (auto& node : nodes) {

		// 检查是否与边界基本平行
		bool is_parallel = false;
		Vec dir = (node.pos - v).normalized();
		for (auto& bdy_dir : bdy_vec) {
			if (dir.dot(bdy_dir) > 0.9) {
				is_parallel = true;
				break;
			}
		}
		if (is_parallel) continue;// 边界奇异点发射的、与边界基本平行的流线，不需要计算，后续直接组织边界即可

		Streamline sline;
		// 设置起点
		sline.beg_sid = point.id;
		sline.points.push_back(v);
		

		std::unordered_map<FH, bool> visited;
		lengthen_initial_streamlines(node, sline, visited);
		//std::cout << "Streamline_V #" << streamlines.size() << ": #P = " << sline.points.size() << std::endl;
		streamlines.push_back(sline);
	}

}

void QuadMesher::generate_initial_streamlines_e(const SingularPoint& point) {

	EH eh(point.handle);
	auto& edge = trimesh.edges(eh);
	std::unordered_set<FH> fhs = trimesh.NeighborFh(eh);
	
	Vex pos = point.pos;

	std::vector<MarchNode> nodes;
	for (auto& fh : fhs) {
		for (auto& adjeh : trimesh.faces(fh).getEdgeHandle()) {
			if (adjeh == eh) continue;
			std::vector<double> weight_vh1;
			std::vector<Vec> march_dir;
			if (get_first_point_streamline(pos, adjeh, fh, march_dir, weight_vh1)) {
				auto& edge = trimesh.edges(adjeh);
				MarchNode tmp_node;
				for (int j = 0; j < march_dir.size(); ++j) {
					tmp_node.dir = march_dir[j];
					tmp_node.pre_fh = fh;
					tmp_node.eh = adjeh;
					tmp_node.pos = (trimesh.vertices(edge.vh1())) * weight_vh1[j] + (trimesh.vertices(edge.vh2())) * (1 - weight_vh1[j]);
					nodes.push_back(tmp_node);
				}
				std::cout << "get first points size = " << march_dir.size() << std::endl;
			}

		}
	}

	for (auto& node : nodes) {
		Streamline sline;
		sline.beg_sid = point.id;
		sline.points.push_back(point.pos);
		/*std::set<EH> visited;
		morch_on_streamline(node, sline, visited);*/
		std::unordered_map<FH, bool> visited;
		lengthen_initial_streamlines(node, sline, visited);
		//std::cout << "Streamline_E #" << streamlines.size() << ": #P = " << sline.points.size() << std::endl;
		streamlines.push_back(sline);
	}
	
}

void QuadMesher::generate_initial_streamlines_f(const SingularPoint& point) {

	FH fh(point.handle);
	Vex pos = point.pos;
	int degree = 0;
	std::vector<MarchNode> nodes;

	for (auto& eh : trimesh.faces(fh).getEdgeHandle()) {

		std::vector<double> weight_vh1;
		std::vector<Vec> march_dir;
		if (get_first_point_streamline(pos, eh, fh, march_dir, weight_vh1)) {
			auto& edge = trimesh.edges(eh);
			MarchNode tmp_node;
			for (int j = 0; j < march_dir.size(); ++j) {
				tmp_node.dir = march_dir[j];
				tmp_node.pre_fh = fh;
				tmp_node.eh = eh;
				tmp_node.pos = (trimesh.vertices(edge.vh1())) * weight_vh1[j] + (trimesh.vertices(edge.vh2())) * (1 - weight_vh1[j]);
				nodes.push_back(tmp_node);
			}
			std::cout << "get first points size = " << march_dir.size() << std::endl;
		}

	}

	for (auto& node : nodes) {
		Streamline sline;
		sline.beg_sid = point.id;
		sline.points.push_back(point.pos);
		/*std::set<EH> visited;
		morch_on_streamline(node, sline, visited);*/
		std::unordered_map<FH, bool> visited;
		lengthen_initial_streamlines(node, sline, visited);
		//std::cout << "Streamline_F #" << streamlines.size() << ": #P = " << sline.points.size() << std::endl;
		streamlines.push_back(sline);
	}

}

void QuadMesher::init_singular_neighbor_handles() {

	singular_neighbor_fhs.clear();

	for (int i = 0; i < singular_points.size(); ++i) {
		auto& sp = singular_points[i];
		std::vector<VH> vhs;
		if (sp.mode == 0) {
			VH vh(sp.handle);
			vhs = { vh };

		} else if (sp.mode == 1) {
			EH eh(sp.handle);
			auto& edge = trimesh.edges(eh);
			vhs = { edge.vh1(), edge.vh2() };
			
		} else if (sp.mode == 2) {
			FH fh(sp.handle);
			auto& face = trimesh.faces(fh);
			vhs = face.getVertexHandle();
			
		}

		for (auto& vh : vhs) {
			for (auto& fh : trimesh.NeighborFh(vh)) {
				singular_neighbor_fhs[fh] = i;
			}
		}

	}


}

void QuadMesher::lengthen_initial_streamlines(const MarchNode& node, Streamline& sline, std::unordered_map<FH, bool>& visited, int bdy_faces_count) {

	if (singular_neighbor_fhs.count(node.pre_fh)) {// 到达奇异点附近，终止
		int sid = singular_neighbor_fhs[node.pre_fh];
		if (sid != sline.beg_sid || sline.points.size() > 25) {// 循环一圈或者到达其他奇异点附近
			auto& sp = singular_points[sid];
			sline.points.push_back(sp.pos);
			sline.fhs.push_back(node.pre_fh);
			sline.end_sid = sid;
			if (sline.beg_sid == sline.end_sid) std::cout << "Warning: Sline beg_sid == end_sid! ";
			std::cout << "End at singularity " << sid << std::endl;
			return;
		}
	}

	if (trimesh.isOnBoundary(node.eh)) {// 到达边界边，终止
		
		sline.points.push_back(node.pos);
		sline.fhs.push_back(node.pre_fh);
		sline.end_on_boundary = true;
		sline.end_edgehandle = node.eh;
		//std::cout << "End at boundary "<< std::endl;
		return;
	}

	if (visited.count(node.pre_fh)) {
		std::cout << "Warning: visit same face " << std::endl;
	}

	// insert new point
	sline.points.push_back(node.pos);
	sline.fhs.push_back(node.pre_fh);
	visited[node.pre_fh] = true;

	const auto& e_fhs = trimesh.NeighborFh(node.eh);
	FH next_fh(node.pre_fh);// next fh, we use the cross field on this triangle
	for (auto& e_fh : e_fhs) {
		if (e_fh != next_fh) {
			next_fh = e_fh;
			break;
		}
	}

	if (next_fh == node.pre_fh) {// 只有一个邻接面
		std::cerr << "Streamline is lengthen to boundary!!!\n";
		sline.end_on_boundary = true;
		sline.end_edgehandle = node.eh;
		return;
	}

	std::vector<Vec> crossdirs_p = get_cross_dirs(node.pos, node.eh, next_fh);// 在新面片上插值求得标架信息
	Vec dir_fp = get_most_similar_dir(node.dir, crossdirs_p);

	Vex tmp_p;// 辅助交点
	EH tmp_eh(-1);
	for (auto& f_eh : trimesh.faces(next_fh).getEdgeHandle()) {
		if (f_eh == node.eh) continue;
		if (get_intersection_line2(node.pos, dir_fp, f_eh, tmp_p)) {
			tmp_eh = f_eh;
			break;
		}
	}
	if (tmp_eh == -1) {// 应是延伸到了某个顶点
		auto bounding_ef = get_samllest_bounding_ehs_fhs(node.eh);
		for (auto& b_ef : bounding_ef) {
			if (visited.count(b_ef.second)) continue;
			if (get_intersection_line2(node.pos, node.dir, b_ef.first, tmp_p)) {
				tmp_eh = b_ef.first;
				next_fh = b_ef.second;
				break;
			}
		}
		if (tmp_eh == -1) {
			std::cerr << "Error: Can not get intersection -- 0!\n";
			sline.end_on_fault = true;
			return;
		}
	}

	std::vector<Vec> crossdirs_tmp_p = get_cross_dirs(tmp_p, tmp_eh, next_fh);
	Vec dir_ftp = get_most_similar_dir(node.dir, crossdirs_tmp_p);

	Vec dir_next = (dir_fp + dir_ftp).normalized();

	/*if (dir_next.dot(dir_fp) < 0) {
		std::cerr << "Dot of dirs is less than zero!\n";
	}*/

	Vex next_pos;
	EH next_eh(-1);
	for (auto& f_eh : trimesh.faces(next_fh).getEdgeHandle()) {
		if (f_eh == node.eh) continue;
		if (get_intersection_line2(node.pos, dir_next, f_eh, next_pos)) {
			next_eh = f_eh;
			break;
		}
	}
	if (next_eh == -1) {
		auto bounding_ef = get_samllest_bounding_ehs_fhs(node.eh);
		for (auto& b_ef : bounding_ef) {
			if (visited.count(b_ef.second)) continue;
			if (get_intersection_line2(node.pos, dir_next, b_ef.first, next_pos)) {
				next_eh = b_ef.first;
				next_fh = b_ef.second;
				break;
			}
		}
		if (next_eh == -1) {
			std::cerr << "Error: Can not get intersection! -- 1\n";
			sline.end_on_fault = true;
			return;
		}
	}

	MarchNode next_node;
	next_node.dir = dir_next;
	next_node.pos = next_pos;
	next_node.eh = next_eh;
	next_node.pre_fh = next_fh;

	if (trimesh.isOnBoundary(next_fh)) {// 贴着边界延伸，直接返回，后续直接将边界组织成流线即可
		bdy_faces_count++;
		if (bdy_faces_count > 5) {
			sline.end_on_fault = true;
			return;
		}
	}

	// check direction
	if (sline.points.size() > 2) {
		int n = sline.points.size();
		auto& v0 = sline.points[n - 3];
		auto& v1 = sline.points[n - 2];
		auto& v2 = sline.points[n - 1];
		Vec vec01 = (v1 - v0).normalized();
		Vec vec02 = (v2 - v0).normalized();
		double cosine = vec01.dot(vec02);
		if (cosine < 0) {
			std::cerr << "Error: Direction is reversed!\n";
			sline.end_on_fault = true;
			return;
		}
	}

	lengthen_initial_streamlines(next_node, sline, visited, bdy_faces_count);

}

void QuadMesher::merge_neighbor_edgehandle() {
	// 第一步：简化 eh 相同，sid 相同的流线
	{

		std::unordered_map<int, std::vector<int>> eh_to_streamline;
		int sl_cnt = streamlines.size();
		std::vector<bool> need_erase(sl_cnt, false);
		for (int i = 0; i < sl_cnt; ++i) {
			if (streamlines[i].end_on_boundary && streamlines[i].end_edgehandle != -1) {
				eh_to_streamline[streamlines[i].end_edgehandle].push_back(i);
				streamlines[i].init_length_ratio();
			}
		}
		std::vector<Streamline> new_streamlines;

		for (auto& eh_slid : eh_to_streamline) {
			EH eh(eh_slid.first);
			auto& sl_indices = eh_slid.second;// 以该边为终点的所有流线
			for (int i = 0; i < sl_indices.size(); ++i) {
				int slid_0 = sl_indices[i];
				if (need_erase[slid_0]) continue;
				std::vector<std::pair<double, int>> candiate_slid;
				candiate_slid.push_back({ streamlines[slid_0].length.back(), slid_0 });
				for (int j = i + 1; j < sl_indices.size(); ++j) {
					int slid_1 = sl_indices[j];
					if (need_erase[slid_1]) continue;
					if (streamlines[slid_1].beg_sid == streamlines[slid_0].beg_sid) {
						candiate_slid.push_back({ streamlines[slid_1].length.back(), slid_1 });
					}
				}
				//std::cout << "We iniaialize candiate_slid success. The size = " << candiate_slid.size() << "\n";
				if (candiate_slid.size() > 1) {

					std::sort(candiate_slid.begin(), candiate_slid.end(), [&](std::pair<double, int>& s0, std::pair<double, int>& s1) {
						return s0.first > s1.first;
						});
					for (auto& _slid : candiate_slid) {
						/*if (_slid.second >= sl_cnt || _slid.second < 0) {
							std::cout << "Streamlines: Out of range!\n";
							continue;
						}*/
						need_erase[_slid.second] = true;
					}

					int slid0 = candiate_slid[0].second;
					int slid1 = candiate_slid[1].second;
					auto& sl_0 = streamlines[slid0];
					auto& sl_1 = streamlines[slid1];

					/*std::cout << "Merging: beg_sid = " << sl_0.beg_sid << ", end_eh = " << sl_0.end_edgehandle <<
						", is same edge-handle? = " << (sl_0.end_edgehandle == sl_1.end_edgehandle) << std::endl;*/

					Streamline sl_new;
					get_streamline_itrp(sl_0, sl_1, sl_new);
					new_streamlines.push_back(sl_new);

					std::cout << "\tMerging part-one success.\n";

				}
			}

		}

		for (int i = 0; i < streamlines.size(); ++i) {
			if (!need_erase[i]) new_streamlines.push_back(streamlines[i]);
		}

		streamlines = new_streamlines;

	}
	

	// 第二步：简化 eh 邻近，sid 相同的流线
	{
		
		for (auto& sl : streamlines) sl.init_length_ratio();
		int sl_cnt = streamlines.size();
		std::vector<bool> need_erase(sl_cnt, false);
		std::vector<Streamline> streamlines_new;

		for (int i = 0; i < sl_cnt; ++i) {
			if (need_erase[i] || !streamlines[i].end_on_boundary) continue;
			std::vector<std::pair<double, int>> candiate_slid;
			std::vector<EH> candiate_eh;
			candiate_slid.push_back({ streamlines[i].length.back(), i });
			candiate_eh.push_back(EH(streamlines[i].end_edgehandle));
			for (int j = i + 1; j < sl_cnt; ++j) {
				if (need_erase[j] || !streamlines[j].end_on_boundary
					|| streamlines[j].beg_sid != streamlines[i].beg_sid) continue;
				EH cur_eh(streamlines[j].end_edgehandle);
				auto& cur_e = trimesh.edges(cur_eh);
				bool is_neighbor = false;
				for (auto& eh : candiate_eh) {
					auto& nei_e = trimesh.edges(eh);
					if (nei_e.vh1() == cur_e.vh1() || nei_e.vh2() == cur_e.vh1() ||
						nei_e.vh1() == cur_e.vh2() || nei_e.vh2() == cur_e.vh2()) {
						is_neighbor = true;
						break;
					}
					if (trimesh.isConnected(nei_e.vh1(), cur_e.vh1()) || trimesh.isConnected(nei_e.vh2(), cur_e.vh1()) ||
						trimesh.isConnected(nei_e.vh1(), cur_e.vh2()) || trimesh.isConnected(nei_e.vh2(), cur_e.vh2())) {
						is_neighbor = true;
						break;
					}
				}
				if (is_neighbor) {
					candiate_slid.push_back({ streamlines[j].length.back(), j });
					candiate_eh.push_back(cur_eh);
				}
			}
			if (candiate_slid.size() > 1) {

				std::sort(candiate_slid.begin(), candiate_slid.end(), [&](std::pair<double, int>& s0, std::pair<double, int>& s1) {
					return s0.first > s1.first;
					});
				for (auto& _slid : candiate_slid) {
					need_erase[_slid.second] = true;
				}

				int slid0 = candiate_slid[0].second;
				int slid1 = candiate_slid[1].second;
				auto& sl_0 = streamlines[slid0];
				auto& sl_1 = streamlines[slid1];

				std::cout << "Merging: beg_sid = " << sl_0.beg_sid << ", end_eh = " << sl_0.end_edgehandle <<
					", is same? = " << (sl_0.beg_sid == sl_1.beg_sid && sl_0.end_edgehandle == sl_1.end_edgehandle) << std::endl;

				Streamline sl_new;
				get_streamline_itrp(sl_0, sl_1, sl_new);
				streamlines_new.push_back(sl_new);

				std::cout << "\tMerging part-two success.\n";

			}
		}

		for (int i = 0; i < streamlines.size(); ++i) {
			if (!need_erase[i]) streamlines_new.push_back(streamlines[i]);
		}

		streamlines = streamlines_new;

	}

}