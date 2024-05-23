#include "Q-MAT.h"

void Q_MAT::simplification(int target_vertices) {

	if (target_vertices >= mesh.vsize()) {
		std::cout << "Q-MAT: It is no need to simplification.\n";
		return;
	}
	int iter_num = mesh.vsize() - target_vertices;

	v_sqe = std::vector<V_SEQ>(mesh.vsize());
	v_collapse = std::vector<std::vector<Skl_CollapseNode*>>(mesh.vsize(), std::vector<Skl_CollapseNode*>{});

	calc_total_sqe();

	for (auto& ep : mesh.alledges()) {
		calc_collapse_cost(ep.first);
	}

	while (iter_num > 0 && !pri_que.empty()) {

		auto* node = pri_que.top();
		pri_que.pop();
		while (node->is_valid == false && !pri_que.empty()) {
			node = pri_que.top();
			pri_que.pop();
		}
		if (node->is_valid == false) break;
		
		auto& edge = mesh.edges(node->eh);
		auto vh_keep = edge.vh1();
		auto vh_remove = edge.vh2();
		assert(mesh.isValid(vh_keep) && mesh.isValid(vh_remove));
		
		//std::cout << "Q-MAT: Cost = " << node->cost << std::endl;
		auto& v_keep = mesh.vertices(vh_keep);
		auto& v_remove = mesh.vertices(vh_remove);
		
		// 为无效的收缩结点做标记
		for (auto* _node : v_collapse[vh_keep]) {
			_node->is_valid = false;
		}
		for (auto* _node : v_collapse[vh_remove]) {
			_node->is_valid = false;
		}
		v_collapse[vh_keep].clear();
		v_collapse[vh_remove].clear();

		// 更新结点几何
		v_keep.setPosition(node->sphere.x(), node->sphere.y(), node->sphere.z());
		mesh.radius[vh_keep] = node->sphere.w();
		v_sqe[vh_keep].A += v_sqe[vh_remove].A;
		v_sqe[vh_keep].b += v_sqe[vh_remove].b;
		v_sqe[vh_keep].c += v_sqe[vh_remove].c;

		std::vector<std::vector<MeshKernel::iGameVertexHandle>> new_triangles;
		for (auto& eh : mesh.NeighborEh(vh_remove)) {
			auto& adje = mesh.edges(eh);
			VH vh_other(adje.vh1() + adje.vh2() - vh_remove);
			if (vh_other == vh_keep || mesh.isConnected(vh_other, vh_keep)) continue;
			mesh.AddEdge(vh_keep, vh_other);
		}
		for (auto& fh : mesh.NeighborFh(vh_remove)) {
			auto vex = (mesh.faces(fh)).getVertexHandle();
			assert(vex.size() > 2);
			if (std::find(vex.begin(), vex.end(), vh_keep) != vex.end()) continue;
			for (int i = 0; i < vex.size(); ++i) {
				if (vex[i] == vh_remove) {
					vex[i] = vh_keep;
					new_triangles.push_back(vex);
					break;
				}
			}
		}
		mesh.DeleteVertex(vh_remove);

		// 添加三角形
		for (auto& tri : new_triangles) {
			mesh.AddFace(tri);
		}
		
		// 重新计算收缩代价
		for (auto& eh : mesh.NeighborEh(vh_keep)) {
			calc_collapse_cost(eh);
		}
		iter_num--;
	}
	mesh.updateAllHandles();
	//mesh.hole_to_triangle();
	std::cout << "Q-MAT: Simplification success. #V = " << mesh.vsize() << ", #E = " << mesh.esize() << ", #F = " << mesh.fsize() << std::endl;

}

void Q_MAT::calc_total_sqe() {
	
	int f_cnt = 0;
	for (auto& fp : mesh.allfaces()) {
		calc_face_sqe(fp.first);
		f_cnt++;
	}

	int e_cnt = 0;
	for (auto& ep : mesh.alledges()) {
		if (mesh.NeighborFh(ep.first).size() < 2) {
			calc_edge_sqe(ep.first);
			e_cnt++;
		}
	}

	int v_cnt = 0;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.NeighborEh(vp.first).size() == 1) {
			calc_vertex_sqe(vp.first);
			v_cnt++;
		}
	}

	std::cout << "Q-MAT: Initial total SQE success. Total calculate #F = " << f_cnt << ", #E = " << e_cnt << ", #V = " << v_cnt << std::endl;

}

void Q_MAT::calc_face_sqe(FH fh) {
	const auto& face = mesh.faces(fh);
	const auto& vhs = face.getVertexHandle();
	auto& v0 = mesh.vertices(vhs[0]);
	auto& v1 = mesh.vertices(vhs[1]);
	auto& v2 = mesh.vertices(vhs[2]);

	Eigen::Vector3d _n0, _n1;
	Eigen::Vector3d c0(v0.x(), v0.y(), v0.z());
	Eigen::Vector3d c1(v1.x(), v1.y(), v1.z());
	Eigen::Vector3d c2(v2.x(), v2.y(), v2.z());
	double r0 = mesh.radius[vhs[0]], r1 = mesh.radius[vhs[1]], r2 = mesh.radius[vhs[2]];

	if (!get_triangle_from_three_spheres(c0,r0,c1,r1,c2,r2,_n0,_n1)) {
		return;
	}
	Eigen::Vector4d n0(_n0.x(), _n0.y(), _n0.z(), 1);
	Eigen::Vector4d n1(_n1.x(), _n1.y(), _n1.z(), 1);
	Eigen::Matrix4d A = n0 * n0.transpose() + n1 * n1.transpose();
	for (auto& vh : vhs) {
		auto& v = mesh.vertices(vh);
		Eigen::Vector4d m(v.x(), v.y(), v.z(), 1);
		Eigen::Vector4d b = -2 * A * m;
		double c = m.transpose() * A * m;
		v_sqe[vh].A += A;
		v_sqe[vh].b += b;
		v_sqe[vh].c += c;
	}

}

void Q_MAT::calc_edge_sqe(EH eh) {
	const auto& edge = mesh.edges(eh);
	double factor = soft_k * get_edge_stability(edge);
	const auto& e_fhs = mesh.NeighborFh(eh);
	int degree = e_fhs.size();
	VH vh1 = edge.vh1(), vh2 = edge.vh2();
	auto& v1 = mesh.vertices(vh1);
	auto& v2 = mesh.vertices(vh2);
	if (degree == 0) {

		Vec edge_dir = v2 - v1;// 向外
		edge_dir.normalize();
		// 求法向量为 dir 的平面, 平面方程 ax + by + cz + d = 0，edge_dir = (a, b, c)
		if (std::abs(edge_dir.z()) <= 1E-8F) edge_dir.z() = 1E-8F;
		double from_d = -(edge_dir.dot(v1));
		double point_z = (-from_d) / edge_dir.z();// 平面上的某个点(0, 0, point_z)
		Vec tp_vec0(-v1.x(), -v1.y(), point_z - v1.z());// 切平面上的某个向量
		Vec tp_vec1 = tp_vec0.cross(edge_dir);
		tp_vec0.normalize();
		tp_vec1.normalize();
		std::vector<std::vector<int>> move_dirs = {// 在切平面上移动的方向
			{-1, -1}, {-1, 1}, {1, -1}, {1, 1}
		};
		std::vector<Vec> from_points, to_points;// 在端点切平面上的点，八个顶点形成一个棱柱
		for (auto& move_dir : move_dirs) {
			Vec move_vec0 = tp_vec0 * move_dir[0];
			Vec move_vec1 = tp_vec1 * move_dir[1];
			from_points.push_back(v1 + move_vec0 * mesh.radius[vh1] + move_vec1 * mesh.radius[vh1]);
			to_points.push_back(v2 + move_vec0 * mesh.radius[vh2] + move_vec1 * mesh.radius[vh2]);
		}
		for (int i = 0; i < 4; ++i) {// 外移，形成大棱柱
			from_points[i] -= (edge_dir * mesh.radius[vh1]);
			to_points[i] += (edge_dir * mesh.radius[vh2]);
		}
		for (int i = 0; i < 4; ++i) {// 施加棱柱的四个表面约束
			Vec plane_p0 = from_points[i];
			Vec plane_p1 = to_points[i];
			Vec plane_p2 = to_points[(i + 1) % 4];
			Vec plane_normal = (plane_p1 - plane_p0).cross(plane_p2 - plane_p0);
			plane_normal.normalize();
			Eigen::Vector4d n(plane_normal.x(), plane_normal.y(), plane_normal.z(), 1);
			Eigen::Matrix4d A = n * n.transpose() * factor;
			Eigen::Vector4d m(v1.x(), v1.y(), v1.z(), 1);
			Eigen::Vector4d b = -2 * A * m;
			double c = m.transpose() * A * m;
			v_sqe[vh1].A += A;
			v_sqe[vh1].b += b;
			v_sqe[vh1].c += c;
			m = Eigen::Vector4d(v2.x(), v2.y(), v2.z(), 1);
			b = -2 * A * m;
			c = m.transpose() * A * m;
			v_sqe[vh2].A += A;
			v_sqe[vh2].b += b;
			v_sqe[vh2].c += c;
		}
		
	} else if (degree == 1) {
		// 需增加边界边软约束
		auto fh = *(e_fhs.begin());
		auto& face = mesh.faces(fh);
		const auto& f_vhs = face.getVertexHandle();
		Vec vec_e = v1 - v2;
		vec_e.normalize();
		Vec nor_f = (mesh.vertices(f_vhs[1]) - mesh.vertices(f_vhs[0])).cross(
			mesh.vertices(f_vhs[2]) - mesh.vertices(f_vhs[0]));
		nor_f.normalize();
		Vec move_dir = vec_e.cross(nor_f);
		move_dir.normalize();
		Vec edge_mid = (v1 + v2) * 0.5f;
		Vec centroid = mesh.getFaceCenter(fh);
		Vec toward_c = centroid - edge_mid;
		toward_c.normalize();
		if (move_dir.dot(toward_c) > 0.f) move_dir *= -1;
		Vec p0 = v1 + nor_f * mesh.radius[vh1];
		Vec p1 = v2 + nor_f * mesh.radius[vh2];
		Vec p2 = v2 - nor_f * mesh.radius[vh2];
		p0 += (move_dir * mesh.radius[vh1]);
		p1 += (move_dir * mesh.radius[vh2]);
		p2 += (move_dir * mesh.radius[vh2]);
		Vec plane_normal = (p1 - p0).cross(p2 - p0);
		plane_normal.normalize();
		Eigen::Vector4d n(plane_normal.x(), plane_normal.y(), plane_normal.z(), 1);
		Eigen::Matrix4d A = n * n.transpose() * factor;

		Eigen::Vector4d m(v1.x(), v1.y(), v1.z(), 1);
		Eigen::Vector4d b = -2 * A * m;
		double c = m.transpose() * A * m;
		v_sqe[vh1].A += A;
		v_sqe[vh1].b += b;
		v_sqe[vh1].c += c;
		m = Eigen::Vector4d(v2.x(), v2.y(), v2.z(), 1);
		b = -2 * A * m;
		c = m.transpose() * A * m;
		v_sqe[vh2].A += A;
		v_sqe[vh2].b += b;
		v_sqe[vh2].c += c;
	}
}

void Q_MAT::calc_vertex_sqe(VH vh) {

	EH eh(*(mesh.NeighborEh(vh).begin()));
	auto& edge = mesh.edges(eh);
	VH vh_o(edge.vh1() + edge.vh2() - vh);
	auto& v = mesh.vertices(vh);
	auto& v_o = mesh.vertices(vh_o);
	Vec edge_dir = (v - v_o).normalize();
	double factor = soft_k * get_edge_stability(edge);
	
	Eigen::Vector4d n(edge_dir.x(), edge_dir.y(), edge_dir.z(), 1);
	Eigen::Matrix4d A = n * n.transpose() * factor;
	Eigen::Vector4d m(v.x(), v.y(), v.z(), 1);
	Eigen::Vector4d b = -2 * A * m;
	double c = m.transpose() * A * m;
	v_sqe[vh].A += A;
	v_sqe[vh].b += b;
	v_sqe[vh].c += c;

}


void Q_MAT::calc_collapse_cost(EH eh) {

	auto& edge = mesh.edges(eh);
	auto& sqe0 = v_sqe[edge.vh1()];
	auto& sqe1 = v_sqe[edge.vh2()];
	Eigen::Matrix4d A = sqe0.A + sqe1.A;
	Eigen::Vector4d b = sqe0.b + sqe1.b;
	double c = sqe0.c + sqe1.c;// not be affected by sphere_x
	double cost = 99999999999;
	Eigen::Vector4d X;
	if (A.determinant()) {
		X = (-0.5) * (A.inverse() * b);
		cost = X.transpose() * A * X + b.dot(X);
	} else {
		for (double ratio = 0; ratio < 1.01; ratio += 0.5) {
			Vec _pos = mesh.vertices(edge.vh1()) * ratio + mesh.vertices(edge.vh2()) * (1 - ratio);
			Eigen::Vector3d pos(_pos.x(), _pos.y(), _pos.z());
			double r = mesh.radius[edge.vh1()] * ratio + mesh.radius[edge.vh2()] * (1 - ratio);
			Eigen::Vector4d m = Eigen::Vector4d(pos.x(), pos.y(), pos.z(), r);
			double tmp_cost = m.transpose() * A * m + b.dot(m);
			if (tmp_cost < cost) {
				cost = tmp_cost;
				X = m;
			}
		}
	}
	Skl_CollapseNode* node = new Skl_CollapseNode();
	double stability_weight = get_edge_stability(edge);
	stability_weight *= stability_weight;
	//std::cout << "vh1.r = " << mesh.radius[edge.vh1()] << ", vh2.r = " << mesh.radius[edge.vh2()] << "; Solved r = " << X.w() << ", ";
	//X.w() = get_radius_by_interpolation(Vec(X.x(), X.y(), X.z()), edge.vh1(), edge.vh2());
	//std::cout << "Interpolation r = " << X.w() << std::endl;
	node->cost = (cost + c + hard_k) * stability_weight;
	if(std::isnan(node->cost)) return;
	//std::cout << "Geometric weight: " << cost + c << std::endl;
	//std::cout << "Old radius: " << edge->v0->radius << ", " << edge->v1->radius << "; new eadius: " << X.w() << std::endl;
	node->eh = eh;
	node->sphere = X;
	pri_que.push(node);
	v_collapse[edge.vh1()].push_back(node);
	v_collapse[edge.vh2()].push_back(node);

}

double Q_MAT::get_edge_stability(const MeshKernel::iGameEdge& edge) {
	auto& vh1 = edge.vh1();
	auto& vh2 = edge.vh2();
	auto& v1 = mesh.vertices(vh1);
	auto& v2 = mesh.vertices(vh2);
	double diff_radius = std::abs(mesh.radius[vh1] - mesh.radius[vh2]);
	double dis_center = (v1 - v2).norm();
	double dis_sphere = std::max((double)0, dis_center - diff_radius);
	double stability_ratio = dis_sphere / dis_center;
	//std::cout << "Stability: " << stability_ratio << std::endl;
	return stability_ratio;
}

double Q_MAT::get_radius_by_interpolation(const MeshKernel::iGameVertex& v_n, VH vh0, VH vh1) {
	auto& v0 = mesh.vertices(vh0);
	auto& v1 = mesh.vertices(vh1);
	double radius = 0.f;
	double dis_0 = (v_n - v0).norm();
	double dis_1 = (v_n - v1).norm();
	double weight_1 = dis_0 / (dis_0 + dis_1);
	radius = weight_1 * mesh.radius[vh1] + (1 - weight_1) * mesh.radius[vh0];

	return radius;

}

bool Q_MAT::get_triangle_from_three_spheres(Eigen::Vector3d c0, double r0, Eigen::Vector3d c1, double r1, Eigen::Vector3d c2, double r2,
	Eigen::Vector3d& n0, Eigen::Vector3d& n1) {

	Eigen::Vector3d c0c1(c1 - c0), c0c2(c2 - c0), c1c2(c2 - c1);
	double c0c1len(c0c1.norm()), c0c2len(c0c2.norm()), c1c2len(c1c2.norm());
	double dr0r1(std::fabs(r0 - r1)), dr0r2(std::fabs(r0 - r2)), dr1r2(std::fabs(r1 - r2));

	// some spheres are concentric and there are no triangles.
	if ((c0c1len < 1e-8) || (c0c2len < 1e-8) || (c1c2len < 1e-8))
		return false;

	//// some spheres are included in some other spheres 
	//if ((c0c1len - abs(r0 - r1) < 1e-8) || (c0c2len - abs(r0 - r2) < 1e-8) || (c1c2len - abs(r1 - r2) < 1e-8))
	//	return false;

	Eigen::Vector3d norm;
	norm = c0c1.cross(c0c2);
	norm.normalize();

	// equal-radius spheres
	if ((dr0r1 < 1e-8) && (dr0r2 < 1e-8) && (dr1r2 < 1e-8)) {
		Eigen::Vector3d tri0_0 = c0 + norm * r0;
		Eigen::Vector3d tri0_1 = c1 + norm * r1;
		Eigen::Vector3d tri0_2 = c2 + norm * r2;
		n0 = (tri0_1 - tri0_0).cross(tri0_2 - tri0_0);
		n0.normalize();

		Eigen::Vector3d tri1_0 = c0 - norm * r0;
		Eigen::Vector3d tri1_1 = c1 - norm * r1;
		Eigen::Vector3d tri1_2 = c2 - norm * r2;
		n1 = (tri1_1 - tri1_0).cross(tri1_2 - tri1_0);
		n1.normalize();
		n1 *= -1;
		return true;
	} else {
		// two points on the tangent plane
		Eigen::Vector3d apex0, apex1;

		// two spheres are equal-radius
		if (dr0r1 < 1e-8)
		{
			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
		} 		else if (dr0r2 < 1e-8)
		{
			apex0 = (r1 * c0 - r0 * c1) / (r1 - r0);
			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
		} 		else if (dr1r2 < 1e-8)
		{
			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
			apex1 = (r0 * c1 - r1 * c0) / (r0 - r1);
		} 		else
		{
			apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
			apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
		}

		double distc0;
		Eigen::Vector3d fp;
		distance_to_line(c0, apex0, apex1, distc0, fp);
		double sangle = r0 / distc0;
		if (std::fabs(sangle) > 1.)
			return false;
		double cangle = std::sqrt(1. - r0 * r0 / distc0 / distc0);
		Eigen::Vector3d norfpc0(c0 - fp);
		norfpc0.normalize();
		Eigen::Vector3d newnorm[2];
		newnorm[0] = norm * cangle - norfpc0 * sangle;
		newnorm[1] = -norm * cangle - norfpc0 * sangle;
		Eigen::Vector3d tri0_0 = c0 + r0 * newnorm[0];
		Eigen::Vector3d tri0_1 = c1 + r1 * newnorm[0];
		Eigen::Vector3d tri0_2 = c2 + r2 * newnorm[0];
		n0 = (tri0_1 - tri0_0).cross(tri0_2 - tri0_0);
		n0.normalize();

		Eigen::Vector3d tri1_0 = c0 + r0 * newnorm[1];
		Eigen::Vector3d tri1_1 = c1 + r1 * newnorm[1];
		Eigen::Vector3d tri1_2 = c2 + r2 * newnorm[1];
		n1 = (tri1_1 - tri1_0).cross(tri1_2 - tri1_0);
		n1.normalize();
		n1 *= -1;
	}

	return true;

}

bool Q_MAT::distance_to_line(const Eigen::Vector3d& p, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, double& dist, Eigen::Vector3d& fp) {
	Eigen::Vector3d v0v1(v1 - v0), pv0(v0 - p), pv1(v1 - p);
	double area = fabs(v0v1.cross(pv0).norm());
	if (v0v1.norm() > 1e-12) {
		dist = area / v0v1.norm();
		double t = (pv0.dot(pv0) - pv0.dot(pv1)) / (pv0.dot(pv0) + pv1.dot(pv1) - 2 * pv0.dot(pv1));
		fp = (1 - t) * v0 + t * v1;
		return true;
	}
	return false;
}