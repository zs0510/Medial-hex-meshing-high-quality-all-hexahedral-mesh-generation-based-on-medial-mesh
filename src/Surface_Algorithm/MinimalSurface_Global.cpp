#include "MinimalSurface_Global.h"

void MinimalSurface_Global::Execute() {
	initCotWeight();
	auto& vertices = mesh.allvertices();
	auto vcnt = mesh.VertexSize();
	std::vector<Tri> triplets;
	Eigen::VectorXd b_x(vcnt), b_y(vcnt), b_z(vcnt);
	bool boundary_flag = false;
	for (auto& vp : vertices) {
		auto& vh = vp.first;
		auto& v = mesh.vertices(vh);
		if (mesh.isOnBoundary(vh)) {// 边界固定
			triplets.emplace_back(vh, vh, 1.f);
			// or triplets.push_back(Tri(vh, vh, 1.f));
			b_x[vh] = v.x();
			b_y[vh] = v.y();
			b_z[vh] = v.z();
			boundary_flag = true;
		} else {
			double sum_weight = 0.f;
			auto adjE = mesh.NeighborEh(vh);
			for (auto& eh : adjE) {
				auto& e = mesh.edges(eh);
				auto adjvh = e.vh1() + e.vh2() - vh;
				double weight = e.getWeight();
				sum_weight += weight;
				triplets.emplace_back(vh, adjvh, weight);
			}
			triplets.emplace_back(vh, vh, -sum_weight);
			b_x[vh] = 0.f;
			b_y[vh] = 0.f;
			b_z[vh] = 0.f;
		}
	}

	if (!boundary_flag) {
		printf("this model has no boundary!!!\n");
		return;
	}

	// 稀疏矩阵
	Eigen::SparseMatrix<double> A(vcnt, vcnt);
	A.setFromTriplets(triplets.begin(), triplets.end());

	// 建立法方程
	Eigen::SparseMatrix<double> AT = A.transpose();
	Eigen::SparseMatrix<double> ATA = AT * A;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> MatricesCholesky(ATA);

	Eigen::VectorXd AT_bx = AT * b_x;
	Eigen::VectorXd AT_by = AT * b_y;
	Eigen::VectorXd AT_bz = AT * b_z;

	// 解法方程
	Eigen::VectorXd x = MatricesCholesky.solve(AT_bx);
	Eigen::VectorXd y = MatricesCholesky.solve(AT_by);
	Eigen::VectorXd z = MatricesCholesky.solve(AT_bz);

	for (auto& vp : vertices) {
		auto& vh = vp.first;
		auto& v = mesh.vertices(vh);
		v.setPosition(x[vh], y[vh], z[vh]);
	}

	printf("generate minimal surface success.\n");
}

void MinimalSurface_Global::initCotWeight() {
	for (auto ep : mesh.alledges()) {
		auto vh = ep.second.vh1();
		auto uh = ep.second.vh2();
		auto v = mesh.vertices(vh);
		auto u = mesh.vertices(uh);
		float cot_weight = 0.f;
		for (auto face : mesh.NeighborFh(ep.first)) {
			int idx = 0;
			auto wh = mesh.faces(face).vh(idx++);
			while (wh == vh || wh == uh) {
				assert(idx < 3);
				wh = mesh.faces(face).vh(idx++);
			}
			auto w = mesh.vertices(wh);
			Eigen::Vector3f v_pos(v.x(), v.y(), v.z());
			Eigen::Vector3f u_pos(u.x(), u.y(), u.z());
			Eigen::Vector3f w_pos(w.x(), w.y(), w.z());
			Eigen::Vector3f vec1 = (v_pos - w_pos).normalized();
			Eigen::Vector3f vec2 = (u_pos - w_pos).normalized();
			float cos_theta = vec1.dot(vec2);
			float cot_theta = cos_theta / std::sqrt(1 - cos_theta * cos_theta);
			cot_weight += std::max(EPSILON, cot_theta);
		}
		auto& e = mesh.edges(ep.first);
		e.setWeight(cot_weight * 0.5f);
	}
	printf("init cot weight success.\n");
}