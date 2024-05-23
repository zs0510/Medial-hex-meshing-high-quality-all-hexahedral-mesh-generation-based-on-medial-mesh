#include "Parameterization_ARAP.h"

void Parameterization_ARAP::execute(int iter_times) {

	double area_origin = MeshMath::get_area_surface_mesh(mesh);

	// 1. 得到最初的参数化结果
	get_initial_parameterization();
	calc_local_coord();

	// 2. 初始化所有的 cot 权
	calc_all_cot();

	// 3. 建立系数矩阵
	build_coefficient_matrix();

	// 4. local-global iteration
	Lts.resize(fcnt);

	for (int it = 0; it < iter_times; ++it) {

		// local calc Lt
		calc_optimal_ration();

		// global solve
		solve_linear_equations();

	}

	for (int i = 0; i < vcnt; ++i) {
		VH vh(i);
		auto& v = mesh.vertices(vh);
		v = Vex(uv(i, 0), uv(i, 1), 0);
	}

	double area_current = MeshMath::get_area_surface_mesh(mesh);

	std::cout << "[ARAP Parameterization]: Area of surface, before = " << area_origin << ", after = " << area_current << std::endl;

}

void Parameterization_ARAP::calc_optimal_ration() {

#pragma omp parallel for
	for (int i = 0; i < fcnt; ++i) {
		FH fh(i);
		auto& face = mesh.faces(fh);
		const auto& vhs = face.getVertexHandle();

		int i0 = vhs[0], i1 = vhs[1], i2 = vhs[2];
		Eigen::Matrix2d P, S, J;

		P << uv(i1, 0) - uv(i0, 0), uv(i2, 0) - uv(i0, 0), uv(i1, 1) - uv(i0, 1), uv(i2, 1) - uv(i0, 1);
		S << local_coord(i, 2) - local_coord(i, 0), local_coord(i, 4) - local_coord(i, 0),
			local_coord(i, 3) - local_coord(i, 1), local_coord(i, 5) - local_coord(i, 1);

		J = P * S.inverse();

		Eigen::JacobiSVD<Eigen::Matrix2d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Eigen::Matrix2d U = svd.matrixU();
		Eigen::Matrix2d V = svd.matrixV();

		Eigen::Matrix2d R = U * V.transpose();

		if (R.determinant() < 0) {
			U(0, 1) = -U(0, 1);
			U(1, 1) = -U(1, 1);
			R = U * V.transpose();
		}
		Lts[i] = R;
	}

	//std::cout << "Calc optimal ration success.\n";

}

void Parameterization_ARAP::solve_linear_equations() {

	Eigen::VectorXd bu, bv;
	bu.setZero(vcnt);
	bv.setZero(vcnt);

	for (int i = 0; i < fcnt; ++i) {

		FH fh(i);
		auto& face = mesh.faces(fh);
		const auto& vhs = face.getVertexHandle();
		int i0 = vhs[0], i1 = vhs[1], i2 = vhs[2];

		Eigen::Vector2d e0, e1, e2;
		e0 << local_coord(i, 2), local_coord(i, 3);
		e1 << local_coord(i, 4) - local_coord(i, 2), local_coord(i, 5) - local_coord(i, 3);
		e2 << -local_coord(i, 4), -local_coord(i, 5);
		Eigen::Vector2d b0 = faces_cot[i][2] * Lts[i] * e0;
		bu[i0] -= b0[0];
		bv[i0] -= b0[1];
		bu[i1] += b0[0];
		bv[i1] += b0[1];

		Eigen::Vector2d b1 = faces_cot[i][0] * Lts[i] * e1;
		bu[i1] -= b1[0];
		bv[i1] -= b1[1];
		bu[i2] += b1[0];
		bv[i2] += b1[1];

		Eigen::Vector2d b2 = faces_cot[i][1] * Lts[i] * e2;
		bu[i2] -= b2[0];
		bv[i2] -= b2[1];
		bu[i0] += b2[0];
		bv[i0] += b2[1];

	}

	uv.col(0) = solver.solve(bu);
	uv.col(1) = solver.solve(bv);

	//std::cout << "Solve linear equations success.\n";

}

void Parameterization_ARAP::get_initial_parameterization() {

	MeshKernel::SurfaceMesh mesh_tmp = mesh;
	LSCM_Parameterization app_lscm(mesh_tmp);
	app_lscm.Execute();
	
	uv.resize(vcnt, 2);

#pragma omp parallel for
	for (int i = 0; i < vcnt; ++i) {
		VH vh(i);
		auto& v = mesh_tmp.vertices(vh);
		uv.row(i) = Eigen::Vector2d(v.x(), v.y());
	}

	std::cout << "LSCM run success.\n";

}

void Parameterization_ARAP::calc_local_coord() {

	local_coord.resize(fcnt, 6);

	mesh.genAllFacesNormal();

#pragma omp parallel for
	for (int i = 0; i < fcnt; ++i) {
		FH fh(i);
		auto& face = mesh.faces(fh);
		auto n = face.getNormal();
		const auto& f_vhs = face.getVertexHandle();
		auto& v0 = mesh.vertices(f_vhs[0]);
		auto& v1 = mesh.vertices(f_vhs[1]);
		auto& v2 = mesh.vertices(f_vhs[2]);
		Vec vec01 = v1 - v0;
		Vec dir01 = vec01.normalized();
		Vec vec_ = n.cross(dir01);
		Vec vec02 = v2 - v0;

		local_coord.row(i) << 0, 0, vec01.norm(), 0, vec02.dot(dir01), vec02.dot(vec_);

	}

	//std::cout << "Calc local coord success.\n";

}

void Parameterization_ARAP::calc_all_cot() {

	faces_cot.resize(fcnt, vector<double>(3));

//#pragma omp parallel for
	for (int i = 0; i < fcnt; ++i) {
		FH fh(i);
		auto& face = mesh.faces(fh);
		const auto& vhs = face.getVertexHandle();
		for (int j = 0; j < 3; ++j) {
			VH vh0(vhs[j]);// 以 vh0 作为对立顶点
			VH vh1 = vhs[(j + 1) % 3], vh2 = vhs[(j + 2) % 3];
			Vex& v0 = mesh.vertices(vh0);
			Vex& v1 = mesh.vertices(vh1);
			Vex& v2 = mesh.vertices(vh2);
			Vec vec01 = (v1 - v0).normalized();
			Vec vec02 = (v2 - v0).normalized();
			double cosine = vec02.dot(vec01);
			double sine = std::sqrt(1 - cosine * cosine);
			double cot = cosine / sine;
			faces_cot[fh][j] = std::max(1e-6, cot);// 边 vh1-vh2 的权重

			triplets.emplace_back(vh1, vh1, faces_cot[fh][j]);
			triplets.emplace_back(vh2, vh2, faces_cot[fh][j]);
			triplets.emplace_back(vh1, vh2, -faces_cot[fh][j]);
			triplets.emplace_back(vh2, vh1, -faces_cot[fh][j]);

		}
	}

	//std::cout << "Calc cots success.\n";

}

void Parameterization_ARAP::build_coefficient_matrix() {

	Eigen::SparseMatrix<double> A;
	A.resize(vcnt, vcnt);
	A.setFromTriplets(triplets.begin(), triplets.end());
	solver.compute(A);

	//std::cout << "Build cofficient matrix success.\n";

}