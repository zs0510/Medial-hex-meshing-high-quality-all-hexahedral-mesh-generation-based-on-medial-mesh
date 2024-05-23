#include "Denoising_BNF.h"

void Denoising_BNF::Execute() {

	
	initPositions();
	initTrianglesData();
	initSigmaCenter();
	for (int i = 0; i < it_normal; ++i) {
		BilateralNormalFiltering();
		//printf("BNF %d success\n", i);
	}
	for (int i = 0; i < it_vertex; ++i) {
		VertexUpdate();
		//printf("VexUpdate %d success\n", i);
	}
	//VertexUpdate_PreventFlip();

	for (auto& vp : mesh.allvertices()) {
		MeshKernel::iGameVertex& vex = mesh.vertices(vp.first);
		vex.setPosition(positions[vp.first][0], positions[vp.first][1], positions[vp.first][2]);
	}
	
	printf("denoising success\n");
}

void Denoising_BNF::BilateralNormalFiltering() {
	int fcnt = mesh.FaceSize();
	double double_SigmaCenter2 = 2.f * SigmaCenter * SigmaCenter;
	double double_SigmaNormal2 = 2.f * SigmaNormal * SigmaNormal;
	std::vector<Eigen::Vector3d> new_normals;
	new_normals.resize(fcnt, Eigen::Vector3d::Zero());

#pragma omp parallel for
	for (int fid = 0; fid < fcnt; ++fid) {
		FH fh(fid);
		auto adjF = mesh.Neighbor2Fh(fh);// share common vertex
		double sum_weight = 0.f;
		for (auto& fj : adjF) {
			double delta_center = (centers[fid] - centers[fj]).norm();
			double delta_normal = (normals[fid] - normals[fj]).norm();
			double weight_center = std::exp(-(delta_center * delta_center) / double_SigmaCenter2);
			double weight_normal = std::exp(-(delta_normal * delta_normal) / double_SigmaNormal2);
			double wj = areas[fj] * weight_center * weight_normal;
			new_normals[fid] += (wj * normals[fj]);
			sum_weight += wj;
		}
		assert(sum_weight > static_cast<double>(0));
		new_normals[fid] /= sum_weight;
		new_normals[fid].normalize();
	}

	normals = new_normals;

}

void Denoising_BNF::VertexUpdate() {

	int vcnt = mesh.vsize();

#pragma omp parallel for
	for (int vid = 0; vid < vcnt; ++vid) {
		VH vh(vid);
		Eigen::Vector3d move = Eigen::Vector3d::Zero();
		auto adj_faces = mesh.NeighborFh(vh);
		for (auto& fh : adj_faces) {
			double wj = normals[fh].transpose() * (centers[fh] - positions[vh]);
			move += (wj * normals[fh]);
		}
		move /= adj_faces.size();
		positions[vh] += move;
	}

}

void Denoising_BNF::VertexUpdate_PreventFlip() {
	// 传统顶点更新只考虑让顶点跟法向量正交，而没有考虑三角形的翻转问题，此函数可以避免面片翻转
	// Reference: Static/Dynamic Filtering for Mesh Geometry - Vertex Update section

	int vcnt = mesh.vsize(), fcnt = mesh.fsize();
	std::vector<Eigen::Matrix3d> face_P(fcnt);//P_i ∈ R3×3 are auxiliary variables representing the closest projection

	Eigen::Matrix3d mat_mean;
	mat_mean << 0.666, -0.333, -0.333,
		-0.333, 0.666, -0.333,
		-0.333, -0.333, 0.666;
	
	std::vector<Eigen::Triplet<double>> triplets;// 系统矩阵, 在迭代过程中始终固定

	for (int fid = 0; fid < fcnt; ++fid) {
		FH fh(fid);
		const auto& f_vhs = mesh.faces(fh).getVertexHandle();
		triplets.emplace_back(fid * 3, f_vhs[0], 0.666);
		triplets.emplace_back(fid * 3 + 1, f_vhs[0], -0.333);
		triplets.emplace_back(fid * 3 + 2, f_vhs[0], -0.333);
		triplets.emplace_back(fid * 3, f_vhs[1], -0.333);
		triplets.emplace_back(fid * 3 + 1, f_vhs[1], 0.666);
		triplets.emplace_back(fid * 3 + 2, f_vhs[1], -0.333);
		triplets.emplace_back(fid * 3, f_vhs[2], -0.333);
		triplets.emplace_back(fid * 3 + 1, f_vhs[2], -0.333);
		triplets.emplace_back(fid * 3 + 2, f_vhs[2], 0.666);
	}
	Eigen::SparseMatrix<double> K(fcnt * 3, vcnt);
	K.setFromTriplets(triplets.begin(), triplets.end());

	Eigen::SparseMatrix<double> KT = K.transpose();// vcnt x (3 * fcnt)
	Eigen::SparseMatrix<double> KTK = KT * K;

	triplets.clear();
	for (int vid = 0; vid < vcnt; ++vid) {
		triplets.emplace_back(vid, vid, weight);// w = 0.001
	}
	Eigen::SparseMatrix<double> wI(vcnt, vcnt);
	wI.setFromTriplets(triplets.begin(), triplets.end());

	Eigen::SparseMatrix<double> A = wI + KTK;

	// 方程组左边
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> MatricesCholesky(A);

	std::vector<Eigen::VectorXd> wV_0;// w * V_0
	wV_0.resize(3);
	wV_0[0].resize(vcnt);
	wV_0[1].resize(vcnt);
	wV_0[2].resize(vcnt);

	for (int vid = 0; vid < vcnt; ++vid) {
		wV_0[0](vid) = positions[vid].x() * weight;
		wV_0[1](vid) = positions[vid].y() * weight;
		wV_0[2](vid) = positions[vid].z() * weight;
	}

	for (int it = 0; it < 20; ++it) {

		time_t it_beg = clock();

		// First, we fix V and update P
#pragma omp parallel for
		for (int fid = 0; fid < fcnt; ++fid) {

			FH fh(fid);
			auto& face = mesh.faces(fh);

			const auto& f_vhs = face.getVertexHandle();
			auto& v0 = positions[f_vhs[0]];
			auto& v1 = positions[f_vhs[1]];
			auto& v2 = positions[f_vhs[2]];
			Eigen::Vector3d normal_curr = ((v1 - v0).cross(v2 - v0)).normalized();// n_i

			Eigen::Matrix3d Vf;
			Vf << v0.x(), v0.y(), v0.z(),
				v1.x(), v1.y(), v1.z(),
				v2.x(), v2.y(), v2.z();

			Eigen::Matrix3d P_tmp = mat_mean * Vf * (Eigen::Matrix3d::Identity() - normals[fid] * normals[fid].transpose());

			if (normal_curr.dot(normals[fid]) >= 0) {

				face_P[fid] = P_tmp;

			} else {

				Eigen::JacobiSVD<Eigen::Matrix3d> svd(P_tmp, Eigen::ComputeFullU | Eigen::ComputeFullV);
				Eigen::MatrixXd singular_values = svd.singularValues();
				//Eigen::MatrixXd left_singular_vectors = svd.matrixU();
				Eigen::MatrixXd right_singular_vectors = svd.matrixV();

				double sv_max = singular_values(0, 0);
				Eigen::Vector3d rv = right_singular_vectors.col(0);

				if (singular_values(1, 0) > sv_max) {
					sv_max = singular_values(1, 0);
					rv = right_singular_vectors.col(1);
				}

				if (singular_values(2, 0) > sv_max) {
					sv_max = singular_values(2, 0);
					rv = right_singular_vectors.col(2);
				}

				face_P[fid] = P_tmp * (rv * rv.transpose());

			}

		}

		std::vector<Eigen::VectorXd> P(3);
		P[0].resize(fcnt * 3);
		P[1].resize(fcnt * 3);
		P[2].resize(fcnt * 3);
		for (int fid = 0; fid < fcnt; ++fid) {
			FH fh(fid);
			const auto& f_vhs = mesh.faces(fh).getVertexHandle();
			auto& v0 = mesh.vertices(f_vhs[0]);
			auto& v1 = mesh.vertices(f_vhs[1]);
			auto& v2 = mesh.vertices(f_vhs[2]);
			P[0](fid * 3) = v0.x(); 
			P[0](fid * 3 + 1) = v1.x();
			P[0](fid * 3 + 2) = v2.x();
			P[1](fid * 3) = v0.y();
			P[1](fid * 3 + 1) = v1.y();
			P[1](fid * 3 + 2) = v2.y();
			P[2](fid * 3) = v0.z();
			P[2](fid * 3 + 1) = v1.z();
			P[2](fid * 3 + 2) = v2.z();
		}

		std::vector<Eigen::VectorXd> KTP = { KT * P[0], KT * P[1], KT * P[2] };

		std::vector<Eigen::VectorXd> b = { wV_0[0] + KTP[0], wV_0[1] + KTP[1], wV_0[2] + KTP[2] };

		std::vector<Eigen::VectorXd> x = {
			MatricesCholesky.solve(b[0]),
			MatricesCholesky.solve(b[1]),
			MatricesCholesky.solve(b[2])
		};

		// update positions
		for (int vid = 0; vid < vcnt; ++vid) {
			positions[vid] = Eigen::Vector3d(x[0](vid), x[1](vid), x[2](vid));
		}

		time_t it_end = clock();

		std::cout << "It" << it << " cost " << int(it_end - it_beg) << "ms.\n";

	}





}

void Denoising_BNF::initSigmaCenter() {
	SigmaCenter = 0.f;
	int cnt = 0;
	for (auto& fp : mesh.allfaces()) {
		auto fhs = mesh.NeighborFh(fp.first);// share common edge
		for (auto fh : fhs) {
			SigmaCenter += (centers[fp.first] - centers[fh]).norm();
			cnt++;
		}
	}
	SigmaCenter /= cnt;
	printf("SigmaCenter %.2f; SigmaNormal: %.2f\n", SigmaCenter, SigmaNormal);
}

void Denoising_BNF::initPositions() {
	int vcnt = mesh.VertexSize();
	positions.resize(vcnt);
	for (auto& vp : mesh.allvertices()) {
		positions[vp.first] = Eigen::Vector3d(vp.second.x(), vp.second.y(), vp.second.z());
	}
}

void Denoising_BNF::initTrianglesData() {
	int fcnt = mesh.FaceSize();
	centers.resize(fcnt, Eigen::Vector3d::Zero());
	normals.resize(fcnt);
	areas.resize(fcnt);
	for (auto& fp : mesh.allfaces()) {
		int f_idx = fp.first;
		auto v_indices = fp.second.getVertexHandle();
		assert(v_indices.size() == 3 && "isn't a triangle");
		Eigen::Vector3d vec23 = positions[v_indices[2]] - positions[v_indices[1]];
		Eigen::Vector3d vec21 = positions[v_indices[0]] - positions[v_indices[1]];

		for (auto& v_idx : v_indices) {
			centers[f_idx] += positions[v_idx];
		}
		centers[f_idx] /= v_indices.size();

		// 计算法向量 注意方向
		normals[f_idx] = (vec23.cross(vec21)).normalized();

		double len1 = vec23.norm();
		double len2 = vec21.norm();
		vec23.normalize();
		vec21.normalize();
		double cos_theta = vec23.dot(vec21);
		double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
		areas[f_idx] = len1 * len2 * sin_theta * 0.5;

	}
	printf("init triangles data success\n");
}