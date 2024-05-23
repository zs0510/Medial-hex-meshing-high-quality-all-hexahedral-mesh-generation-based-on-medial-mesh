#include "Regularization.h"

using namespace std;
using namespace MeshKernel;

void HexMesh_Regularization::regularization(uint32_t max_iter) {
	while (1) {
		auto vps = mesh.allvertices();
		bool find_singular_vertex = false;
		for (auto& vp : vps) {
			auto adjehs = mesh.NeighborEh(vp.first);
			if (adjehs.empty()) {
				mesh.DeleteVertex(vp.first);
				find_singular_vertex = true;
			}
		}
		if (find_singular_vertex) {
			printf("find singular vertex!!!\n");
		} else {
			break;
		}
	}
	

	auto& allcells = mesh.allcells();
	auto& alledges = mesh.alledges();
	vector<iGameCellHandle> chs;
	vector<iGameEdgeHandle> ehs;
	chs.reserve(allcells.size());
	ehs.reserve(alledges.size());
	for (auto& cp : allcells) 
		chs.push_back(cp.first);
	for (auto& ep : alledges)
		ehs.push_back(ep.first);

	for (auto& ch : chs) {
		auto& c = mesh.cells(ch);
		auto fhs = c.faces_size();
		// support pure hexahedron model only
		assert(c.faces_size() == 6 && "input isn't a pure hexahedron model!!!");
	}

	printf("check success!\n");


	CalculateCellQuality(chs, 0);
	std::cout << "Esteem cell success" << std::endl;
	time_t time1 = clock();
	for (uint32_t it = 0; it < max_iter; ++it) {

		Preprocessing(ehs, chs);// preprocessing
		std::cout << "Preprocessing success" << std::endl;
		// Step 1: Local regularization of hexahedral meshes
		for (int i = 0; i < chs.size(); ++i) {
			LocalRegularization(chs[i]);
		}
		printf("local regularization success\n");

		// compute weight of edges
		for (auto& ch : chs) {
			auto vhs = (mesh.cells(ch)).getVertexHandle();
			assert(vhs.size() == 8);
			unordered_map<int, iGameVertexHandle> id;
			id[0] = vhs[3], id[1] = vhs[2], id[2] = vhs[1], id[3] = vhs[0];
			id[4] = vhs[7], id[5] = vhs[6], id[6] = vhs[5], id[7] = vhs[4];
			for (auto& idx : tetFaces) {
				ComputeTetrahedronLaplacianWeight({ id[idx[0]], id[idx[1]],id[idx[2]], id[idx[3]] });
			}
		}
		printf("compute tetrahedron Laplacian success\n");

		if (check_weight_flag) {
			for (auto& cp : weight_cell) {
				printf("weight of cell %d: %.8f\n", cp.first, cp.second);
				assert(cp.second >= 0.f);
			}
			for (auto& wp : weight_edge) {
				printf("weight of edge %d: %.8f\n", wp.first, wp.second);
				assert(wp.second >= 0.f);
			}
		}

		// compute ideal node position
		for (auto& pa : transforms) {
			auto& trans = pa.second;
			Eigen::Vector3d sum_pos(0, 0, 0);
			double sum_weight = 0;
			for (auto& t : trans) {
				double weight = weight_cell[t.ch];
				sum_pos += weight * t.pos;
				sum_weight += weight;
			}
			sum_pos /= sum_weight;
			ideal_position[pa.first] = sum_pos;
			/*auto& v = mesh.vertices(iGameVertexHandle(pa.first));
			printf("vh: %d, degree: %d, initial (%.3f, %.3f, %.3f), ideal (%.3f, %.3f, %.3f)\n", 
				pa.first, trans.size(), v.x(), v.y(), v.z(), sum_pos.x(), sum_pos.y(), sum_pos.z());
			v.setX(sum_pos.x());
			v.setY(sum_pos.y());
			v.setZ(sum_pos.z());*/
		}

		// Step 2: Global optimization of hexahedral meshes
		if (use_exterior_only) {
			GlobalOptimization_Exterior();
		} else {
			GlobalOptimization();
		}
		
		//printf("\nsmoothing cost %d ms\n", int(clock() - time1));
		CalculateCellQuality(chs, 1);
		//printf("iter%d smoothing success. mesh quality min: %.2f, avg: %.2f\n", it, out_min_scaled_jacobian, out_avg_scaled_jacobian);
		//if (out_min_scaled_jacobian > 0.8f) break;
	}

	// CalculateCellQuality(chs, 1);

	if (use_exterior_only) printf("hexMesh regularization success(faster).\n");
	else printf("hexMesh regularization success.\n");
	
	

	printf("mesh quality before min: %.2f, avg: %.2f; after min: %.2f, avg: %.2f\n", 
		in_min_scaled_jacobian, in_avg_scaled_jacobian, out_min_scaled_jacobian, out_avg_scaled_jacobian);

}

void HexMesh_Regularization::Preprocessing(vector<iGameEdgeHandle>& ehs, vector<iGameCellHandle>& chs) {

	// clear all old data
	transforms.clear();
	// other datas just need update

	// init weight_edge
	for (int i = 0; i < ehs.size(); ++i) {
		weight_edge[ehs[i]] = 0;
	}

	// compute cell weight
	for (int i = 0; i < chs.size(); ++i) {
		auto& ch = chs[i];
		auto& c = mesh.cells(ch);
		double min_det = 2.f;// the value of Scaled Jacobian is the	determinant of the matrix Jacobian
		unordered_map<int, vector<int>> neighbor;
		const auto& _ehs = c.getEdgeHandle();
		for (auto& _eh : _ehs) {
			auto& _e = mesh.edges(_eh);
			auto& vh1 = _e.vh1();
			auto& vh2 = _e.vh2();
			neighbor[vh1].push_back(vh2);
			neighbor[vh2].push_back(vh1);
		}
		const auto& _vhs = c.getVertexHandle();
		for (auto& _vh : _vhs) {
			auto& v = mesh.vertices(_vh);
			auto adjvh = neighbor[_vh];
			assert(adjvh.size() == 3);
			auto v1 = mesh.vertices(iGameVertexHandle(adjvh[0]));
			auto v2 = mesh.vertices(iGameVertexHandle(adjvh[1]));
			auto v3 = mesh.vertices(iGameVertexHandle(adjvh[2]));
			auto center = (v1 + v2 + v3) / 3;
			auto vec = (center - v).normalized();

			// ajust oridering
			auto vec12 = v2 - v1, vec13 = v3 - v1;
			auto normal = (vec12 % vec13).normalized();
			double cosine = vec * normal;
			if (cosine < 0) swap(v2, v3);

			vector<Eigen::Vector3d> ev(3);
			auto tmp = v1 - v;
			ev[0] = (Eigen::Vector3d(tmp.x(), tmp.y(), tmp.z())).normalized();
			tmp = v2 - v;
			ev[1] = (Eigen::Vector3d(tmp.x(), tmp.y(), tmp.z())).normalized();
			tmp = v3 - v;
			ev[2] = (Eigen::Vector3d(tmp.x(), tmp.y(), tmp.z())).normalized();
			
			Eigen::Matrix3d J;
			J << ev[0][0], ev[1][0], ev[2][0],
				ev[0][1], ev[1][1], ev[2][1],
				ev[0][2], ev[1][2], ev[2][2];
			//std::cout << J << endl;
			/*
			For an element, the full range of the Scaled Jacobian value is from −1 to + 1. 
			And fscaled_jacobian = 1 if the element is an ideal element,
			fscaled_jacobian = −1 if the element is a worst distorted element.
			Where the value is less than 0, the element is inverted!!!
			*/
			min_det = min(min_det, J.determinant());
		}
		//printf("weight of cell %d: %.4f\n", ch, min_det);
		weight_cell[ch] = max(1 - min_det, KEEPOPSITIVE);
	}

	if (geometric_constraint_flag) {
		// compute vertex normal
		for (auto& vp : mesh.allvertices()) {
			if (!mesh.isOnBoundary(vp.first)) continue;
			iGameVertex sum_normal(0.f, 0.f, 0.f);
			auto v = vp.second;
			const auto& adjfhs = mesh.NeighborFh(vp.first);
			double sum_area = 0.f;
			for (auto fh : adjfhs) {
				auto adjch = mesh.NeighborCh(fh);
				if (adjch.size() == 2) continue;
				sum_normal = sum_normal + mesh.getQuadNormal(fh) * mesh.getQuadArea(fh);
			}
			normal_vertex[vp.first] = Eigen::Vector3d(sum_normal.x(), sum_normal.y(), sum_normal.z());
			normal_vertex[vp.first].normalize();
			double norm = normal_vertex[vp.first].norm();
			if (abs(norm - 1.f) > 1E-3f) normal_vertex.erase(vp.first);
		}

	}
	

}

void HexMesh_Regularization::LocalRegularization(const iGameCellHandle& ch) {

	auto& cell = mesh.cells(ch);
	auto fhs = cell.getFaceHandle();
	auto vhs = cell.getVertexHandle();
	auto ehs = cell.getEdgeHandle();

//#define WHERE_ARE_YOU

#if DEBUG
	//std::cout << "DEBUG" << std::endl;
	assert(fhs.size() == 6 && vhs.size() == 8);
#else
	//std::cout << "Release" << std::endl;
	if (fhs.size() != 6 || vhs.size() != 8) {
		std::cerr << "size is wrong!" << std::endl;
		return;
	}
#endif
	
	std::unordered_map<int, bool> cell_faces;// 记录这个 cell 的面
	std::unordered_map<int, bool> cell_vertices;// 记录这个 cell 的点
	std::unordered_map<int, Eigen::Vector3d> dual_node;// 建立面到 dual 顶点的映射
	std::unordered_map<int, Eigen::Vector3d> hex_node;
	std::unordered_map<int, int> local_id;// vh 在这个函数内的局部 id
	std::unordered_map<int, Eigen::Vector3d> ideal_pos;// save ideal positions

#ifdef WHERE_ARE_YOU
	std::cout << "you are here 1" << std::endl;
#endif // WHERE_ARE_YOU
	int tmp_cnt = 0;
	for (auto& vh : vhs) {
		cell_vertices[vh] = true;
		auto& v = mesh.vertices(vh);
		hex_node[vh] = Eigen::Vector3d(v.x(), v.y(), v.z());
		local_id[vh] = tmp_cnt++;// 用local_id.size()赋值可能会出错
	}

	for (auto& fh : fhs) {
		cell_faces[fh] = true;
	}

	for (int it = 0; it < 3; ++it) {
		for (auto& fh : fhs) {
			dual_node[fh] = Eigen::Vector3d::Zero();
			auto& face = mesh.faces(fh);
			auto f_vhs = face.getVertexHandle();
			for (auto& vh : f_vhs) {
				dual_node[fh] += hex_node[vh];
			}
			dual_node[fh] /= f_vhs.size();
		}
		for (auto& vh : vhs) {
			auto& v = hex_node[vh];
			vector<int> dual_id;
			for (auto& fh : mesh.NeighborFh(vh)) {
				if (cell_faces.count(fh)) {
					dual_id.push_back(fh);
				}
			}
#if __DEBUG__
			assert(dual_id.size() == 3);
#else
			if (dual_id.size() != 3) {
				std::cerr << "size is wrong!" << std::endl;
				return;
			}
#endif
			
			auto& v1 = dual_node[dual_id[0]];
			auto& v2 = dual_node[dual_id[1]];
			auto& v3 = dual_node[dual_id[2]];
			auto center = (v1 + v2 + v3) / 3;
			auto out_dir = (v - center).normalized();
			auto normal = ((v2 - v1).cross(v3 - v1)).normalized();
			if (out_dir.dot(normal) < 0.f) normal *= -1;
			ideal_pos[vh] = center + scaling_factor * normal;
		}
		hex_node = ideal_pos;
	}
#ifdef WHERE_ARE_YOU
	std::cout << "you are here 4" << std::endl;
#endif // WHERE_ARE_YOU
	

	Eigen::VectorXd B;
	B.resize(60);// 8 vertices and 12 edges, 3 * 8 + 3 * 12

	std::vector<Eigen::Triplet<double>> triplets;
	
	for (int i = 0; i < 12; ++i) {
		auto& eh = ehs[i];
		auto& edge = mesh.edges(eh);
		auto& vh1 = edge.vh1();
		auto& vh2 = edge.vh2();
		vector<int> l_id = { local_id[vh1] , local_id[vh2] };

		Eigen::Vector3d e = ideal_pos[vh2] - ideal_pos[vh1];
		e *= weight_e;
		for (int j = 0; j < 2; ++j) {
			auto id = l_id[j];
			if (j == 1) e *= -1;
			triplets.emplace_back(i * 3, id * 3 + 1, e.z());
			triplets.emplace_back(i * 3, id * 3 + 2, -e.y());
			triplets.emplace_back(i * 3 + 1, id * 3, e.z());
			triplets.emplace_back(i * 3 + 1, id * 3 + 2, -e.x());
			triplets.emplace_back(i * 3 + 2, id * 3, e.y());
			triplets.emplace_back(i * 3 + 2, id * 3 + 1, -e.x());
		}
		B[i * 3] = B[i * 3 + 1] = B[i * 3 + 2] = 0;
	}

#ifdef WHERE_ARE_YOU
	std::cout << "you are here 5" << std::endl;
#endif // WHERE_ARE_YOU

	for (auto& vh : vhs) {// the lower half one of vector B
		int i = local_id[vh];
		auto& v = mesh.vertices(vh);
		int row = 3 * i + 36;
		int col = 3 * i;
		B[row] = v.x();
		B[row + 1] = v.y();
		B[row + 2] = v.z();
		triplets.emplace_back(row, col, 1);
		triplets.emplace_back(row + 1, col + 1, 1);
		triplets.emplace_back(row + 2, col + 2, 1);
		//std::cout << "vh: " << row << ", " << col << std::endl;
	}
#ifdef WHERE_ARE_YOU
	std::cout << "you are here 6" << std::endl;
#endif // WHERE_ARE_YOU
	
	// 稀疏矩阵
	Eigen::SparseMatrix<double> A(60, 24);// 8 vertices and 12 edges, 3 * 8 + 3 * 12; 8 vertices
	A.setFromTriplets(triplets.begin(), triplets.end());

#ifdef WHERE_ARE_YOU
	std::cout << "you are here 7" << std::endl;
#endif // WHERE_ARE_YOU

	// 建立法方程
	Eigen::SparseMatrix<double> AT = A.transpose();
	Eigen::SparseMatrix<double> ATA = AT * A;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> MatricesCholesky(ATA);
	Eigen::VectorXd ATB = AT * B;

#ifdef WHERE_ARE_YOU
	std::cout << "you are here 8" << std::endl;
#endif // WHERE_ARE_YOU

	Eigen::VectorXd X = MatricesCholesky.solve(ATB);

#ifdef WHERE_ARE_YOU
	std::cout << "you are here 9" << std::endl;
#endif // WHERE_ARE_YOU

	for (int i = 0; i < 8; ++i) {
		auto& vh = vhs[i];
		auto& v = mesh.vertices(vh);
		assert(local_id[vh] == i);
		Eigen::Vector3d ideal_pos(X[i * 3], X[i * 3 + 1], X[i * 3 + 2]);
		transforms[vh].push_back(transform_node(ch, ideal_pos));
		/*printf("ch: %d, vh: %d, initial_pos: (%.2f, %.2f, %.2f), ideal_pos: (%.2f, %.2f, %.2f)\n",
			ch, vh, v.x(), v.y(), v.z(), ideal_pos.x(), ideal_pos.y(), ideal_pos.z());*/
	}

#ifdef WHERE_ARE_YOU
	std::cout << "you are here 10" << std::endl;
#endif // WHERE_ARE_YOU

}

void HexMesh_Regularization::GlobalOptimization() {


	auto& vps = mesh.allvertices();
	auto vcnt = mesh.VertexSize();
	unordered_map<int, int> vh_idx;
	vector<int> idx_vh;
	idx_vh.reserve(vcnt);
	for (auto& vp : vps) {
		vh_idx[vp.first] = idx_vh.size();
		idx_vh.push_back(vp.first);
	}

	assert(idx_vh.size() == vcnt && vh_idx.size() == vcnt);

	int normal_cnt = geometric_constraint_flag ? normal_vertex.size() : 0;

	Eigen::VectorXd B(6 * vcnt + normal_cnt);

	for (int i = 0; i < normal_cnt; ++i) {
		B[6 * vcnt + i] = 0;
	}

	for (int i = 0; i < vcnt; ++i) {
		auto& vh = idx_vh[i];
		auto& init_pos = mesh.vertices(iGameVertexHandle(vh));
		auto& ideal_pos = ideal_position[vh];
		B[i * 3] = weight_d * (ideal_pos.x() - init_pos.x());
		B[i * 3 + 1] = weight_d * (ideal_pos.y() - init_pos.y());
		B[i * 3 + 2] = weight_d * (ideal_pos.z() - init_pos.z());
		B[i * 6] = B[i * 6 + 1] = B[i * 6 + 2] = 0.f;
	}

	vector<Eigen::Triplet<double>> triplets;

	for (int i = 0; i < vcnt; ++i) {
		triplets.emplace_back(i * 3, i * 3, weight_d);
		triplets.emplace_back(i * 3 + 1, i * 3 + 1, weight_d);
		triplets.emplace_back(i * 3 + 2, i * 3 + 2, weight_d);
	}

	// init Laplacian matrix
	for (auto& vp : vps) {
		auto vh = vp.first;
		auto i = vh_idx[vh];
		const auto& adjehs = mesh.NeighborEh(vh);
		assert(adjehs.size() != 0);
		double sum_weight = 0.f;
		if (check_weight_flag) printf("vex: %d, weight: ", vp.first);
		for (auto& adjeh : adjehs) {
			auto& adje = mesh.edges(adjeh);
			double weight = weight_edge[adjeh];// tetrahedron cotangent weight
			int adjvh = adje.vh1() + adje.vh2() - vh;
			//weight = (mesh.vertices(vh) - mesh.vertices(iGameVertexHandle(adjvh))).norm();// length weight
			//weight = 1;// uniform weight
			auto j = vh_idx[adjvh];
			triplets.emplace_back(i * 6, j * 3, weight);
			triplets.emplace_back(i * 6 + 1, j * 3 + 1, weight);
			triplets.emplace_back(i * 6 + 2, j * 3 + 2, weight);
			if (check_weight_flag) printf("%.6f ", weight);
			sum_weight -= weight;
		}
		triplets.emplace_back(i * 6, i * 3, sum_weight);
		triplets.emplace_back(i * 6 + 1, i * 3 + 1, sum_weight);
		triplets.emplace_back(i * 6 + 2, i * 3 + 2, sum_weight);
		if (check_weight_flag) printf("; sum_weight: %.6f\n", sum_weight);
	}

	if (geometric_constraint_flag) {
		// init geometric constraint matrix
		int i = 0;
		int offset = 6 * vcnt;
		for (auto& np : normal_vertex) {
			int idx = vh_idx[np.first];
			auto& N = np.second;
			auto len = N.norm();
			assert(fabs(len - 1.f) < 1E-3F && i < normal_cnt);
			triplets.emplace_back(offset + i, idx * 3, weight_g * N.x());
			triplets.emplace_back(offset + i, idx * 3 + 1, weight_g * N.y());
			triplets.emplace_back(offset + i, idx * 3 + 2, weight_g * N.z());
			i++;
		}
		
	}
	
	time_t solve_start = clock();
	// 稀疏矩阵
	Eigen::SparseMatrix<double> A(6 * vcnt + normal_cnt, 3 * vcnt);
	A.setFromTriplets(triplets.begin(), triplets.end());

	// 建立法方程
	Eigen::SparseMatrix<double> AT = A.transpose();
	Eigen::SparseMatrix<double> ATA = AT * A;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> MatricesCholesky(ATA);

	Eigen::VectorXd ATB = AT * B;

	Eigen::VectorXd delta_X = MatricesCholesky.solve(ATB);

	time_t solve_end = clock();

	printf("Matrix %d rows, %d cols, solve cost: %dms\n", 6 * vcnt + normal_cnt, 3 * vcnt, int(solve_end - solve_start));

	for (int i = 0; i < vcnt; ++i) {
		auto vh = idx_vh[i];
		auto& v = mesh.vertices(iGameVertexHandle(vh));
		if (check_coordinate_flag) printf("v_pos: (%.2f, %.2f, %.2f) --> ", v.x(), v.y(), v.z());
		v.setX(v.x() + delta_X[i * 3]);
		v.setY(v.y() + delta_X[i * 3 + 1]);
		v.setZ(v.z() + delta_X[i * 3 + 2]);
		if (check_coordinate_flag) printf("(%.2f, %.2f, %.2f)\n", v.x(), v.y(), v.z());
	}

}

void HexMesh_Regularization::GlobalOptimization_Exterior() {
	// only consider the vertex who on the boundary
	auto& vps = mesh.allvertices();
	
	unordered_map<int, int> vh_idx;
	vector<int> idx_vh;
	for (auto& vp : vps) {
		if (mesh.isOnBoundary(vp.first)) {
			vh_idx[vp.first] = idx_vh.size();
			idx_vh.push_back(vp.first);
		}
	}
	assert(idx_vh.size() == vh_idx.size());

	auto exterior_vcnt = idx_vh.size();
	int normal_cnt = normal_vertex.size();

	Eigen::VectorXd B(6 * exterior_vcnt + normal_cnt);

	for (int i = 0; i < normal_cnt; ++i) {
		B[6 * exterior_vcnt + i] = 0;
	}

	for (int i = 0; i < exterior_vcnt; ++i) {
		auto& vh = idx_vh[i];
		auto& init_pos = mesh.vertices(iGameVertexHandle(vh));
		auto& ideal_pos = ideal_position[vh];
		B[i * 3] = weight_d * (ideal_pos.x() - init_pos.x());
		B[i * 3 + 1] = weight_d * (ideal_pos.y() - init_pos.y());
		B[i * 3 + 2] = weight_d * (ideal_pos.z() - init_pos.z());
		B[i * 6] = B[i * 6 + 1] = B[i * 6 + 2] = 0.f;
	}

	vector<Eigen::Triplet<double>> triplets;

	for (int i = 0; i < exterior_vcnt; ++i) {
		triplets.emplace_back(i * 3, i * 3, weight_d);
		triplets.emplace_back(i * 3 + 1, i * 3 + 1, weight_d);
		triplets.emplace_back(i * 3 + 2, i * 3 + 2, weight_d);
	}

	for (int i = 0; i < exterior_vcnt; ++i) {
		MeshKernel::iGameVertexHandle vh(idx_vh[i]);
		assert(mesh.isOnBoundary(vh));
		double sum_weight = 0.f;
		int count = 0;
		for (auto& adjeh : mesh.NeighborEh(vh)) {
			auto& e = mesh.edges(adjeh);
			int adjvh = e.vh1() + e.vh2() - vh;
			if (vh_idx.count(adjvh)) {
				int j = vh_idx[adjvh];
				double weight = weight_edge[adjeh];
				triplets.emplace_back(i * 6, j * 3, weight);
				triplets.emplace_back(i * 6 + 1, j * 3 + 1, weight);
				triplets.emplace_back(i * 6 + 2, j * 3 + 2, weight);
				sum_weight -= weight;
				count++;
			}
		}
		assert(count != 0);
		triplets.emplace_back(i * 6, i * 3, sum_weight);
		triplets.emplace_back(i * 6 + 1, i * 3 + 1, sum_weight);
		triplets.emplace_back(i * 6 + 2, i * 3 + 2, sum_weight);
	}

	if (geometric_constraint_flag) {
		// init geometric constraint matrix
		int i = 0;
		int offset = 6 * exterior_vcnt;
		for (auto& np : normal_vertex) {
			if (!vh_idx.count(np.first)) continue;
			int idx = vh_idx[np.first];
			auto& N = np.second;
			auto len = N.norm();
			assert(fabs(len - 1.f) < 1E-3F && i < normal_cnt);
			triplets.emplace_back(offset + i, idx * 3, weight_g * N.x());
			triplets.emplace_back(offset + i, idx * 3 + 1, weight_g * N.y());
			triplets.emplace_back(offset + i, idx * 3 + 2, weight_g * N.z());
			i++;
		}

	}


	time_t solve_start = clock();
	// 稀疏矩阵
	Eigen::SparseMatrix<double> A(6 * exterior_vcnt + normal_cnt, 3 * exterior_vcnt);
	A.setFromTriplets(triplets.begin(), triplets.end());

	// 建立法方程
	Eigen::SparseMatrix<double> AT = A.transpose();
	Eigen::SparseMatrix<double> ATA = AT * A;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> MatricesCholesky(ATA);

	Eigen::VectorXd ATB = AT * B;

	Eigen::VectorXd delta_X = MatricesCholesky.solve(ATB);

	time_t solve_end = clock();

	printf("Matrix %d rows, %d cols, solve cost: %dms\n", 6 * exterior_vcnt + normal_cnt, 3 * exterior_vcnt, int(solve_end - solve_start));

	for (int i = 0; i < exterior_vcnt; ++i) {
		auto vh = idx_vh[i];
		auto& v = mesh.vertices(iGameVertexHandle(vh));
		if (check_coordinate_flag) printf("exterior v_pos: (%.2f, %.2f, %.2f) --> ", v.x(), v.y(), v.z());
		v.setX(v.x() + delta_X[i * 3]);
		v.setY(v.y() + delta_X[i * 3 + 1]);
		v.setZ(v.z() + delta_X[i * 3 + 2]);
		if (check_coordinate_flag) printf("(%.2f, %.2f, %.2f)\n", v.x(), v.y(), v.z());
	}

	LaplacianSmoothing_Interior();

}

void HexMesh_Regularization::LaplacianSmoothing_Interior() {

	vector<iGameVertexHandle> interior_vhs;
	for (auto& vp : mesh.allvertices()) {
		if (!mesh.isOnBoundary(vp.first)) {
			interior_vhs.push_back(vp.first);
		}
	}

	int interior_vcnt = interior_vhs.size();
	vector<iGameVertex> positions(interior_vcnt, iGameVertex(0, 0, 0));

	for (int i = 0; i < interior_vcnt; ++i) {
		auto adjvhs = mesh.NeighborVh(interior_vhs[i]);
		for (auto& adjvh : adjvhs) {
			auto& adjv = mesh.vertices(adjvh);
			positions[i] = positions[i] + adjv;
		}
		positions[i] = positions[i] / adjvhs.size();
	}


	for (int i = 0; i < interior_vcnt; ++i) {
		auto& vh = interior_vhs[i];
		auto& v = mesh.vertices(vh);
		if (check_coordinate_flag) printf("interior v_pos: (%.2f, %.2f, %.2f) --> ", v.x(), v.y(), v.z());
		v = positions[i];
		if (check_coordinate_flag) printf("(%.2f, %.2f, %.2f)\n", v.x(), v.y(), v.z());
	}

}

void HexMesh_Regularization::ComputeTetrahedronLaplacianWeight(vector<iGameVertexHandle> vhs) {
	// Paper: Gradient, field based inhomogeneous volumetric mesh deformation for maxillofacial surgery simulation.Comput Graph 2009
	assert(vhs.size() == 4);

	auto& v0 = mesh.vertices(vhs[0]);
	auto& v1 = mesh.vertices(vhs[1]);
	auto& v2 = mesh.vertices(vhs[2]);
	auto& v3 = mesh.vertices(vhs[3]);
	auto eh1 = mesh.getEdgeHandle(vhs[0], vhs[1]);
	auto eh2 = mesh.getEdgeHandle(vhs[0], vhs[2]);
	auto eh3 = mesh.getEdgeHandle(vhs[0], vhs[3]);
	assert(eh1 != -1 && eh2 != -1 && eh3 != -1);

	auto vec12 = v2 - v1;
	auto vec13 = v3 - v1;
	auto vec03 = v3 - v0;
	auto vec02 = v2 - v0;
	auto vec01 = v1 - v0;
	auto normal_123 = (vec12 % vec13).normalized();
	auto normal_032 = (vec03 % vec02).normalized();
	auto normal_021 = (vec02 % vec01).normalized();
	auto normal_013 = (vec01 % vec03).normalized();

	// note: the consine of normals equal to -consine of faces
	double cosine = normal_123 * normal_032 * -1;
	double cotangent = max(cosine / sqrt(1 - cosine * cosine), KEEPOPSITIVE);
	//printf("cot: %.2f\n", cotangent);
	weight_edge[eh1] += (factor_e * cotangent * (v3 - v2).norm());

	cosine = normal_123 * normal_013 * -1;
	cotangent = max(cosine / sqrt(1 - cosine * cosine), KEEPOPSITIVE);
	weight_edge[eh2] += (factor_e * cotangent * (v3 - v1).norm());
	//printf("cot: %.2f\n", cotangent);
	cosine = normal_123 * normal_021 * -1;
	cotangent = max(cosine / sqrt(1 - cosine * cosine), KEEPOPSITIVE);
	weight_edge[eh3] += (factor_e * cotangent * (v2 - v1).norm());
	//printf("cot: %.2f\n", cotangent);
}

void HexMesh_Regularization::CalculateCellQuality(vector<iGameCellHandle>& chs, int mode) {
	// compute result cell quality
	double sum_sj = 0, min_sj = 2.f;
	for (int i = 0; i < chs.size(); ++i) {
		auto& ch = chs[i];
		auto& c = mesh.cells(ch);
		double min_det = 2.f;// the value of Scaled Jacobian is the	determinant of the matrix Jacobian
		unordered_map<int, vector<int>> neighbor;
		const auto& _ehs = c.getEdgeHandle();
		for (auto& _eh : _ehs) {
			auto& _e = mesh.edges(_eh);
			auto& vh1 = _e.vh1();
			auto& vh2 = _e.vh2();
			neighbor[vh1].push_back(vh2);
			neighbor[vh2].push_back(vh1);
		}
		const auto& _vhs = c.getVertexHandle();
		for (auto& _vh : _vhs) {
			auto& v = mesh.vertices(_vh);
			auto adjvh = neighbor[_vh];
			assert(adjvh.size() == 3);
			auto v1 = mesh.vertices(iGameVertexHandle(adjvh[0]));
			auto v2 = mesh.vertices(iGameVertexHandle(adjvh[1]));
			auto v3 = mesh.vertices(iGameVertexHandle(adjvh[2]));
			auto center = (v1 + v2 + v3) / 3;
			auto vec = (center - v).normalized();

			// ajust oridering
			auto vec12 = v2 - v1, vec13 = v3 - v1;
			auto normal = (vec12 % vec13).normalized();
			double cosine = vec * normal;
			if (cosine < 0) swap(v2, v3);

			vector<Eigen::Vector3d> ev(3);
			auto tmp = v1 - v;
			ev[0] = (Eigen::Vector3d(tmp.x(), tmp.y(), tmp.z())).normalized();
			tmp = v2 - v;
			ev[1] = (Eigen::Vector3d(tmp.x(), tmp.y(), tmp.z())).normalized();
			tmp = v3 - v;
			ev[2] = (Eigen::Vector3d(tmp.x(), tmp.y(), tmp.z())).normalized();

			Eigen::Matrix3d J;
			J << ev[0][0], ev[1][0], ev[2][0],
				ev[0][1], ev[1][1], ev[2][1],
				ev[0][2], ev[1][2], ev[2][2];
			//std::cout << J << endl;
			/*
			For an element, the full range of the Scaled Jacobian value is from −1 to + 1.
			And fscaled_jacobian = 1 if the element is an ideal element,
			fscaled_jacobian = −1 if the element is a worst distorted element
			*/
			min_det = min(min_det, J.determinant());
		}
		//printf("weight of cell %d: %.4f\n", ch, min_det);
		sum_sj += min_det;
		min_sj = min(min_det, min_sj);
	}
	if (mode == 0) {
		in_avg_scaled_jacobian = sum_sj / chs.size();
		in_min_scaled_jacobian = min_sj;
	} else if (mode == 1) {
		out_avg_scaled_jacobian = sum_sj / chs.size();
		out_min_scaled_jacobian = min_sj;
	}
	
}