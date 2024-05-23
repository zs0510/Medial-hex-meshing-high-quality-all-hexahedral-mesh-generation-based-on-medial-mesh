#include "Smoothing.h"

HexMesh_Smoothing::HexMesh_Smoothing(MeshKernel::VolumeMesh& _mesh) : mesh(_mesh) {

}

void HexMesh_Smoothing::Laplacian_Smoothing_HC(std::vector<int> feature_vhs, int iter_max) {

	double alpha = 0, beta = 0.5;
	double volume_origin = MeshMath::get_volume_hexahedral_mesh(mesh);
	int n_vertices = mesh.vsize();
	exterior_vhs.clear();
	interior_vhs.clear();
	vector<iGameVertex> vertices_origin;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) exterior_vhs.push_back(vp.first);
		else interior_vhs.push_back(vp.first);
		vertices_origin.push_back(vp.second);
	}

	int exterior_cnt = exterior_vhs.size(), interior_cnt = interior_vhs.size();

	for (int it = 0; it < iter_max; ++it) {

		vector<iGameVertex> differ(n_vertices);
		vector<iGameVertex> vertices_smoothed(n_vertices);

#pragma omp parallel for
		for (int i = 0; i < exterior_cnt; ++i) {
			auto vh = exterior_vhs[i];
			int ex_cnt = 0;
			iGameVertex adj_v(0, 0, 0);
			for (auto& adjvh : mesh.NeighborVh(vh)) {
				if (mesh.isOnBoundary(adjvh)) {
					adj_v += mesh.vertices(adjvh);
					ex_cnt++;
				}
			}
			vertices_smoothed[vh] = adj_v / ex_cnt;
			differ[vh] = vertices_smoothed[vh] - (vertices_origin[vh] * alpha + mesh.vertices(vh) * (1 - alpha));
		}

#pragma omp parallel for
		for (int i = 0; i < exterior_cnt; ++i) {
			auto vh = exterior_vhs[i];
			int ex_cnt = 0;
			iGameVertex adj_diff(0, 0, 0);
			for (auto& adjvh : mesh.NeighborVh(vh)) {
				if (mesh.isOnBoundary(adjvh)) {
					adj_diff += differ[adjvh];
					ex_cnt++;
				}
			}
			adj_diff /= ex_cnt;
			vertices_smoothed[vh] = vertices_smoothed[vh] - (differ[vh] * beta + adj_diff * (1 - beta));
			auto& v = mesh.vertices(vh);
			v = vertices_smoothed[vh];
		}

#pragma omp parallel for
		for (int i = 0; i < interior_cnt; ++i) {
			auto vh = interior_vhs[i];
			const auto& adjvhs = mesh.NeighborVh(vh);
			iGameVertex adj_v(0, 0, 0);
			for (auto& adjvh : adjvhs) {
				adj_v += mesh.vertices(adjvh);
			}
			vertices_smoothed[vh] = adj_v / adjvhs.size();
		}

#pragma omp parallel for
		for (int i = 0; i < interior_cnt; ++i) {
			auto vh = interior_vhs[i];
			auto& v = mesh.vertices(vh);
			v = vertices_smoothed[vh];
		}

	}

	double volume_current = MeshMath::get_volume_hexahedral_mesh(mesh);

	std::cout << "Smooth(HC): Volume_origin = " << volume_origin << ", Volume_current = " << volume_current << std::endl;

}

void HexMesh_Smoothing::Laplacian_Smoothing_Diffusion(int max_iter, std::vector<int> feature_cells, std::vector<int> feature_faces,
	std::vector<int> feature_edges, std::vector<int> feature_vertices) {

	double original_volume = MeshMath::get_volume_hexahedral_mesh(mesh);
	int n_vertices = mesh.vsize();
	vector<bool> is_feature(n_vertices, false);
	for (auto& ch : feature_cells) {
		iGameCellHandle _ch(ch);
		auto& cell = mesh.cells(_ch);
		for (auto& vh : cell.getVertexHandle()) {
			is_feature[vh] = true;
		}
	}
	for (auto& fh : feature_faces) {
		iGameFaceHandle _fh(fh);
		auto& face = mesh.faces(_fh);
		for (auto& vh : face.getVertexHandle()) {
			is_feature[vh] = true;
		}
	}
	for (auto& eh : feature_edges) {
		iGameEdgeHandle _eh(eh);
		auto& edge = mesh.edges(_eh);
		is_feature[edge.vh1()] = true;
		is_feature[edge.vh2()] = true;
	}
	for (auto& vh : feature_vertices) {
		is_feature[iGameVertexHandle(vh)] = true;
	}

	vector<double> weight_s(n_vertices, 1.f);// 采用渐进扩散，越靠近特征顶点该权重越小
	vector<bool> visited(n_vertices, false);
	queue<int> que;
	int feature_cnt = 0;
	for (int i = 0; i < n_vertices; ++i) {
		if (is_feature[i]) {
			feature_cnt++;
			visited[i] = true;
			weight_s[i] = 0.f;
			que.push(i);
		}
	}
	while (!que.empty()) {
		auto vid = que.front();
		que.pop();
		MeshKernel::iGameVertexHandle vh(vid);
		double weight_cur = weight_s[vid];
		for (auto& adjvh : mesh.NeighborVh(vh)) {
			if (visited[adjvh]) continue;
			visited[adjvh] = true;
			weight_s[adjvh] = std::min((double)1.f, weight_cur + 0.01f);
			if (weight_s[adjvh] < 1.f - (1E-4F)) {
				que.push(adjvh);
			}
		}
	}
	/*for (int i = 0; i < n_vertices; ++i) {
		if (weight_s[i] < 1.f)
			std::cout << "#V " << i << ": " << weight_s[i] << std::endl;
	}*/


	unordered_map<int, bool> is_on_boundary;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			exterior_vhs.push_back(vp.first);
			is_on_boundary[vp.first] = true;
		} else {
			interior_vhs.push_back(vp.first);
		}
	}
	int exterior_vcnt = exterior_vhs.size(), interior_vcnt = interior_vhs.size();

	for (int it = 0; it < max_iter; ++it) {

		vector<iGameVertex> total_vertices;
		vector<iGameVertex> positions(exterior_vcnt);

#pragma omp parallel for
		for (int i = 0; i < exterior_vcnt; ++i) {
			auto& vh = exterior_vhs[i];
			iGameVertex pos(0, 0, 0);
			auto adjvhs = mesh.NeighborVh(vh);
			int vcnt = 0;
			for (auto& adjvh : adjvhs) {
				if (is_on_boundary.count(adjvh)) {
					pos = pos + mesh.vertices(adjvh);
					vcnt++;
				}
			}
			positions[i] = pos / vcnt;
		}
		total_vertices.insert(total_vertices.end(), positions.begin(), positions.end());

#pragma omp parallel for
		for (int i = 0; i < exterior_vcnt; ++i) {
			if (is_feature[exterior_vhs[i]]) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(exterior_vhs[i]);
			v = v * (1 - weight_s[exterior_vhs[i]]) + positions[i] * weight_s[exterior_vhs[i]];
		}

		positions.resize(interior_vcnt);

#pragma omp parallel for
		for (int i = 0; i < interior_vcnt; ++i) {
			auto& vh = interior_vhs[i];
			iGameVertex pos(0, 0, 0);
			auto adjvhs = mesh.NeighborVh(vh);
			for (auto& adjvh : adjvhs) {
				pos = pos + mesh.vertices(adjvh);
			}
			positions[i] = pos / adjvhs.size();
		}
		total_vertices.insert(total_vertices.end(), positions.begin(), positions.end());

#pragma omp parallel for
		for (int i = 0; i < interior_vcnt; ++i) {
			if (is_feature[interior_vhs[i]]) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(interior_vhs[i]);
			v = v * (1 - weight_s[interior_vhs[i]]) + positions[i] * weight_s[interior_vhs[i]];
		}

		// rescale
		double current_volume = MeshMath::get_volume_hexahedral_mesh(mesh);
		double scale = std::pow(original_volume / current_volume, 1.f / 3);
#pragma omp parallel for
		for (int i = 0; i < mesh.vsize(); ++i) {
			if (is_feature[i]) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(MeshKernel::iGameVertexHandle(i));
			v *= scale;
		}

	}

	std::cout << "Smoothing success. Feature vertices size = " << feature_cnt << std::endl;
}

void HexMesh_Smoothing::Laplacian_Smoothing(int max_iter, std::vector<int> feature_cells, std::vector<int> feature_faces, 
	std::vector<int> feature_edges, std::vector<int> feature_vertices) {

	if (feature_cells.empty() && feature_edges.empty() && feature_faces.empty() && feature_vertices.empty()) {
		update_bbox_center();
	}
	

	double volume_origin = MeshMath::get_volume_hexahedral_mesh(mesh);

	unordered_map<int, bool> feature_vhs;
	for (auto& ch : feature_cells) {
		iGameCellHandle _ch(ch);
		auto& cell = mesh.cells(_ch);
		for (auto& vh : cell.getVertexHandle()) {
			feature_vhs[vh] = true;
		}
	}
	for (auto& fh : feature_faces) {
		iGameFaceHandle _fh(fh);
		auto& face = mesh.faces(_fh);
		for (auto& vh : face.getVertexHandle()) {
			feature_vhs[vh] = true;
		}
	}
	for (auto& eh : feature_edges) {
		iGameEdgeHandle _eh(eh);
		auto& edge = mesh.edges(_eh);
		feature_vhs[edge.vh1()] = true;
		feature_vhs[edge.vh2()] = true;
	}
	for (auto& vh : feature_vertices) {
		feature_vhs[iGameVertexHandle(vh)] = true;
	}

	unordered_map<int, bool> is_on_boundary;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			exterior_vhs.push_back(vp.first);
			is_on_boundary[vp.first] = true;
		} else {
			interior_vhs.push_back(vp.first);
		}
	}
	int exterior_vcnt = exterior_vhs.size(), interior_vcnt = interior_vhs.size();

	for (int it = 0; it < max_iter; ++it) {

		vector<iGameVertex> total_vertices;
		vector<iGameVertex> positions(exterior_vcnt);

#pragma omp parallel for
		for (int i = 0; i < exterior_vcnt; ++i) {
			auto& vh = exterior_vhs[i];
			iGameVertex pos(0, 0, 0);
			auto adjvhs = mesh.NeighborVh(vh);
			int vcnt = 0;
			for (auto& adjvh : adjvhs) {
				if (is_on_boundary.count(adjvh)) {
					pos = pos + mesh.vertices(adjvh);
					vcnt++;
				}
			}
			positions[i] = pos / vcnt;
		}
		total_vertices.insert(total_vertices.end(), positions.begin(), positions.end());

#pragma omp parallel for
		for (int i = 0; i < exterior_vcnt; ++i) {
			if (feature_vhs.count(exterior_vhs[i])) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(exterior_vhs[i]);
			v = positions[i];
		}

		positions.resize(interior_vcnt);

#pragma omp parallel for
		for (int i = 0; i < interior_vcnt; ++i) {
			auto& vh = interior_vhs[i];
			iGameVertex pos(0, 0, 0);
			auto adjvhs = mesh.NeighborVh(vh);
			for (auto& adjvh : adjvhs) {
				pos = pos + mesh.vertices(adjvh);
			}
			positions[i] = pos / adjvhs.size();
		}
		total_vertices.insert(total_vertices.end(), positions.begin(), positions.end());

#pragma omp parallel for
		for (int i = 0; i < interior_vcnt; ++i) {
			if (feature_vhs.count(interior_vhs[i])) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(interior_vhs[i]);
			v = positions[i];
		}

		// rescale
		double current_volume = MeshMath::get_volume_hexahedral_mesh(mesh);
		double scale = std::pow(volume_origin / current_volume, 1.f / 3);
#pragma omp parallel for
		for (int i = 0; i < mesh.vsize(); ++i) {
			if (feature_vhs.count(i)) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(MeshKernel::iGameVertexHandle(i));
			v *= scale;
		}

	}
	double volume_current = MeshMath::get_volume_hexahedral_mesh(mesh);

	if (feature_cells.empty() && feature_edges.empty() && feature_faces.empty() && feature_vertices.empty()) {
		move_to_original_center();
	}
	

	std::cout << "Smooth: Volume_origin = " << volume_origin << ", Volume_current = " << volume_current << std::endl;

	//std::cout << "Smoothing success. Feature vertices size = " << feature_vhs.size() << std::endl;
}

void HexMesh_Smoothing::Laplacian_Level(int iter_max, bool smooth_boundary) {

	double volume_origin = MeshMath::get_volume_hexahedral_mesh(mesh);

	int vcnt = mesh.vsize();
	vector<bool> visited(vcnt, false);
	vector<vector<VH>> level_vhs;// 从外至内
	vector<VH> boundary_vhs;
	queue<VH> que;

	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			visited[vp.first] = true;
			boundary_vhs.emplace_back(vp.first);
			que.emplace(vp.first);
		}
	}

	if (smooth_boundary) {

		vector<Vex> positions(boundary_vhs.size());

		for (int it = 0; it < iter_max; ++it) {

#pragma omp parallel for
			for (int i = 0; i < boundary_vhs.size(); ++i) {
				auto& vh = boundary_vhs[i];
				const auto& adjvhs = mesh.NeighborVh(vh);
				positions[i] = Vex(0, 0, 0);
				int bdy_cnt = 0;
				for (auto& adjvh : adjvhs) {
					if (visited[adjvh]) {// 此时只有边界顶点为 visited
						auto& v = mesh.vertices(adjvh);
						positions[i] += v;
						bdy_cnt++;
					}
				}
				if (bdy_cnt == 0) positions[i] = mesh.vertices(vh);
				else positions[i] /= bdy_cnt;
			}

#pragma omp parallel for
			for (int i = 0; i < boundary_vhs.size(); ++i) {
				auto& v = mesh.vertices(boundary_vhs[i]);
				v = positions[i];
			}


		}




	}

	while (!que.empty()) {// 层序遍历
		vector<VH> vhs;
		int sz = que.size();
		for (int i = 0; i < sz; ++i) {
			VH vh = que.front();
			que.pop();
			const auto& adjvhs = mesh.NeighborVh(vh);
			for (auto& adjvh : adjvhs) {
				if (visited[adjvh]) continue;
				visited[adjvh] = true;
				vhs.emplace_back(adjvh);
				que.emplace(adjvh);
			}
		}
		level_vhs.emplace_back(vhs);
	}
	
	for (int it = 0; it < iter_max; ++it) {

		for (auto& current_level_vhs : level_vhs) {

			int c_vcnt = current_level_vhs.size();
			vector<Vex> positions(c_vcnt, Vex(0, 0, 0));

#pragma omp parallel for
			for (int i = 0; i < c_vcnt; ++i) {
				VH vh(current_level_vhs[i]);
				positions[i] = Vex(0, 0, 0);
				const auto& adjvhs = mesh.NeighborVh(vh);
				for (auto& adjvh : adjvhs) {
					auto& adjv = mesh.vertices(adjvh);
					positions[i] += adjv;
				}

				if (adjvhs.empty()) positions[i] = mesh.vertices(vh);
				else positions[i] /= adjvhs.size();

			}

#pragma omp parallel for
			for (int i = 0; i < c_vcnt; ++i) {
				auto& v = mesh.vertices(current_level_vhs[i]);
				v = positions[i];
			}

		}

	}

	// rescale
	double volume_current = MeshMath::get_volume_hexahedral_mesh(mesh);
	double scale = std::pow(volume_origin / volume_current, 1.f / 3);

#pragma omp parallel for
	for (int i = 0; i < mesh.vsize(); ++i) {
		auto& v = mesh.vertices(MeshKernel::iGameVertexHandle(i));
		v *= scale;
	}


}

double HexMesh_Smoothing::CalculateCellQuality(MeshKernel::iGameCellHandle ch) {
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
	return min_det;
}

void SurfaceMesh_Smoothing::Laplacian_smoothing(int iter_times) {

	int vcnt = mesh.vsize();
	std::vector<Vex> positions(vcnt);

	for (int it = 0; it < iter_times; ++it) {

#pragma omp parallel for
		for (int vi = 0; vi < vcnt; ++vi) {
			VH vh(vi);
			const auto& adjvhs = mesh.NeighborVh(vh);
			Vex pos_sum(0, 0, 0);
			for (auto& adjvh : adjvhs) {
				auto& adjv = mesh.vertices(adjvh);
				pos_sum += adjv;
			}
			if (adjvhs.empty()) positions[vi] = mesh.vertices(vh);
			else positions[vi] = pos_sum / adjvhs.size();
		}

#pragma omp parallel for
		for (int vi = 0; vi < vcnt; ++vi) {
			VH vh(vi);
			auto& v = mesh.vertices(vh);
			v = positions[vi];
		}

	}



}