#include "Remeshing.h"

Remeshing::Remeshing(MeshKernel::SurfaceMesh& _mesh, bool _uniform_flag) : mesh(_mesh), uniform_flag(_uniform_flag) {

	initAABBTree();
	initTargetEdgeLength();

}

void Remeshing::Execute() {

	time_t exe_start = clock();
	time_t before, after;

	if (!uniform_flag) {

		for (int it = 0; it < 5; ++it) {

			convergence_flag = true;

			int k = mesh.FaceSize() * 0.2;
			printf("\niter %d: vertices_size = %d, faces_size = %d, k = %d\n", it, mesh.VertexSize(), mesh.FaceSize(), k);

			before = clock();
			removeLargeAngle(k);
			after = clock();
			printf("iter %d: remove large angle success, cost: %dms\n", it, int(after - before));

			before = clock();
			equalizeValence();
			after = clock();
			printf("iter %d: equalize valence success, cost: %dms\n", it, int(after - before));

			before = clock();
			tangentialRelaxation(3);
			after = clock();
			printf("iter %d: tangential relaxation success, cost: %dms\n", it, int(after - before));

			before = clock();
			removeSmallAngle(k);
			after = clock();
			printf("iter %d: remove small angle success, cost: %dms\n", it, int(after - before));

			before = clock();
			equalizeValence();
			after = clock();
			printf("iter %d: equalize valence success, cost: %dms\n", it, int(after - before));

			before = clock();
			tangentialRelaxation(3);
			after = clock();
			printf("iter %d: tangential relaxation success, cost: %dms\n", it, int(after - before));

			while (1) {
				int singular_cnt = 0;
				auto vertices = mesh.allvertices();
				for (auto& vp : vertices) {
					auto vh = vp.first;
					if (mesh.NeighborVh(vh).empty()) {
						mesh.DeleteVertex(vh);
						singular_cnt++;
					}
				}
				if (singular_cnt == 0) break;
				else printf("iter %d: find %d singular vertices\n", it, singular_cnt);
			}
			if (convergence_flag) break;
		}

	} else {

		for (int it = 0; it < 5; ++it) {

			convergence_flag = true;

			before = clock();
			splitLongEdge();
			after = clock();
			printf("iter %d: split long edge success, cost: %dms\n", it, int(after - before));

			before = clock();
			collapseShortEdge();
			after = clock();
			printf("iter %d: collapse short edge success, cost: %dms\n", it, int(after - before));

			before = clock();
			equalizeValence();
			after = clock();
			printf("iter %d: equalize valence success, cost: %dms\n", it, int(after - before));


			before = clock();
			tangentialRelaxation(3);
			after = clock();
			printf("iter %d: tangential relaxation success, cost: %dms\n", it, int(after - before));

			while (1) {
				int singular_cnt = 0;
				auto vertices = mesh.allvertices();
				for (auto& vp : vertices) {
					auto vh = vp.first;
					if (mesh.NeighborVh(vh).empty()) {
						mesh.DeleteVertex(vh);
						singular_cnt++;
					}
				}
				if (singular_cnt == 0) break;
				else printf("iter %d: find %d singular vertices\n", it, singular_cnt);
			}
			if (convergence_flag) break;
		}

	}
	
	
	
	time_t exe_end = clock();
	printf("remeshing success, cost: %dms\n", int(exe_end - exe_start));

	mesh.updateAllHandles();

}

void Remeshing::initAABBTree() {
	// bulid aabb tree
	std::vector<Vector3f> vertices;
	for (auto& fp : mesh.allfaces()) {
		auto vhs = fp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& v = mesh.vertices(vh);
			vertices.push_back(Vector3f(v.x(), v.y(), v.z()));
		}
	}
	ab_tree = new AABB_Tree(vertices);
}

void Remeshing::initTargetEdgeLength() {
	// get target edge length
	double area_sum = 0.f;
	auto& fps = mesh.allfaces();
	for (auto& fp : fps) {
		auto vhs = fp.second.getVertexHandle();
		auto& v1 = mesh.vertices(vhs[0]);
		auto& v2 = mesh.vertices(vhs[1]);
		auto& v3 = mesh.vertices(vhs[2]);
		auto vec12 = v2 - v1;
		auto vec13 = v3 - v1;
		double area = (vec12 % vec13).norm() * 0.5f;
		area_sum += area;
	}

	target_length = sqrt(area_sum / (1.5 * mesh.VertexSize()));// target number of vertices of output mesh
	target_length *= (2 / pow(3, 0.25f));
	target_length *= 0.5;
	printf("area_average_length: %.2f\n", target_length);
	lower_length = target_length * 4 / 5;
	upper_length = target_length * 4 / 3;
}

void Remeshing::initEdgesNorm2() {
	edge_norm2.clear();
	for (auto& ep : mesh.alledges()) {
		calEdgeNorm2(ep.first);
	}
}

void Remeshing::calEdgeNorm2(MeshKernel::iGameEdgeHandle eh) {
	auto& edge = mesh.edges(eh);
	auto& v1 = mesh.vertices(edge.vh1());
	auto& v2 = mesh.vertices(edge.vh2());
	edge_norm2[eh] = (v2 - v1).norm2();
}

void Remeshing::removeLargeAngle(int count) {
	int remove_cnt = 0;
	int flip_cnt = 0;
	auto oldFaces = mesh.allfaces();
	initEdgesNorm2();
	for (auto& fp : oldFaces) {
		if (!mesh.isValid(fp.first)) continue;
		auto fh = fp.first;
		auto& face = mesh.faces(fh);
		auto ehs = face.getEdgeHandle();
		auto vhs = face.getVertexHandle();
		assert(edge_norm2.count(ehs[0]) && edge_norm2.count(ehs[1]) && edge_norm2.count(ehs[2]));
		MeshKernel::iGameEdgeHandle longest_eh = ehs[0];
		if (edge_norm2[ehs[1]] > edge_norm2[longest_eh]) longest_eh = ehs[1];
		if (edge_norm2[ehs[2]] > edge_norm2[longest_eh]) longest_eh = ehs[2];

		auto& longest_e = mesh.edges(longest_eh);
		auto vh1 = longest_e.vh1();
		auto vh2 = longest_e.vh2();
		MeshKernel::iGameVertexHandle vh3(vhs[0] + vhs[1] + vhs[2] - vh1 - vh2);
		auto& v1 = mesh.vertices(vh1);
		auto& v2 = mesh.vertices(vh2);
		auto& v3 = mesh.vertices(vh3);
		auto vec31 = (v1 - v3).normalized();
		auto vec32 = (v2 - v3).normalized();
		double angle = std::acos(vec31 * vec32)  * rad2angle;
		if (angle < upper_angle) continue;// the angle is in the bound
		MeshKernel::iGameVertex new_v = (v1 + v2) / 2;
		
		auto adjfhs = mesh.NeighborFh(longest_eh);
		if (adjfhs.size() == 1) continue;
		std::vector<MeshKernel::iGameEdgeHandle> quadehs;
		for (auto& fh : adjfhs) {
			for (auto& f_eh : mesh.faces(fh).getEdgeHandle()) {
				if (f_eh != longest_eh)
					quadehs.push_back(f_eh);
			}
		}
		auto new_vh = mesh.AddVertex(new_v);
		for (auto& fh : adjfhs) {// adjacent triangles are bisected
			auto newFace = mesh.faces(fh).getVertexHandle();
			if (newFace[0] == vh1) newFace[0] = new_vh;
			else if (newFace[1] == vh1) newFace[1] = new_vh;
			else if (newFace[2] == vh1) newFace[2] = new_vh;
			mesh.AddFace(newFace);
			newFace = mesh.faces(fh).getVertexHandle();
			if (newFace[0] == vh2) newFace[0] = new_vh;
			else if (newFace[1] == vh2) newFace[1] = new_vh;
			else if (newFace[2] == vh2) newFace[2] = new_vh;
			mesh.AddFace(newFace);
		}
		mesh.DeleteEdge(longest_eh);
		edge_norm2.erase(longest_eh);
		for (auto& new_eh : mesh.NeighborEh(new_vh)) {
			calEdgeNorm2(new_eh);
		}
		MeshKernel::iGameEdgeHandle flip_eh(-1);
		double min_square = Double_MAX;
		for (auto& q_eh : quadehs) {
			assert(mesh.isValid(q_eh));
			if (mesh.isOnBoundary(q_eh)) continue;
			if (flip_check_flag) {
				// check dihedral angle
				auto tmp_fhs = mesh.NeighborFh(q_eh);
				std::vector<MeshKernel::iGameVertex> normals;
				if (tmp_fhs.size() != 2) continue;
				for (auto& tmp_fh : tmp_fhs) {
					mesh.genNormal(tmp_fh);
					auto& tmp_f = mesh.faces(tmp_fh);
					auto n = MeshKernel::iGameVertex(tmp_f.getNormalX(), tmp_f.getNormalY(), tmp_f.getNormalZ());
					n = n.normalized();
					normals.push_back(n);
				}
				double angle = std::acos(normals[0] * normals[1]) * rad2angle;
				if (angle > flip_ok_angle) continue;// dihedral angle is too small(maybe a feature)
			}
			double square = get_square_difference(q_eh);
			if (square < min_square) {
				min_square = square;
				flip_eh = q_eh;
			}
		}
		if (is_flip_ok(flip_eh)) {
			flip(flip_eh);
			//printf("the square difference is: %.2f\n", min_square);
			flip_cnt++;
		}


		//printf("remove angle: %.2f\n", angle);
		if (++remove_cnt == count) break;
	}
	if (remove_cnt != 0) convergence_flag = false;
	printf("remove %d large angles, flip %d edges\n", remove_cnt, flip_cnt);

}

void Remeshing::removeSmallAngle(int count) {
	int remove_cnt = 0;
	auto oldFaces = mesh.allfaces();
	initEdgesNorm2();
	for (auto& fp : oldFaces) {
		if (!mesh.isValid(fp.first)) continue;
		auto fh = fp.first;
		auto& face = mesh.faces(fh);
		auto ehs = face.getEdgeHandle();
		auto vhs = face.getVertexHandle();
		if (!edge_norm2.count(ehs[0])) calEdgeNorm2(ehs[0]);
		if (!edge_norm2.count(ehs[1])) calEdgeNorm2(ehs[1]);
		if (!edge_norm2.count(ehs[2])) calEdgeNorm2(ehs[2]);
		assert(edge_norm2.count(ehs[0]) && edge_norm2.count(ehs[1]) && edge_norm2.count(ehs[2]));
		MeshKernel::iGameEdgeHandle shortest_eh = ehs[0];
		if (edge_norm2[ehs[1]] < edge_norm2[shortest_eh]) shortest_eh = ehs[1];
		if (edge_norm2[ehs[2]] < edge_norm2[shortest_eh]) shortest_eh = ehs[2];
		auto& shortest_e = mesh.edges(shortest_eh);
		auto vh1 = shortest_e.vh1();
		auto vh2 = shortest_e.vh2();
		MeshKernel::iGameVertexHandle vh3(vhs[0] + vhs[1] + vhs[2] - vh1 - vh2);
		auto& v1 = mesh.vertices(vh1);
		auto& v2 = mesh.vertices(vh2);
		auto& v3 = mesh.vertices(vh3);
		
		auto vec31 = (v1 - v3).normalized();
		auto vec32 = (v2 - v3).normalized();
		double angle = std::acos(vec31 * vec32)  * rad2angle;
		if (angle > lower_angle) continue;// the angle is in the bound
		auto vh_keep = vh1;
		auto vh_remove = vh2;
		auto& vex_keep = mesh.vertices(vh_keep);
		auto& vex_remove = mesh.vertices(vh_remove);
		// preserve boundary
		if (mesh.isOnBoundary(vh_keep)) {
			// do nothing
		} else if (mesh.isOnBoundary(vh_remove)) {
			vex_keep = vex_remove;
		} else {
			vex_keep = (vex_keep + vex_remove) / 2;
		}
		// erase invalid edge norm2
		auto old_ehs1 = mesh.NeighborEh(vh_keep);
		auto old_ehs2 = mesh.NeighborEh(vh_remove);
		for (auto& old_eh : old_ehs1) 
			edge_norm2.erase(old_eh);
		for (auto& old_eh : old_ehs2)
			edge_norm2.erase(old_eh);

		std::vector<std::vector<MeshKernel::iGameVertexHandle>> new_triangles;
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
		for (auto& new_eh : mesh.NeighborEh(vh_keep)) {
			calEdgeNorm2(new_eh);
		}
		//printf("remove angle: %.2f\n", angle);
		if (++remove_cnt == count) break;
	}
	if (remove_cnt != 0) convergence_flag = false;
	printf("remove %d small angles\n", remove_cnt);

}

void Remeshing::equalizeValence() {
	int flip_cnt = 0;
	auto oldEdges = mesh.alledges();
	for (auto& ep : oldEdges) {
		auto eh = ep.first;
		if (!mesh.isValid(eh) || mesh.isOnBoundary(eh)) continue;
		auto& e = mesh.edges(eh);

		if (flip_check_flag) {
			// check dihedral angle
			auto tmp_fhs = mesh.NeighborFh(eh);
			std::vector<MeshKernel::iGameVertex> normals;
			if (tmp_fhs.size() != 2) continue;
			for (auto& tmp_fh : tmp_fhs) {
				mesh.genNormal(tmp_fh);
				auto& tmp_f = mesh.faces(tmp_fh);
				auto n = MeshKernel::iGameVertex(tmp_f.getNormalX(), tmp_f.getNormalY(), tmp_f.getNormalZ());
				n = n.normalized();
				normals.push_back(n);
			}
			double angle = std::acos(normals[0] * normals[1])  * rad2angle;
			if (angle > flip_ok_angle) continue;// dihedral angle is too small(maybe a feature)
		}
		
		if (!is_flip_ok(eh)) continue;
		auto vh1 = e.vh1();
		auto vh2 = e.vh2();
		MeshKernel::iGameVertexHandle vh3(-1), vh4(-1);
		auto adjfhs = mesh.NeighborFh(eh);
		//if (adjfhs.size() != 2 || !mesh.isValid(*adjfhs.begin()) || !mesh.isValid(*(adjfhs.begin()++))) continue;
		//assert(adjfhs.size() == 2 && !mesh.isValid(*adjfhs.begin()) && !mesh.isValid(*(adjfhs.begin()++)));
		for (auto fh : adjfhs) {
			auto vex = (mesh.faces(fh)).getVertexHandle();
			assert(vex.size() > 2);
			int tmp = vex[0] + vex[1] + vex[2] - vh1 - vh2;
			if (vh3 == -1) vh3 = MeshKernel::iGameVertexHandle(tmp);
			else vh4 = MeshKernel::iGameVertexHandle(tmp);
		}
		//if (vh3 == vh4 || vh3 == -1 || vh4 == -1) continue;
		int target_degree1 = (mesh.isOnBoundary(vh1)) ? 4 : 6;
		int target_degree2 = (mesh.isOnBoundary(vh2)) ? 4 : 6;
		int target_degree3 = (mesh.isOnBoundary(vh3)) ? 4 : 6;
		int target_degree4 = (mesh.isOnBoundary(vh4)) ? 4 : 6;
		int diff1 = mesh.NeighborVh(vh1).size() - target_degree1;
		int diff2 = mesh.NeighborVh(vh2).size() - target_degree2;
		int diff3 = mesh.NeighborVh(vh3).size() - target_degree3;
		int diff4 = mesh.NeighborVh(vh4).size() - target_degree4;
		int old_diff = diff1 * diff1 + diff2 * diff2 + diff3 * diff3 + diff4 * diff4;
		diff1--, diff2--;
		diff3++, diff4++;
		int new_diff = diff1 * diff1 + diff2 * diff2 + diff3 * diff3 + diff4 * diff4;
		if (new_diff < old_diff) {
			flip_cnt++;
			flip(eh);
		}
	}
	if (flip_cnt != 0) convergence_flag = false;
	printf("flip %d edges\n", flip_cnt);

}

void Remeshing::tangentialRelaxation(int iter) {
	for (int i = 0; i < iter; ++i) {
		std::unordered_map<int, Eigen::Vector3d> centroid;// centroid of face
		std::unordered_map<int, Eigen::Vector3d> avg;// average of vertex's 1-ring face centroid
		std::unordered_map<int, Eigen::Matrix3d> N_face;// prepare ofr error-aware vertex smoothing
		mesh.genAllVerticesNormal();
		for (auto& fp : mesh.allfaces()) {
			Eigen::Vector3d tmp = Eigen::Vector3d::Zero();
			auto vhs = fp.second.getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = mesh.vertices(vh);
				tmp += Eigen::Vector3d(v.x(), v.y(), v.z());
			}
			tmp /= vhs.size();
			centroid[fp.first] = tmp;
		}
		for (auto& vp : mesh.allvertices()) {
			if (mesh.NeighborVh(vp.first).empty() || mesh.isOnBoundary(vp.first)) continue;
			Eigen::Vector3d tmp = Eigen::Vector3d::Zero();
			auto fhs = mesh.NeighborFh(vp.first);
			for (auto& fh : fhs) {
				tmp += centroid[fh];
			}
			tmp /= fhs.size();
			avg[vp.first] = tmp;

		}
		for (auto& vp : mesh.allvertices()) {
			if (mesh.NeighborVh(vp.first).empty()) continue;
			if (mesh.isOnBoundary(vp.first)) continue;
			auto vh = vp.first;
			auto& v = mesh.vertices(vh);
			Eigen::Vector3d N(v.getNormalX(), v.getNormalY(), v.getNormalZ());// the normal of vertex
			Eigen::Vector3d p(v.x(), v.y(), v.z());
			Eigen::Vector3d u = avg[vh] + N * (N.transpose() * (p - avg[vh]));// the update vector u
			v.setPosition(u[0], u[1], u[2]);
		}

		// After each smoothing iteration, we project the vertices in M0 to the original input surface M
		// to maintain a high approximation fidelity.
		std::vector<MeshKernel::iGameVertexHandle> vhs;
		vhs.reserve(mesh.VertexSize());
		for (auto& vp : mesh.allvertices()) {
			vhs.push_back(vp.first);
		}

#pragma omp parallel for
		for (int i = 0; i < vhs.size(); ++i) {// ab_tree is more faster
			if (!mesh.isValid(vhs[i]) || mesh.NeighborVh(vhs[i]).empty()) continue;
			if (mesh.isOnBoundary(vhs[i])) continue;
			MeshKernel::iGameVertex& v = mesh.vertices(vhs[i]);
			Vector3f v_pos(v.x(), v.y(), v.z()), nv;
			ab_tree->findNearstPoint(v_pos, nv);
			v.setPosition(nv[0], nv[1], nv[2]);
		}

	}

}

void Remeshing::projectToSurface() {

	std::vector<MeshKernel::iGameVertexHandle> vhs;
	vhs.reserve(mesh.VertexSize());
	for (auto& vp : mesh.allvertices()) {
		vhs.push_back(vp.first);
	}

#pragma omp parallel for
	for (int i = 0; i < vhs.size(); ++i) {// ab_tree is more faster
		if (!mesh.isValid(vhs[i]) || mesh.NeighborVh(vhs[i]).empty()) continue;
		if (mesh.isOnBoundary(vhs[i])) continue;
		MeshKernel::iGameVertex& v = mesh.vertices(vhs[i]);
		Vector3f v_pos(v.x(), v.y(), v.z()), nv;
		ab_tree->findNearstPoint(v_pos, nv);
		v.setPosition(nv[0], nv[1], nv[2]);
	}

}

void Remeshing::splitLongEdge() {
	int split_cnt = 0;
	auto oldEhs = mesh.alledges();
	for (auto& ep : oldEhs) {
		if (!mesh.isValid(ep.first)) continue;
		auto& e = mesh.edges(ep.first);
		auto vh1 = ep.second.vh1();
		auto vh2 = ep.second.vh2();
		auto& v1 = mesh.vertices(vh1);
		auto& v2 = mesh.vertices(vh2);
		double len = (v1 - v2).norm();
		if (len < upper_length) continue;
		MeshKernel::iGameVertex new_v = (v1 + v2) / 2;
		auto adjfhs = mesh.NeighborFh(ep.first);
		auto new_vh = mesh.AddVertex(new_v);
		for (auto& fh : adjfhs) {// adjacent triangles are bisected
			auto newFace = mesh.faces(fh).getVertexHandle();
			if (newFace[0] == vh1) newFace[0] = new_vh;
			else if (newFace[1] == vh1) newFace[1] = new_vh;
			else if (newFace[2] == vh1) newFace[2] = new_vh;
			mesh.AddFace(newFace);
			newFace = mesh.faces(fh).getVertexHandle();
			if (newFace[0] == vh2) newFace[0] = new_vh;
			else if (newFace[1] == vh2) newFace[1] = new_vh;
			else if (newFace[2] == vh2) newFace[2] = new_vh;
			mesh.AddFace(newFace);
		}
		mesh.DeleteEdge(ep.first);
		split_cnt++;
		

	}

	if (split_cnt != 0) convergence_flag = false;
	printf("split %d edges\n", split_cnt);

}

void Remeshing::collapseShortEdge() {

	int collapse_cnt = 0;
	auto oldEhs = mesh.alledges();
	for (auto& ep : oldEhs) {
		if (!mesh.isValid(ep.first)) continue;
		auto& e = mesh.edges(ep.first);
		auto vh1 = ep.second.vh1();
		auto vh2 = ep.second.vh2();
		auto& v1 = mesh.vertices(vh1);
		auto& v2 = mesh.vertices(vh2);
		double len = (v1 - v2).norm();
		if (len > lower_length) continue;
		auto vh_keep = vh1;
		auto vh_remove = vh2;
		auto& vex_keep = mesh.vertices(vh_keep);
		auto& vex_remove = mesh.vertices(vh_remove);
		// preserve boundary
		if (mesh.isOnBoundary(vh_keep)) {
			// do nothing
		} else if (mesh.isOnBoundary(vh_remove)) {
			vex_keep = vex_remove;
		} else {
			vex_keep = (vex_keep + vex_remove) / 2;
		}

		std::vector<std::vector<MeshKernel::iGameVertexHandle>> new_triangles;
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

		collapse_cnt++;

	}

	if (collapse_cnt != 0) convergence_flag = false;
	printf("collapse %d edges\n", collapse_cnt);

}

bool Remeshing::is_flip_ok(MeshKernel::iGameEdgeHandle eh) {
	if (!mesh.isValid(eh) || mesh.isOnBoundary(eh)) return false;
	auto& e = mesh.edges(eh);
	auto vh1 = e.vh1();
	auto vh2 = e.vh2();
	MeshKernel::iGameVertexHandle vh3(-1), vh4(-1);
	auto adjfhs = mesh.NeighborFh(eh);
	if (adjfhs.size() != 2 || !mesh.isValid(*adjfhs.begin()) || !mesh.isValid(*(adjfhs.begin()++))) 
		return false;

	std::vector<MeshKernel::iGameVertex> face_normal;
	for (auto fh : adjfhs) {
		auto vex = (mesh.faces(fh)).getVertexHandle();
		assert(vex.size() > 2);
		int tmp = vex[0] + vex[1] + vex[2] - vh1 - vh2;
		if (vh3 == -1) vh3 = MeshKernel::iGameVertexHandle(tmp);
		else vh4 = MeshKernel::iGameVertexHandle(tmp);
	}
	if (vh3 == vh4 || vh3 == -1 || vh4 == -1) 
		return false;

	return true;
}

double Remeshing::get_square_difference(MeshKernel::iGameEdgeHandle eh) {
	if (!mesh.isValid(eh)) return Double_MAX;
	auto& e = mesh.edges(eh);
	auto& vh1 = e.vh1();
	auto& vh2 = e.vh2();
	auto adjfhs = mesh.NeighborFh(eh);
	MeshKernel::iGameVertexHandle vh3(-1), vh4(-1);
	for (auto& fh : adjfhs) {
		auto vhs = (mesh.faces(fh)).getVertexHandle();
		if (vh3 == -1) {
			vh3 = MeshKernel::iGameVertexHandle(vhs[0] + vhs[1] + vhs[2] - vh1 - vh2);
		} else {
			vh4 = MeshKernel::iGameVertexHandle(vhs[0] + vhs[1] + vhs[2] - vh1 - vh2);
		}
	}
	if (vh3 == vh4 || vh3 == -1 || vh4 == -1) return Double_MAX;
	auto& v1 = mesh.vertices(vh1);
	auto& v2 = mesh.vertices(vh2);
	auto& v3 = mesh.vertices(vh3);
	auto& v4 = mesh.vertices(vh4);
	auto vec13 = (v3 - v1).normalized();
	auto vec14 = (v4 - v1).normalized();
	auto vec34 = (v4 - v3).normalized();
	auto vec23 = (v3 - v2).normalized();
	auto vec24 = (v4 - v2).normalized();
	double angle314 = std::acos(vec13 * vec14) * rad2angle;// angle314
	double angle134 = std::acos((vec13 * -1) * vec34) * rad2angle;// angle134
	double angle341 = std::acos(vec34 * vec14) * rad2angle;// angle341
	double angle324 = std::acos(vec23 * vec24) * rad2angle;// angle324
	double angle234 = std::acos((vec23 * -1) * vec34)  * rad2angle;// angle234
	double angle243 = std::acos(vec24 * vec34) * rad2angle;// angle243

	double diff1 = angle314 - 60, diff2 = angle134 - 60, diff3 = angle341 - 60,
		diff4 = angle324 - 60, diff5 = angle234 - 60, diff6 = angle243 - 60;

	return (diff1 * diff1 + diff2 * diff2 + diff3 * diff3 + diff4 * diff4 + diff5 * diff5 + diff6 * diff6);

}

void Remeshing::flip(MeshKernel::iGameEdgeHandle eh) {
	if (!mesh.isValid(eh)) return;
	auto& e = mesh.edges(eh);
	auto vh1 = e.vh1();
	auto vh2 = e.vh2();
	MeshKernel::iGameVertexHandle vh3(-1), vh4(-1);
	auto adjfhs = mesh.NeighborFh(eh);
	//if (adjfhs.size() != 2 || !mesh.isValid(*adjfhs.begin()) || !mesh.isValid(*(adjfhs.begin()++))) continue;
	//assert(adjfhs.size() == 2 && !mesh.isValid(*adjfhs.begin()) && !mesh.isValid(*(adjfhs.begin()++)));
	for (auto fh : adjfhs) {
		auto vex = (mesh.faces(fh)).getVertexHandle();
		assert(vex.size() > 2);
		int tmp = vex[0] + vex[1] + vex[2] - vh1 - vh2;
		if (vh3 == -1) vh3 = MeshKernel::iGameVertexHandle(tmp);
		else vh4 = MeshKernel::iGameVertexHandle(tmp);
	}
	for (auto& fh : adjfhs) {
		auto vhs = (mesh.faces(fh)).getVertexHandle();
		int i;
		for (i = 0; i < 3; ++i) {
			int j = (i + 1) % 3;
			if ((vhs[i] == vh1 && vhs[j] == vh2) || (vhs[i] == vh2 && vhs[j] == vh1)) {
				int tmp = vhs[0] + vhs[1] + vhs[2] - vh1 - vh2;
				if (tmp == vh3) vhs[i] = vh4;
				else vhs[i] = vh3;
				assert(vhs[0] != vhs[1] && vhs[1] != vhs[2] && vhs[2] != vhs[0]);
				break;
			}
		}
		assert(i < 4);
		auto new_fh = mesh.AddFace(vhs);
	}
	edge_norm2.erase(eh);
	mesh.DeleteEdge(eh);
	// printf("flip edge %d success\n", eh);
}