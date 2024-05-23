#include "FunctionTester.h"

void FunctionTester::test_kdtree6d(MeshKernel::SurfaceMesh& trimesh, std::vector<int>& fhs) {

	trimesh.genAllVerticesNormal();

	KDTree6d kdtree6d(trimesh);

	double diagonal_length = (trimesh.BBoxMax - trimesh.BBoxMin).norm();
	double weight_n = diagonal_length * 0.2;

	vector<Vector6d> vertices;
	vertices.reserve(trimesh.vsize());
	for (auto& vp : trimesh.allvertices()) {
		auto& v = trimesh.vertices(vp.first);
		const auto& N = v.getNormal();
		vertices.emplace_back(Vector6d(v.x(), v.y(), v.z(), N[0] * weight_n, N[1] * weight_n, N[2] * weight_n));
	}

	for (auto& fp : trimesh.allfaces()) {
		if (std::rand() % 100 == 0) {
			Vex face_center = trimesh.getFaceCenter(fp.first);
			Vex face_normal = trimesh.getFaceNormal(fp.first);
			Vector6d vec6(face_center.x(), face_center.y(), face_center.z(),
				face_normal.x(), face_normal.y(), face_normal.z());
			auto nn = kdtree6d.nearest(vec6);
			Vector6d nearest = nn.nearest;
			//std::cout << "Input vec6: " << vec6 << std::endl;
			//std::cout << "Nearest vec6: " << nearest << std::endl;
			//std::cout << "Dist = " << nn.dist << ", FH = " << nn.fh << std::endl;
			//std::cout << "Input FH = " << fp.first << ", nearest fh = " << nn.fh << std::endl;
			double dist_min = 99999999;
			Vector6d nearest_trace;
			for (auto& v6d : vertices) {
				double dist = (v6d - vec6).Length();
				if (dist < dist_min) {
					dist_min = dist;
					nearest_trace = v6d;
				}
			}
			std::cout << "\nKD Tree result: " << nearest << std::endl;
			std::cout << "Actually result: " << nearest_trace << std::endl;
			/*Vector3d diff_vec(nearest[0] - nearest_trace[0], nearest[1] - nearest_trace[1], nearest[2] - nearest_trace[2]);
			double diff_length = diff_vec.Length();
			std::cout << "The diff_vec length between KD Tree result and trace result: " << diff_length << std::endl;*/
			fhs.emplace_back(fp.first);
			fhs.emplace_back(nn.fh);
		}
	}

	//for (auto& vp : trimesh.allvertices()) {
	//	if (std::rand() % 100 == 0) {
	//		Vector6d vec6(vp.second.x(), vp.second.y(), vp.second.z(), 
	//			vp.second.getNormalX(), vp.second.getNormalY(), vp.second.getNormalZ());
	//		auto nn = kdtree6d.nearest(vec6);
	//		/*Vector6d nearest = nn.nearest;
	//		std::cout << "Input vec6: " << vec6 << std::endl;
	//		std::cout << "Nearest vec6: " << nearest << std::endl;
	//		std::cout << "Dist = " << nn.dist << std::endl;*/
	//		auto& face = trimesh.faces(nn.fh);
	//		const auto& f_vhs = face.getVertexHandle();
	//		if (std::find(f_vhs.begin(), f_vhs.end(), vp.first) == f_vhs.end()) {
	//			std::cout << "VH is not included in FH\n";
	//		} else {
	//			std::cout << "VH is included in FH\n";
	//		}
	//	}
	//}


}

void FunctionTester::test_pca(SkeletalMesh& sklmesh, std::vector<std::vector<Vex>>& pca_res) {

	pca_res.clear();

	unordered_map<EH, bool> curve_skel_ehs;
	double edge_len_avg = 0;

	for (auto& ep : sklmesh.alledges()) {
		if (sklmesh.NeighborFh(ep.first).empty()) {
			curve_skel_ehs[ep.first] = true;// 曲线边的度为0
			edge_len_avg += sklmesh.getLength(ep.first);
		}
	}
	if (curve_skel_ehs.empty()) return;
	edge_len_avg /= curve_skel_ehs.size();

	for (auto& vp : sklmesh.allvertices()) {
		VH vh = vp.first;
		auto& v = vp.second;
		std::vector<Eigen::Vector3d> vectors;
		const auto& adjehs = sklmesh.NeighborEh(vh);
		for (auto& adjeh : adjehs) {
			if (curve_skel_ehs.count(adjeh)) {
				auto& adje = sklmesh.edges(adjeh);
				VH adjvh(adje.vh1() + adje.vh2() - vh);
				auto& adjv = sklmesh.vertices(adjvh);
				Vec vec = (adjv - v).normalized();
				vectors.emplace_back(Eigen::Vector3d(vec.x(), vec.y(), vec.z()));
			}
		}
		if (vectors.size() > 2) {
			std::vector<Eigen::Vector3d> result;
			Math_PCA::get_principal_components(vectors, result);
			for (auto& dir : result) {
				Vec dir_(dir.x(), dir.y(), dir.z());
				dir_.normalize();
				std::vector<Vex> dd = {
					v, v + dir_ * edge_len_avg * 0.1
				};
				pca_res.emplace_back(dd);
				break;
			}
		}
	}

}

void FunctionTester::flip_same_orientation_faces(MeshKernel::SurfaceMesh& mesh, FH fh_src, vector<int>& selected_fhs) {

	auto fhs = mesh.allfaces();

	for (auto& fp : fhs) {
		auto vhs = mesh.faces(fp.first).getVertexHandle();
		std::reverse(vhs.begin(), vhs.end());
		mesh.AddFace(vhs);
	}

	/*if (!mesh.isValid(fh_src)) return;

	mesh.genAllFacesNormal();

	std::unordered_set<FH> fhs;
	std::queue<FH> que;
	que.emplace(fh_src);
	fhs.insert(fh_src);

	while (!que.empty()) {

		auto fh = que.front();
		que.pop();
		auto normal0 = mesh.faces(fh).getNormal();

		for (auto& adjfh : mesh.NeighborFh(fh)) {
			if (fhs.count(adjfh)) continue;
			auto normal1 = mesh.faces(adjfh).getNormal();
			if (normal0.dot(normal1) > 0) {
				que.emplace(adjfh);
				fhs.insert(adjfh);
			}
		}


	}

	std::cout << "fhs size = " << fhs.size() << std::endl;

	selected_fhs.clear();
	selected_fhs = vector<int>(fhs.begin(), fhs.end());*/

	/*for (auto& fh : fhs) {
		auto vhs = mesh.faces(fh).getVertexHandle();
		mesh.DeleteFace(fh);
		std::reverse(vhs.begin(), vhs.end());
		mesh.AddFace(vhs);
	}

	mesh.updateAllHandles();*/

}