#include "CVIF_Bezier_Unsmoothed.h"

CVIF_Bezier_Unsmoothed::CVIF_Bezier_Unsmoothed(MeshKernel::VolumeMesh& hexmesh, MeshKernel::SurfaceMesh& target_mesh): refmesh(target_mesh) {

	generate_control_mesh(hexmesh);

}

void CVIF_Bezier_Unsmoothed::cvif_bezier(vector<vector<vector<Vector3d>>>& MutiBezierPatch) {


	time_t beg = clock();

	int it = 0;
	ctrl_subdivision_iter = 2;

	for (; it < 100; ++it) {

		//if (it == 5) 

		if (is_converged) {
			std::cout << "[CVIF]: Result is converaged!\n";
			break;
		}

		// 2. generate subdivision quad mesh
		QuadMeshSubdivision_CC app(control_mesh);
		app.subdivision(ctrl_subdivision_iter);
		app.get_subdivision_mesh(subdivision_mesh);
		app.get_control_weights(control_weights);
#ifdef OUTPUT_ITER_INFO
		printf("Iter %d: Generate subdivision mesh success.\n", it);
#endif // OUTPUT_ITER_INFO



#ifdef OUTPUT_SUBD_MESH
		MeshKernel::IO io_sub;
		io_sub.WriteObjFile(subdivision_mesh, "C:/My Files/Graphics/model_data/!test_data/CVIF/subdivision_mesh.obj");
#undef OUTPUT_SUBD_MESH
#endif // OUTPUT_SUBD_MESH

		// 3. generate kd_tree using subdivision quad mesh
		init_kdtree();
#ifdef OUTPUT_ITER_INFO
		printf("Iter %d: Initialize KD-Tree success.\n", it);
#endif // OUTPUT_ITER_INFO

		// 4. calculate movement for each vertex on control quad mesh
		calc_movement();
#ifdef OUTPUT_ITER_INFO
		printf("Iter %d: Calculate movement success.\n", it);
#endif // OUTPUT_ITER_INFO


		// 5. update control quad mesh
		update_control_mesh();
#ifdef OUTPUT_ITER_INFO
		printf("Iter %d: Update control mesh success.\n", it);
#endif // OUTPUT_ITER_INFO


		// 6. smoothing if iter is not last once
		if (it + 1 != 100) {// add candition [&& it < 5] is original paper
			SurfaceMesh_Smoothing app_s(control_mesh);
			app_s.Laplacian_smoothing(1);
			//#ifdef OUTPUT_ITER_INFO
			//			printf("Iter %d: Smoothing control mesh success.\n", it);
			//#endif // OUTPUT_ITER_INFO
		}

	}

	extract_bezier_surface(MutiBezierPatch);

	time_t end = clock();

	std::cout << "\n[CVIF]: Iterations = " << (it - 1) << ", cost " << int(end - beg) << "ms.\n" << std::endl;

}

void CVIF_Bezier_Unsmoothed::init_kdtree() {

	// 1. generate triangle mesh using subdivision quad mesh
	MeshKernel::SurfaceMesh trimesh;
	trisub2quadsub.clear();
	int vcnt = subdivision_mesh.vsize();
	for (int i = 0; i < vcnt; ++i) {
		VH vh(i);
		auto& v = subdivision_mesh.vertices(vh);
		VH tri_vh = trimesh.AddVertex(v);
		assert(tri_vh == i);
	}
	for (auto& fp : subdivision_mesh.allfaces()) {
		const auto& quad_vhs = fp.second.getVertexHandle();
		FH fh1 = trimesh.AddFace({ quad_vhs[0], quad_vhs[1], quad_vhs[2] });
		FH fh2 = trimesh.AddFace({ quad_vhs[0], quad_vhs[2], quad_vhs[3] });
		trisub2quadsub[fh1] = fp.first;
		trisub2quadsub[fh2] = fp.first;
	}

	delete kdtree6d;
	kdtree6d = nullptr;
	delete kd_tree;
	kd_tree = nullptr;

	if (use_kdtree6d_flag) {
		kdtree6d = new KDTree6d(trimesh);
	} else {
		kd_tree = new MeshKernel::iGameKdTree(trimesh);
	}
	//#ifdef USE_KD_TREE_6D
	//	kdtree6d = new KDTree6d(trimesh);
	//#else
	//	kd_tree = new MeshKernel::iGameKdTree(trimesh);
	//#endif // USE_KD_TREE_6D

#ifdef OUTPUT_KDTREE_MESH
	MeshKernel::IO io_trimesh;
	io_trimesh.WriteObjFile(trimesh, "C:/My Files/Graphics/model_data/!test_data/CVIF/triangle_mesh.obj");
#undef OUTPUT_KDTREE_MESH
#endif // OUTPUT_KDTREE_MESH

	subdivision_tri_mesh = trimesh;

}

void CVIF_Bezier_Unsmoothed::calc_movement() {

	int ref_vcnt = refmesh.vsize();
	nearest_deltas.resize(ref_vcnt);
	nearest_subvhs.resize(ref_vcnt);

	double weight_n;
	if (use_kdtree6d_flag) {
		weight_n = kdtree6d->diagonal_length * kdtree6d->weight_normal;
	}

#ifdef USE_KD_TREE_6D
	weight_n = kdtree6d->diagonal_length * kdtree6d->weight_normal;
#endif // USE_KD_TREE_6D



#pragma omp parallel for
	for (int vi = 0; vi < ref_vcnt; ++vi) {// 并行以加快速度
		VH vh_ref(vi);
		auto& v_ref = refmesh.vertices(vh_ref);

		VH vh_sub;// KD 树输出的 点

		if (use_kdtree6d_flag) {
			Vector6d input_vec6(v_ref.x(), v_ref.y(), v_ref.z(), v_ref.getNormalX() * weight_n, v_ref.getNormalY() * weight_n, v_ref.getNormalZ() * weight_n);
			auto nn = kdtree6d->nearest(input_vec6);
			auto& face_sub = subdivision_tri_mesh.faces(nn.fh);
			Vector3d output_vec3(nn.nearest[0], nn.nearest[1], nn.nearest[2]);
			double dist_min = DOUBLE_MAX;
			for (auto& vh_ : face_sub.getVertexHandle()) {
				auto& v_sub_ = subdivision_tri_mesh.vertices(vh_);
				Vector3d vec_sub_(v_sub_.x(), v_sub_.y(), v_sub_.z());
				double dis = (output_vec3 - vec_sub_).Length();
				if (dis < dist_min) {
					dist_min = dis;
					vh_sub = vh_;
				}
			}
		} else {
			auto nearest = kd_tree->nearest(v_ref);
			auto fh_sub = trisub2quadsub[nearest.fh];
			auto& face_sub = subdivision_mesh.faces(fh_sub);
			double dist_min = DOUBLE_MAX;
			for (auto& vh_ : face_sub.getVertexHandle()) {
				auto& v_sub_ = subdivision_mesh.vertices(vh_);
				double dis = (v_ref - v_sub_).norm();
				if (dis < dist_min) {
					dist_min = dis;
					vh_sub = vh_;// 细分曲面上离当前点最近的一个点
				}
			}

		}

		//#ifdef USE_KD_TREE_6D
		//		// KD Tree 6d
		//		Vector6d input_vec6(v_ref.x(), v_ref.y(), v_ref.z(), v_ref.getNormalX() * weight_n, v_ref.getNormalY() * weight_n, v_ref.getNormalZ() * weight_n);
		//		auto nn = kdtree6d->nearest(input_vec6);
		//		auto& face_sub = subdivision_tri_mesh.faces(nn.fh);
		//		Vector3d output_vec3(nn.nearest[0], nn.nearest[1], nn.nearest[2]);
		//		double dist_min = DOUBLE_MAX;
		//		for (auto& vh_ : face_sub.getVertexHandle()) {
		//			auto& v_sub_ = subdivision_tri_mesh.vertices(vh_);
		//			Vector3d vec_sub_(v_sub_.x(), v_sub_.y(), v_sub_.z());
		//			double dis = (output_vec3 - vec_sub_).Length();
		//			if (dis < dist_min) {
		//				dist_min = dis;
		//				vh_sub = vh_;
		//			}
		//		}
		//
		//
		//#else
		//		// KD Tree 3d
		//		auto nearest = kd_tree->nearest(v_ref);
		//		auto& face_sub = nearest.face;
		//		double dist_min = DOUBLE_MAX;
		//		for (auto& vh_ : face_sub.getVertexHandle()) {
		//			auto& v_sub_ = subdivision_mesh.vertices(vh_);
		//			double dis = (v_ref - v_sub_).norm();
		//			if (dis < dist_min) {
		//				dist_min = dis;
		//				vh_sub = vh_;
		//			}
		//		}
		//#endif // USE_KD_TREE_6D


		Vec delta = v_ref - subdivision_mesh.vertices(vh_sub);// δ = p−vcl,γ
		nearest_deltas[vi] = delta;
		nearest_subvhs[vi] = vh_sub;
	}

	double new_rms_error = 0;
	movements.clear();
	for (int vi = 0; vi < ref_vcnt; ++vi) {// 此处不可并行
		VH vh_sub = nearest_subvhs[vi];
		for (auto& wp : control_weights[vh_sub]) {// 所有与细分曲面顶点有关的控制网格都会受到影响
			movements[wp.first].emplace_back(wp.second, nearest_deltas[vi]);
		}
		new_rms_error += nearest_deltas[vi].norm2();
	}
	new_rms_error = std::sqrt(new_rms_error);
	new_rms_error /= ref_vcnt;

	if (rms_error > 0) {
		double ratio = new_rms_error / rms_error;
		double absdiff = std::abs(1.0 - ratio);
		if (absdiff < epsilon_error) {
			is_converged = true;
			//std::cout << "[CVIF]: Result is converaged!\n";
		}
		printf("RMS_ERROR: prev = %.6f, curr = %.6f, abs(diff) = %.6f\n", rms_error, new_rms_error, absdiff);
	}

	rms_error = new_rms_error;

}

void CVIF_Bezier_Unsmoothed::update_control_mesh() {


	int ctrl_vcnt = control_mesh.vsize();
	bool cvif_nan_error = false;
#pragma omp parallel for
	for (int vi = 0; vi < ctrl_vcnt; ++vi) {
		VH vh(vi);
		double weight_sum = 0;
		Vec move_sum(0, 0, 0);
		for (auto& mp : movements[vh]) {
			move_sum += mp.second * mp.first;
			weight_sum += mp.first;
		}
		Vec move_v = move_sum / weight_sum;
		/*if (std::isnan(weight_sum)) {
			std::cerr << "[CVIF Error]: Weight sum is not a number." << std::endl;
			continue;
		}*/
		if (std::isnan(move_v.x()) || std::isnan(move_v.y()) || std::isnan(move_v.z())) {
			cvif_nan_error = true;
			continue;
		}
		auto& v = control_mesh.vertices(vh);
		v += move_v;
	}

	if (cvif_nan_error) {
		std::cerr << "[CVIF Error]: Exist number is not a number." << std::endl;
	}

}

void CVIF_Bezier_Unsmoothed::extract_bezier_surface(vector<vector<vector<Vector3d>>>& MutiBezierPatch) {

	MutiBezierPatch.clear();

	for (auto& single_path : mutil_bezier_patch) {
		vector<vector<Vector3d>> controlPointsSurface;
		for (auto& single_line : single_path) {
			vector<Vector3d> controlPointsLine;
			for (auto& vh : single_line) {
				auto& v = control_mesh.vertices(vh);
				controlPointsLine.push_back(Vector3d(v.x(), v.y(), v.z()));
			}
			controlPointsSurface.emplace_back(controlPointsLine);
		}
		MutiBezierPatch.emplace_back(controlPointsSurface);
	}

}

void CVIF_Bezier_Unsmoothed::generate_control_mesh(MeshKernel::VolumeMesh& hexmesh) {

	control_mesh.destory();

	unordered_map<VH, VH> hex2ctrl;// 输入六面体 VH 与其控制网格 VH 的映射关系
	unordered_map<VH, VH> ctrl2hex;

	for (auto& vp : hexmesh.allvertices()) {
		if (hexmesh.isOnBoundary(vp.first)) {
			VH qvh = control_mesh.AddVertex(vp.second);
			hex2ctrl[vp.first] = qvh;
			ctrl2hex[qvh] = vp.first;
		}
	}

	for (auto& fp : hexmesh.allfaces()) {
		if (hexmesh.isOnBoundary(fp.first)) {
			FH hex_fh = fp.first;
			Vex face_center = hexmesh.getFaceCenter(hex_fh);
			CH ch = *(hexmesh.NeighborCh(fp.first).begin());
			Vec face_normal = hexmesh.getFaceNormal(hex_fh);
			Vex cell_center = hexmesh.getCellCenter(ch);
			Vec dir_out = (face_center - cell_center).normalized();
			auto vhs = fp.second.getVertexHandle();
			if (dir_out.dot(face_normal) < 0) std::reverse(vhs.begin(), vhs.end());
			for (auto& vh : vhs) {
				vh = hex2ctrl[vh];
			}
			control_mesh.AddFace(vhs);
		}
	}

	sub_to_bezier_ctrlmesh();

}

void CVIF_Bezier_Unsmoothed::sub_to_bezier_ctrlmesh() {

	int ecnt = control_mesh.esize();
	int fcnt_before = control_mesh.fsize();

	vector<vector<VH>> edge_vertices;// 每条边生成两个顶点
	edge_vertices.reserve(ecnt);

	// 1. 生成所有边点
	for (int eid = 0; eid < ecnt; ++eid) {

		EH eh(eid);
		auto& edge = control_mesh.edges(eh);
		VH vh1 = edge.vh1(), vh2 = edge.vh2();
		auto& v1 = control_mesh.vertices(vh1);
		auto& v2 = control_mesh.vertices(vh2);

		vector<VH> vh_new;
		Vex v1_new = v1 * 0.666666 + v2 * 0.333333;
		Vex v2_new = v1 * 0.333333 + v2 * 0.666666;
		vh_new.emplace_back(control_mesh.AddVertex(v1_new));
		vh_new.emplace_back(control_mesh.AddVertex(v2_new));

		edge_vertices.emplace_back(vh_new);

	}


	// 2. 生成所有面点并添加面
	auto faces_origin = control_mesh.allfaces();
	auto edges_origin = control_mesh.alledges();
	mutil_bezier_patch.clear();

	for (auto& fp : faces_origin) {

		auto& face = control_mesh.faces(fp.first);
		auto f_vhs = face.getVertexHandle();
		vector<EH> f_ehs = { control_mesh.getEdgeHandle(f_vhs[0], f_vhs[1]), control_mesh.getEdgeHandle(f_vhs[1], f_vhs[2]),
		control_mesh.getEdgeHandle(f_vhs[2], f_vhs[3]), control_mesh.getEdgeHandle(f_vhs[3], f_vhs[0]) };
		for (int i = 0; i < 4; ++i) {
			auto& edge = control_mesh.edges(f_ehs[i]);
			//if (edge.vh1() != f_vhs[i]) {// 确保顺序不错
			//	std::swap(edge_vertices[f_ehs[i]][0], edge_vertices[f_ehs[i]][1]);
			//}
			auto& ev = control_mesh.vertices(edge_vertices[f_ehs[i]][0]);
			auto& v_cur = control_mesh.vertices(f_vhs[i]);
			auto& v_next = control_mesh.vertices(f_vhs[(i + 1) % 4]);
			double dist0 = (v_cur - ev).norm2();
			double dist1 = (v_next - ev).norm2();
			if (dist0 > dist1) {
				std::swap(edge_vertices[f_ehs[i]][0], edge_vertices[f_ehs[i]][1]);
			}
		}
		Vex v0_new = control_mesh.vertices(edge_vertices[f_ehs[0]][0]) * 0.666666 + control_mesh.vertices(edge_vertices[f_ehs[2]][1]) * 0.333333;
		Vex v3_new = control_mesh.vertices(edge_vertices[f_ehs[0]][0]) * 0.333333 + control_mesh.vertices(edge_vertices[f_ehs[2]][1]) * 0.666666;
		Vex v1_new = control_mesh.vertices(edge_vertices[f_ehs[0]][1]) * 0.666666 + control_mesh.vertices(edge_vertices[f_ehs[2]][0]) * 0.333333;
		Vex v2_new = control_mesh.vertices(edge_vertices[f_ehs[0]][1]) * 0.333333 + control_mesh.vertices(edge_vertices[f_ehs[2]][0]) * 0.666666;
		vector<VH> v_vhs = {
			control_mesh.AddVertex(v0_new), control_mesh.AddVertex(v1_new), control_mesh.AddVertex(v2_new), control_mesh.AddVertex(v3_new)
		};
		// 沿着所有顶点添加面
		control_mesh.AddFace({ f_vhs[0], edge_vertices[f_ehs[0]][0], v_vhs[0], edge_vertices[f_ehs[3]][1] });
		control_mesh.AddFace({ f_vhs[1], edge_vertices[f_ehs[1]][0], v_vhs[1], edge_vertices[f_ehs[0]][1] });
		control_mesh.AddFace({ f_vhs[2], edge_vertices[f_ehs[2]][0], v_vhs[2], edge_vertices[f_ehs[1]][1] });
		control_mesh.AddFace({ f_vhs[3], edge_vertices[f_ehs[3]][0], v_vhs[3], edge_vertices[f_ehs[2]][1] });
		// 沿着所有边添加面
		control_mesh.AddFace({ edge_vertices[f_ehs[0]][0], edge_vertices[f_ehs[0]][1], v_vhs[1], v_vhs[0] });
		control_mesh.AddFace({ edge_vertices[f_ehs[1]][0], edge_vertices[f_ehs[1]][1], v_vhs[2], v_vhs[1] });
		control_mesh.AddFace({ edge_vertices[f_ehs[2]][0], edge_vertices[f_ehs[2]][1], v_vhs[3], v_vhs[2] });
		control_mesh.AddFace({ edge_vertices[f_ehs[3]][0], edge_vertices[f_ehs[3]][1], v_vhs[0], v_vhs[3] });

		// 在面中心添加面
		control_mesh.AddFace(v_vhs);

		vector<vector<VH>> single_patch = {
			{ f_vhs[0], edge_vertices[f_ehs[0]][0], edge_vertices[f_ehs[0]][1], f_vhs[1] },
			{ edge_vertices[f_ehs[3]][1], v_vhs[0], v_vhs[1], edge_vertices[f_ehs[1]][0] },
			{ edge_vertices[f_ehs[3]][0], v_vhs[3], v_vhs[2], edge_vertices[f_ehs[1]][1] },
			{ f_vhs[3], edge_vertices[f_ehs[2]][1], edge_vertices[f_ehs[2]][0], f_vhs[2] }
		};

		mutil_bezier_patch.push_back(single_patch);

	}

	/*for (auto& fp : faces_origin) {
		control_mesh.DeleteFace(fp.first);
	}*/

	for (auto& ep : edges_origin) {
		control_mesh.DeleteEdge(ep.first);
	}

	//control_mesh.updateAllHandles();
	int fcnt_after = control_mesh.fsize();

	std::cout << "[CVIF_Bezier]: Bezier generater, faces size " << fcnt_before << " --> " << fcnt_after << std::endl;

}