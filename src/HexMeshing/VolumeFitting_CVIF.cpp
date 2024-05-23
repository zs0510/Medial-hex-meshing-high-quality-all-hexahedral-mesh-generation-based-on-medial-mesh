#include "VolumeFitting_CVIF.h"

void VolumeFitting_CVIF::volume_fitting(int subd_times, 
	bool project_to_surface_flag, 
	int iter_times,
	unordered_map<VH, Vex> _corner_features_hex,
	unordered_map<VH, int> _line_features_hex,
	vector<vector<EH>> _line_features_ref) {

//#define USE_KD_TREE_6D
//#define OUTPUT_CTRL_MESH
//#define OUTPUT_SUBD_MESH
//#define OUTPUT_KDTREE_MESH
//#define OUTPUT_ITER_INFO

#ifdef USE_KD_TREE_6D
	use_kdtree6d_flag = true;
#endif // USE_KD_TREE_6D

	use_refmesh_face = true;
	corner_features_hex = _corner_features_hex;
	line_features_hex = _line_features_hex;
	line_features_ref = _line_features_ref;

	time_t beg = clock();
	ctrl_subdivision_iter = subd_times;

	vector<Vector3f> aabb_vertices;
	for (auto& fp : refmesh.allfaces()) {
		const auto& vhs = fp.second.getVertexHandle();
		for (int i = 2; i < vhs.size(); ++i) {
			vector<VH> trivhs = { vhs[0], vhs[i - 1], vhs[i] };
			for (auto& vh : trivhs) {
				auto& v = refmesh.vertices(vh);
				aabb_vertices.emplace_back(Vector3f(v.x(), v.y(), v.z()));
			}
		}
	}
	aabb_tree = new AABB_Tree(aabb_vertices);

	init_refmesh_data();

	// 1. extract control quad mesh
	init_control_mesh();

#ifdef OUTPUT_CTRL_MESH
	MeshKernel::IO io_ctrl;
	io_ctrl.WriteObjFile(control_mesh, "C:/My Files/Graphics/model_data/!test_data/CVIF/control_mesh.obj");
#endif // OUTPUT_CTRL_MESH


	//use_kdtree6d_flag = true;// 第五次迭代后切换至六维的 KD Tree
	int it = 0;
	for ( ; it < iter_times; ++it) {

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
		//calc_movement_face_center();
#ifdef OUTPUT_ITER_INFO
		printf("Iter %d: Calculate movement success.\n", it);
#endif // OUTPUT_ITER_INFO
		
		
		// 5. update control quad mesh
		update_control_mesh();
		
#ifdef OUTPUT_ITER_INFO
		printf("Iter %d: Update control mesh success.\n", it);
#endif // OUTPUT_ITER_INFO
		
		/*string filename = "C:/My Files/Graphics/model_data/CVIF/CVIF_Result" + to_string(it + 1) + ".off";
		MeshKernel::IO io;
		io.WriteOffFile(control_mesh, filename);*/

		feature_recovery();

		// 6. smoothing if iter is not last once
		if (it + 1 != iter_times && !is_converged) {// add candition [&& it < 5] is original paper
			if (project_to_surface_flag) {
				project_to_surface();
			}
			SurfaceMesh_Smoothing app_s(control_mesh);
			app_s.Laplacian_smoothing(1);
			feature_recovery();
		}

		
		
	}

	// 7. update surface of hexmesh and diffuse the movement to inner vertices by using level Laplacian operation
	update_hexmesh();

	time_t end = clock();

	std::cout << "\n[CVIF]: Iterations = " << (it - 1) << ", cost " << int(end - beg) << "ms.\n" << std::endl;

}

void VolumeFitting_CVIF::feature_recovery() {

	// 特征恢复
	for (auto& line_feature_hex : line_features_hex) {
		VH vh = hex2ctrl[line_feature_hex.first];// 四边形网格(control_mesh)上的顶点
		auto& ehs_ref = line_features_ref[line_feature_hex.second];
		Vex& v = control_mesh.vertices(vh);
		double dist_nearest = DOUBLE_MAX;
		Vex pos_nearest;
		for (auto& eid : ehs_ref) {
			EH eh(eid);
			Vex p;
			double dist = dist_point_line_segment(v, refmesh.vertices(refmesh.edges(eh).vh1()), refmesh.vertices(refmesh.edges(eh).vh2()), p);
			if (dist < dist_nearest) {
				dist_nearest = dist;
				pos_nearest = p;
			}
		}
		v = pos_nearest;
	}

	for (auto& corner_feature_hex : corner_features_hex) {
		VH vh = hex2ctrl[corner_feature_hex.first];// 四边形网格(control_mesh)上的顶点
		Vex& v = control_mesh.vertices(vh);
		v = corner_feature_hex.second;
	}

}

void VolumeFitting_CVIF::init_refmesh_data() {

	if (use_refmesh_face) {

		refmesh_points.reserve(refmesh.fsize());
		refmesh_normals.reserve(refmesh.fsize());
		for (auto& fp : refmesh.allfaces()) {
			auto v = refmesh.getFaceCenter(fp.first);
			refmesh_points.emplace_back(v);
			refmesh_normals.emplace_back(Vector3d(fp.second.getNormalX(), fp.second.getNormalY(), fp.second.getNormalZ()));
		}
		std::cout << "\n[CVIF]: We sample on faces of reference mesh.\n" << std::endl;

	} else {

		refmesh_points.reserve(refmesh.vsize());
		refmesh_normals.reserve(refmesh.vsize());
		for (auto& vp : refmesh.allvertices()) {
			auto& v = vp.second;
			refmesh_points.emplace_back(v);
			refmesh_normals.emplace_back(Vector3d(v.getNormalX(), v.getNormalY(), v.getNormalZ()));
		}
		std::cout << "\n[CVIF]: We sample on vertices of reference mesh.\n" << std::endl;

	}

	

}

void VolumeFitting_CVIF::init_control_mesh() {

	control_mesh.destory();

	
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



}

void VolumeFitting_CVIF::init_kdtree() {

	// 1. generate triangle mesh using subdivision quad mesh
	trisub2quadsub.clear();
	MeshKernel::SurfaceMesh trimesh;
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

void VolumeFitting_CVIF::calc_movement() {

	int ref_vcnt = refmesh_points.size();
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
		//VH vh_ref(vi);
		//auto& v_ref = refmesh.vertices(vh_ref);
		auto& v_ref = refmesh_points[vi];
		auto& normal_ref = refmesh_normals[vi];

		VH vh_sub;// KD 树输出的 点

		if (use_kdtree6d_flag) {
			Vector6d input_vec6(v_ref.x(), v_ref.y(), v_ref.z(), normal_ref.x() * weight_n, normal_ref.y() * weight_n, normal_ref.z() * weight_n);
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

void VolumeFitting_CVIF::calc_movement_face_center() {

	int ref_fcnt = refmesh.fsize();
	nearest_deltas.resize(ref_fcnt);
	nearest_subvhs.resize(ref_fcnt);

	double weight_n;
	if (use_kdtree6d_flag) {
		weight_n = kdtree6d->diagonal_length * kdtree6d->weight_normal;
	}

#ifdef USE_KD_TREE_6D
	weight_n = kdtree6d->diagonal_length * kdtree6d->weight_normal;
#endif // USE_KD_TREE_6D



#pragma omp parallel for
	for (int fi = 0; fi < ref_fcnt; ++fi) {// 并行以加快速度
		FH fh_ref(fi);
		auto v_ref = refmesh.getFaceCenter(fh_ref);

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
			auto fh_sub = nearest.fh;
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

		
		Vec delta = v_ref - subdivision_mesh.vertices(vh_sub);// δ = p−vcl,γ
		nearest_deltas[fi] = delta;
		nearest_subvhs[fi] = vh_sub;
	}

	double new_rms_error = 0;
	movements.clear();
	for (int fi = 0; fi < ref_fcnt; ++fi) {// 此处不可并行
		VH vh_sub = nearest_subvhs[fi];
		for (auto& wp : control_weights[vh_sub]) {// 所有与细分曲面顶点有关的控制网格都会受到影响
			movements[wp.first].emplace_back(wp.second, nearest_deltas[fi]);
		}
		new_rms_error += nearest_deltas[fi].norm2();
	}
	new_rms_error = std::sqrt(new_rms_error);
	new_rms_error /= ref_fcnt;

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

void VolumeFitting_CVIF::update_control_mesh() {

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

void VolumeFitting_CVIF::update_hexmesh() {

	int ctrl_vcnt = control_mesh.vsize();
#pragma omp parallel for
	for (int vi = 0; vi < ctrl_vcnt; ++vi) {
		VH vh_ctrl(vi);
		auto& v_ctrl = control_mesh.vertices(vh_ctrl);
		auto& v_hex = hexmesh.vertices(ctrl2hex[vh_ctrl]);
		v_hex = v_ctrl;
	}

	HexMesh_Smoothing app(hexmesh);
	app.Laplacian_Level(15, false);// 只光滑内部, 不光滑边界

}

void Coarse_CVIF::cvif_coarse() {

	auto mesh_sub = hexmesh;
	Subdivision_Linear app_sub(mesh_sub);
	app_sub.subdivision(2, false);// 保证了顶点 id 不变

	VolumeFitting_CVIF app_cvif(mesh_sub, trimesh);
	app_cvif.volume_fitting(1);

	for (auto& vp : hexmesh.allvertices()) {
		auto& v = hexmesh.vertices(vp.first);
		auto& v_cvif = mesh_sub.vertices(vp.first);
		v = v_cvif;
	}

}

void VolumeFitting_CVIF::project_to_surface() {

	int ctrl_vcnt = control_mesh.vsize();
	int iters = 3;
	for (int it = 0; it < iters; ++it) {
	
#pragma omp parallel for
		for (int vi = 0; vi < ctrl_vcnt; ++vi) {
			VH vh(vi);
			auto& v = control_mesh.vertices(vh);
			Vector3f pos(v.x(), v.y(), v.z()), npos;
			aabb_tree->findNearstPoint(pos, npos);
			v.setPosition(npos[0], npos[1], npos[2]);
		}
		if (it + 1 != iters) {
			SurfaceMesh_Smoothing app_s(control_mesh);
			app_s.Laplacian_smoothing(1);
		}
		
	}



}

double VolumeFitting_CVIF::dist_point_line_segment(iGameVertex& v, iGameVertex& v0, iGameVertex& v1, iGameVertex& nearest_vertex) {
	iGameVertex vec1 = v - v0;
	iGameVertex vec2 = v1 - v0;
	double t = vec2 * vec2;
	iGameVertex min_v = v0;
	if (t > 1E-6F) {
		t = vec1 * vec2 / t;
		if (t > 1.f) {
			min_v = v1;
			vec1 = v - min_v;
		} else if (t > 0.f) {
			min_v = v0 + vec2 * t;
			vec1 = v - min_v;
		}
	}
	nearest_vertex = min_v;
	return vec1.norm();
}