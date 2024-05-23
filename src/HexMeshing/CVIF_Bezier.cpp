#include "CVIF_Bezier.h"

CVIF_Bezier::CVIF_Bezier(MeshKernel::VolumeMesh& hexmesh, MeshKernel::SurfaceMesh& target_mesh): refmesh(target_mesh) {

	generate_surface_mesh(hexmesh);

}

void CVIF_Bezier::cvif_bezier(vector<vector<vector<Vector3d>>>& MutiBezierPatch) {


	time_t beg = clock();

	int it = 0;

	for (; it < 100; ++it) {

		if (is_converged) {
			std::cout << "[CVIF]: Result is converaged!\n";
			break;
		}

		// generate subdivision quad mesh

		generate_bezier_mesh();// 生成 Bezier 控制网格

		QuadMeshSubdivision_CC app(bezier_mesh);// 将 Bezier 网格细分
		app.subdivision(1);
		app.get_subdivision_mesh(subdivision_mesh);
		app.get_control_weights(subdivision_to_bezier_weights);

		// generate kd_tree using subdivision quad mesh
		init_kdtree();

		// calculate movement for each vertex on control quad mesh
		calc_movement();

		// update control quad mesh
		update_surface_mesh();

		// smoothing if iter is not last once
		if (!is_converged) {// add candition [&& it < 5] is original paper
			SurfaceMesh_Smoothing app_s(surface_mesh);
			app_s.Laplacian_smoothing(1);
		}

	}

	generate_bezier_mesh(true);
	extract_bezier_surface(MutiBezierPatch);

	time_t end = clock();

	std::cout << "\n[CVIF_Bezier]: Iterations = " << (it - 1) << ", cost " << int(end - beg) << "ms.\n" << std::endl;

}

void CVIF_Bezier::init_kdtree() {

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

void CVIF_Bezier::calc_movement() {

	int ref_vcnt = refmesh.vsize();
	nearest_deltas.resize(ref_vcnt);
	nearest_subvhs.resize(ref_vcnt);

	double weight_n;
	if (use_kdtree6d_flag) {
		weight_n = kdtree6d->diagonal_length * kdtree6d->weight_normal;
	}

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

		Vec delta = v_ref - subdivision_mesh.vertices(vh_sub);// δ = p−vcl,γ
		nearest_deltas[vi] = delta;
		nearest_subvhs[vi] = vh_sub;
	}

	// 计算Bezier曲面的移动向量
	bezier_movements.clear();
	double new_rms_error = 0;
	for (int vi = 0; vi < ref_vcnt; ++vi) {// 此处不可并行
		VH vh_sub = nearest_subvhs[vi];
		for (auto& wp : subdivision_to_bezier_weights[vh_sub]) {// 所有与细分曲面顶点有关的控制网格都会受到影响
			bezier_movements[wp.first].emplace_back(wp.second, nearest_deltas[vi]);
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
		//printf("RMS_ERROR: prev = %.6f, curr = %.6f, abs(diff) = %.6f\n", rms_error, new_rms_error, absdiff);
		printf("RMS_ERROR: abs(diff) = %.6f\n", absdiff);
	}
	rms_error = new_rms_error;

	int bezier_vcnt = bezier_mesh.vsize();
	vector<Vec> bezier_move(bezier_vcnt, Vec(0, 0, 0));
#pragma omp parallel for
	for (int vi = 0; vi < bezier_vcnt; ++vi) {
		VH vh(vi);
		double weight_sum = 0;
		Vec move_sum(0, 0, 0);
		for (auto& mp : bezier_movements[vh]) {
			move_sum += mp.second * mp.first;
			weight_sum += mp.first;
		}
		Vec move_v = move_sum / weight_sum;
		if (std::isnan(move_v.x()) || std::isnan(move_v.y()) || std::isnan(move_v.z())) {
			continue;
		}
		bezier_move[vi] = move_v;
	}

	// 计算控制曲面的移动向量
	control_movements.clear();
	
	for (int vi = 0; vi < bezier_vcnt; ++vi) {
		VH vh(vi);
		for (auto& wp : bezier_to_control_weights[vh]) {
			control_movements[wp.first].emplace_back(wp.second, bezier_move[vi]);
		}
	}

}

void CVIF_Bezier::calc_movement_face_center() {

	int ref_fcnt = refmesh.fsize();
	nearest_deltas.resize(ref_fcnt);
	nearest_subvhs.resize(ref_fcnt);

	double weight_n;
	if (use_kdtree6d_flag) {
		weight_n = kdtree6d->diagonal_length * kdtree6d->weight_normal;
	}
	
#pragma omp parallel for
	for (int fi = 0; fi < ref_fcnt; ++fi) {// 并行以加快速度
		FH fh_ref(fi);
		Vex v_ref = refmesh.getFaceCenter(fh_ref);

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

	// 计算Bezier曲面的移动向量
	bezier_movements.clear();
	double new_rms_error = 0;
	for (int fi = 0; fi < ref_fcnt; ++fi) {// 此处不可并行
		VH vh_sub = nearest_subvhs[fi];
		for (auto& wp : subdivision_to_bezier_weights[vh_sub]) {// 所有与细分曲面顶点有关的控制网格都会受到影响
			bezier_movements[wp.first].emplace_back(wp.second, nearest_deltas[fi]);
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
		//printf("RMS_ERROR: prev = %.6f, curr = %.6f, abs(diff) = %.6f\n", rms_error, new_rms_error, absdiff);
		printf("RMS_ERROR: abs(diff) = %.6f\n", absdiff);
	}
	rms_error = new_rms_error;

	int bezier_vcnt = bezier_mesh.vsize();
	vector<Vec> bezier_move(bezier_vcnt, Vec(0, 0, 0));
#pragma omp parallel for
	for (int vi = 0; vi < bezier_vcnt; ++vi) {
		VH vh(vi);
		double weight_sum = 0;
		Vec move_sum(0, 0, 0);
		for (auto& mp : bezier_movements[vh]) {
			move_sum += mp.second * mp.first;
			weight_sum += mp.first;
		}
		Vec move_v = move_sum / weight_sum;
		if (std::isnan(move_v.x()) || std::isnan(move_v.y()) || std::isnan(move_v.z())) {
			continue;
		}
		bezier_move[vi] = move_v;
	}

	// 计算控制曲面的移动向量
	control_movements.clear();

	for (int vi = 0; vi < bezier_vcnt; ++vi) {
		VH vh(vi);
		for (auto& wp : bezier_to_control_weights[vh]) {
			control_movements[wp.first].emplace_back(wp.second, bezier_move[vi]);
		}
	}

}

void CVIF_Bezier::update_surface_mesh() {


	int ctrl_vcnt = surface_mesh.vsize();
	bool cvif_nan_error = false;
#pragma omp parallel for
	for (int vi = 0; vi < ctrl_vcnt; ++vi) {
		VH vh(vi);
		double weight_sum = 0;
		Vec move_sum(0, 0, 0);
		for (auto& mp : control_movements[vh]) {
			move_sum += mp.second * mp.first;
			weight_sum += mp.first;
		}
		Vec move_v = move_sum / weight_sum;
		if (std::isnan(move_v.x()) || std::isnan(move_v.y()) || std::isnan(move_v.z())) {
			cvif_nan_error = true;
			continue;
		}
		auto& v = surface_mesh.vertices(vh);
		v += move_v;
	}

	if (cvif_nan_error) {
		std::cerr << "[CVIF Error]: Exist number is not a number." << std::endl;
	}

}

void CVIF_Bezier::extract_bezier_surface(vector<vector<vector<Vector3d>>>& MutiBezierPatch) {

	MutiBezierPatch.clear();

	for (auto& single_path : mutil_bezier_patch) {
		vector<vector<Vector3d>> controlPointsSurface;
		for (auto& single_line : single_path) {
			vector<Vector3d> controlPointsLine;
			for (auto& vh : single_line) {
				auto& v = bezier_mesh.vertices(vh);
				controlPointsLine.push_back(Vector3d(v.x(), v.y(), v.z()));
			}
			controlPointsSurface.emplace_back(controlPointsLine);
		}
		MutiBezierPatch.emplace_back(controlPointsSurface);
	}

}

void CVIF_Bezier::generate_surface_mesh(MeshKernel::VolumeMesh& hexmesh) {

	surface_mesh.destory();

	unordered_map<VH, VH> hex2ctrl;// 输入六面体 VH 与其控制网格 VH 的映射关系
	unordered_map<VH, VH> ctrl2hex;

	for (auto& vp : hexmesh.allvertices()) {
		if (hexmesh.isOnBoundary(vp.first)) {
			VH qvh = surface_mesh.AddVertex(vp.second);
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
			surface_mesh.AddFace(vhs);
		}
	}

}

void CVIF_Bezier::generate_bezier_mesh(bool extract_patch) {
	// 根据当前控制顶点重新生成 Bezier 控制网格
	bezier_mesh.destory();
	bezier_to_control_weights.clear();
	v_to_corner_point.clear();
	fv_to_inter_point.clear();
	ev_to_edge_point.clear();

	if (extract_patch) {
		mutil_bezier_patch.clear();
	}
	
	for (auto& vp : surface_mesh.allvertices()) {
		calc_corner_point(vp.first);
	}

	for (auto& ep : surface_mesh.alledges()) {
		calc_edge_point(ep.first, ep.second.vh1());
		calc_edge_point(ep.first, ep.second.vh2());
	}

	for (auto& fp : surface_mesh.allfaces()) {
		const auto& f_vhs = fp.second.getVertexHandle();
		calc_inter_point(fp.first, f_vhs[0]);
		calc_inter_point(fp.first, f_vhs[1]);
		calc_inter_point(fp.first, f_vhs[2]);
		calc_inter_point(fp.first, f_vhs[3]);
	}

	for (auto& fp : surface_mesh.allfaces()) {

		const auto& f_vhs = fp.second.getVertexHandle();
		vector<EH> f_ehs = { surface_mesh.getEdgeHandle(f_vhs[0], f_vhs[1]), surface_mesh.getEdgeHandle(f_vhs[1], f_vhs[2]), 
		surface_mesh.getEdgeHandle(f_vhs[2], f_vhs[3]), surface_mesh.getEdgeHandle(f_vhs[3], f_vhs[0]), };
		// Bezier 曲面的数据
		vector<vector<VH>> _e_vhs = { ev_to_edge_point[f_ehs[0]], ev_to_edge_point[f_ehs[1]], ev_to_edge_point[f_ehs[2]], ev_to_edge_point[f_ehs[3]] };
		vector<VH> _v_vhs = { v_to_corner_point[f_vhs[0]], v_to_corner_point[f_vhs[1]], v_to_corner_point[f_vhs[2]], v_to_corner_point[f_vhs[3]] };
		vector<VH> _f_vhs = fv_to_inter_point[fp.first];

		for (int i = 0; i < 4; ++i) {
			auto& v1 = surface_mesh.vertices(f_vhs[i]);
			auto& _v1 = bezier_mesh.vertices(_e_vhs[i][0]);
			auto& _v2 = bezier_mesh.vertices(_e_vhs[i][1]);
			double norm2_1 = (v1 - _v1).norm2();
			double norm2_2 = (v1 - _v2).norm2();
			if (norm2_2 < norm2_1) {
				std::swap(_e_vhs[i][0], _e_vhs[i][1]);
			}
		}

		// 沿着所有顶点添加面
		bezier_mesh.AddFace({ _v_vhs[0], _e_vhs[0][0], _f_vhs[0], _e_vhs[3][1] });
		bezier_mesh.AddFace({ _v_vhs[1], _e_vhs[1][0], _f_vhs[1], _e_vhs[0][1] });
		bezier_mesh.AddFace({ _v_vhs[2], _e_vhs[2][0], _f_vhs[2], _e_vhs[1][1] });
		bezier_mesh.AddFace({ _v_vhs[3], _e_vhs[3][0], _f_vhs[3], _e_vhs[2][1] });
		// 沿着所有边添加面
		bezier_mesh.AddFace({ _e_vhs[0][0], _e_vhs[0][1], _f_vhs[1], _f_vhs[0] });
		bezier_mesh.AddFace({ _e_vhs[1][0], _e_vhs[1][1], _f_vhs[2], _f_vhs[1] });
		bezier_mesh.AddFace({ _e_vhs[2][0], _e_vhs[2][1], _f_vhs[3], _f_vhs[2] });
		bezier_mesh.AddFace({ _e_vhs[3][0], _e_vhs[3][1], _f_vhs[0], _f_vhs[3] });

		// 在面中心添加面
		bezier_mesh.AddFace(_f_vhs);

		if (extract_patch) {
			// 记录 Bezier 曲面控制顶点
			vector<vector<VH>> single_patch = {
				{ _v_vhs[0], _e_vhs[0][0], _e_vhs[0][1], _v_vhs[1] },
				{ _e_vhs[3][1], _f_vhs[0], _f_vhs[1], _e_vhs[1][0] },
				{ _e_vhs[3][0], _f_vhs[3], _f_vhs[2], _e_vhs[1][1] },
				{ _v_vhs[3], _e_vhs[2][1], _e_vhs[2][0], _v_vhs[2] }
			};

			mutil_bezier_patch.push_back(single_patch);
		}

	}



}

void CVIF_Bezier::calc_corner_point(VH vh) {

	unordered_map<VH, double> weight;
	auto& v = surface_mesh.vertices(vh);
	const auto& adjvhs = surface_mesh.NeighborVh(vh);
	int n = adjvhs.size();// 顶点的度
	double weight_sum = 0.0;

	Vex v_new = v * n * n;// 计算新顶点
	weight[vh] += n * n;// 保存权值
	weight_sum += n * n;// 计算总权值

	for (auto& adjvh : adjvhs) {
		auto& adjv = surface_mesh.vertices(adjvh);
		v_new += adjv * 4.0;
		weight[adjvh] += 4.0;
		weight_sum += 4.0;
	}

	const auto& adjfhs = surface_mesh.NeighborFh(vh);

	for (auto& adjfh : adjfhs) {
		auto& face = surface_mesh.faces(adjfh);
		const auto& f_vhs = face.getVertexHandle();
		VH vh_oppo(-1);
		for (auto& f_vh : f_vhs) {
			if (f_vh == vh || adjvhs.count(f_vh)) continue;
			vh_oppo = f_vh;
			break;
		}
		assert(vh_oppo != -1);
		auto& v_oppo = surface_mesh.vertices(vh_oppo);
		v_new += v_oppo;
		weight[vh_oppo] += 1.0;
		weight_sum += 1.0;
	}

	v_new /= weight_sum;

	VH vh_new = bezier_mesh.AddVertex(v_new);

	v_to_corner_point[vh] = vh_new;

	for (auto& wp : weight) {
		bezier_to_control_weights[vh_new][wp.first] = wp.second / weight_sum;// 归一化
	}

}

void CVIF_Bezier::calc_inter_point(FH fh, VH vh) {

	unordered_map<VH, double> weight;
	auto& v = surface_mesh.vertices(vh);
	const auto& adjvhs = surface_mesh.NeighborVh(vh);
	int n = adjvhs.size();// 顶点的度
	double weight_sum = 0.0;

	Vex v_new = v * n;
	weight[vh] += n;
	weight_sum += n;

	auto& face = surface_mesh.faces(fh);
	const auto& f_vhs = face.getVertexHandle();
	for (auto& f_vh : f_vhs) {
		if (f_vh == vh) continue;
		double w = 2.0;
		if (!adjvhs.count(f_vh)) {// 对立顶点
			w = 1.0;
		}
		auto& f_v = surface_mesh.vertices(f_vh);
		v_new += f_v * w;
		weight[f_vh] += w;
		weight_sum += w;
	}

	v_new /= weight_sum;

	VH vh_new = bezier_mesh.AddVertex(v_new);

	fv_to_inter_point[fh].emplace_back(vh_new);

	for (auto& wp : weight) {
		bezier_to_control_weights[vh_new][wp.first] = wp.second / weight_sum;// 归一化
	}

}

void CVIF_Bezier::calc_edge_point(EH eh, VH vh) {

	unordered_set<VH> adjvhs_e;
	const auto& adjfhs_e = surface_mesh.NeighborFh(eh);
	for (auto& adjfh : adjfhs_e) {
		auto& face = surface_mesh.faces(adjfh);
		const auto& f_vhs = face.getVertexHandle();
		adjvhs_e.insert(f_vhs.begin(), f_vhs.end());
	}

	unordered_map<VH, double> weight;
	auto& v = surface_mesh.vertices(vh);
	const auto& adjvhs = surface_mesh.NeighborVh(vh);
	int n = adjvhs.size();// 顶点的度
	double weight_sum = 0.0;

	Vex v_new = v * 2.0 * n;// 计算新顶点
	weight[vh] += 2.0 * n;// 保存权值
	weight_sum += 2.0 * n;// 计算总权值

	auto& edge = surface_mesh.edges(eh);
	VH vh_op(edge.vh1() + edge.vh2() - vh);
	auto& v_op = surface_mesh.vertices(vh_op);

	v_new += v_op * 4.0;
	weight[vh_op] += 4.0;
	weight_sum += 4.0;

	const auto& adjvhs_op = surface_mesh.NeighborVh(vh_op);

	for (auto& adjvh : adjvhs) {
		if (adjvh == vh_op || !adjvhs_e.count(adjvh)) continue;
		auto& adjv = surface_mesh.vertices(adjvh);
		v_new += adjv * 2.0;
		weight[adjvh] += 2.0;
		weight_sum += 2.0;
	}

	for (auto& adjvh : adjvhs_op) {
		if (adjvh == vh || !adjvhs_e.count(adjvh)) continue;
		auto& adjv = surface_mesh.vertices(adjvh);
		v_new += adjv * 1.0;
		weight[adjvh] += 1.0;
		weight_sum += 1.0;
	}

	v_new /= weight_sum;

	VH vh_new = bezier_mesh.AddVertex(v_new);

	ev_to_edge_point[eh].emplace_back(vh_new);

	for (auto& wp : weight) {
		bezier_to_control_weights[vh_new][wp.first] = wp.second / weight_sum;// 归一化
	}


}