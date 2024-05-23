#include "CAD_Reconstruction.h"


void CAD_Reconstruction::get_bezier_control_mesh(MeshKernel::VolumeMesh& _hexmesh, vector<vector<vector<Vector3d>>>& MutiBezierPatch) {

	generate_surface_mesh(_hexmesh);
	get_bezier_control_mesh(MutiBezierPatch);

}

void CAD_Reconstruction::get_bezier_control_mesh(MeshKernel::SurfaceMesh& _quadmesh, vector<vector<vector<Vector3d>>>& MutiBezierPatch) {

	surface_mesh = _quadmesh;
	get_bezier_control_mesh(MutiBezierPatch);

}

void  CAD_Reconstruction::get_bezier_control_mesh(MeshKernel::SurfaceMesh& _quadmesh, vector<vector<vector<Vector3d>>>& MutiBezierPatch,
	MeshKernel::SurfaceMesh& _refmesh, unordered_set<EH>& feature_ehs, unordered_set<VH>& feature_vhs) {

	surface_mesh = _quadmesh;
	get_bezier_control_mesh(MutiBezierPatch, _refmesh, feature_ehs, feature_vhs);

}

void CAD_Reconstruction::get_bezier_control_mesh(vector<vector<vector<Vector3d>>>& MutiBezierPatch) {

	// 根据当前控制顶点重新生成 Bezier 控制网格
	bezier_mesh.destory();
	v_to_corner_point.clear();
	fv_to_inter_point.clear();
	ev_to_edge_point.clear();

	for (auto& vp : surface_mesh.allvertices()) {
		calc_corner_point(vp.first);
	}

	for (auto& ep : surface_mesh.alledges()) {
		calc_edge_point(ep.first, ep.second.vh1());
		calc_edge_point(ep.first, ep.second.vh2());
	}

	for (auto& fp : surface_mesh.allfaces()) {
		const auto& f_vhs = fp.second.getVertexHandle();
		/*auto vp = surface_mesh.vertices(f_vhs[0]);
		cout << "edges_1v: " << vp.x() << "," << vp.y() << "," << vp.z() << "," << endl;
		vp = surface_mesh.vertices(f_vhs[1]);
		cout << "edges_2v: " << vp.x() << "," << vp.y() << "," << vp.z() << "," << endl;
		vp = surface_mesh.vertices(f_vhs[2]);
		cout << "edges_3v: " << vp.x() << "," << vp.y() << "," << vp.z() << "," << endl;
		vp = surface_mesh.vertices(f_vhs[3]);
		cout << "edges_4v: " << vp.x() << "," << vp.y() << "," << vp.z() << "," << endl;*/

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

		vector<vector<VH>> single_patch = {
				{ _v_vhs[0], _e_vhs[0][0], _e_vhs[0][1], _v_vhs[1] },
				{ _e_vhs[3][1], _f_vhs[0], _f_vhs[1], _e_vhs[1][0] },
				{ _e_vhs[3][0], _f_vhs[3], _f_vhs[2], _e_vhs[1][1] },
				{ _v_vhs[3], _e_vhs[2][1], _e_vhs[2][0], _v_vhs[2] }
		};

		mutil_bezier_patch.push_back(single_patch);

	}

	for (auto& single_path : mutil_bezier_patch) {
		//std::cout << "Bezier patch: \n";
		vector<vector<Vector3d>> controlPointsSurface;
		for (auto& single_line : single_path) {
			vector<Vector3d> controlPointsLine;
			//std::cout << "\tSingle line: ";
			for (auto& vh : single_line) {
				auto& v = bezier_mesh.vertices(vh);
				controlPointsLine.push_back(Vector3d(v.x(), v.y(), v.z()));
				//std::cout << "(" << v.x() << "," << v.y() << "," << v.z() << "), ";
			}
			//std::cout << "\n";
			controlPointsSurface.emplace_back(controlPointsLine);
		}
		//std::cout << "\n";
		MutiBezierPatch.emplace_back(controlPointsSurface);
	}

}


void CAD_Reconstruction::get_bezier_control_mesh(vector<vector<vector<Vector3d>>>& MutiBezierPatch,
	MeshKernel::SurfaceMesh& _refmesh, unordered_set<EH>& feature_ehs, unordered_set<VH>& feature_vhs) {


	// 根据当前控制顶点重新生成 Bezier 控制网格
	bezier_mesh.destory();
	v_to_corner_point.clear();
	fv_to_inter_point.clear();
	ev_to_edge_point.clear();

	for (auto& vp : surface_mesh.allvertices()) {
		if (!feature_vhs.count(vp.first)) {
			calc_corner_point(vp.first);

		} else {// 特征顶点不变
			auto& v = _refmesh.vertices(vp.first);
			VH vh_new = bezier_mesh.AddVertex(v);
			v_to_corner_point[vp.first] = vh_new;

		}
		
	}
	std::cout << "Calculate corner points success.\n";

	for (auto& ep : surface_mesh.alledges()) {
		calc_edge_point(ep.first, ep.second.vh1());
		calc_edge_point(ep.first, ep.second.vh2());

		if (feature_ehs.count(ep.first)) {// 将特征边的边点插值在边上
			auto& v1 = _refmesh.vertices(ep.second.vh1());
			auto& v2 = _refmesh.vertices(ep.second.vh2());
			for (auto& vh : ev_to_edge_point[ep.first]) {
				auto& v = bezier_mesh.vertices(vh);
				auto nearest_v = v;
				dist_point_line_segment(v, v1, v2, nearest_v);
				v = nearest_v;
			}

		}
	}
	std::cout << "Calculate edge points success.\n";

	for (auto& fp : surface_mesh.allfaces()) {
		const auto& f_vhs = fp.second.getVertexHandle();
		calc_inter_point(fp.first, f_vhs[0]);
		calc_inter_point(fp.first, f_vhs[1]);
		calc_inter_point(fp.first, f_vhs[2]);
		calc_inter_point(fp.first, f_vhs[3]);
	}
	std::cout << "Calculate face points success.\n";

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

		vector<vector<VH>> single_patch = {
				{ _v_vhs[0], _e_vhs[0][0], _e_vhs[0][1], _v_vhs[1] },
				{ _e_vhs[3][1], _f_vhs[0], _f_vhs[1], _e_vhs[1][0] },
				{ _e_vhs[3][0], _f_vhs[3], _f_vhs[2], _e_vhs[1][1] },
				{ _v_vhs[3], _e_vhs[2][1], _e_vhs[2][0], _v_vhs[2] }
		};

		mutil_bezier_patch.push_back(single_patch);

	}

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


void CAD_Reconstruction::generate_surface_mesh(MeshKernel::VolumeMesh& hexmesh) {

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

void CAD_Reconstruction::calc_corner_point(VH vh) {

	unordered_map<VH, double> weight;
	auto& v = surface_mesh.vertices(vh);
	const auto& adjvhs = surface_mesh.NeighborVh(vh);
	int n = adjvhs.size();// 顶点的度
	double weight_sum = 0.0;

	Vex v_new = v * n * n;// 计算新顶点
	weight[vh] += n * n;// 保存权值
	weight_sum += n * n;// 计算总权值
	//triplets.emplace_back(vh, vh, n * n);

	for (auto& adjvh : adjvhs) {
		auto& adjv = surface_mesh.vertices(adjvh);
		v_new += adjv * 4.0;
		weight[adjvh] += 4.0;
		weight_sum += 4.0;
		//triplets.emplace_back(vh, adjvh, 4);
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
		//triplets.emplace_back(vh, vh_oppo, 1);
	}

	//corner_weights_sum[vh] = weight_sum;

	v_new /= weight_sum;

	VH vh_new = bezier_mesh.AddVertex(v_new);

	v_to_corner_point[vh] = vh_new;



	/*std::cout << "V surface: (" << v.x() << "," << v.y() << "," << v.z() << ") \t---> "
		<< "V bezier: (" << v_new.x() << "," << v_new.y() << "," << v_new.z() << ")"
		<< std::endl;*/

}

void CAD_Reconstruction::calc_inter_point(FH fh, VH vh) {

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

	//for (auto& wp : weight) {
	//	bezier_to_control_weights[vh_new][wp.first] = wp.second / weight_sum;// 归一化
	//}

}

void CAD_Reconstruction::calc_edge_point(EH eh, VH vh) {

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

	//for (auto& wp : weight) {
	//	bezier_to_control_weights[vh_new][wp.first] = wp.second / weight_sum;// 归一化
	//}


}

void CAD_Reconstruction::solve_aux_surface_mesh(MeshKernel::VolumeMesh& hexmesh, MeshKernel::SurfaceMesh& auxmesh) {

	generate_surface_mesh(hexmesh);
	solve_aux_surface_mesh(surface_mesh, auxmesh);

}

void CAD_Reconstruction::solve_aux_surface_mesh(MeshKernel::SurfaceMesh& quadmesh, MeshKernel::SurfaceMesh& auxmesh) {

	if (surface_mesh.vsize() == 0) surface_mesh = quadmesh;// 否则 surface_mesh 已被六面体网格初始化, 无需再次赋值
	auxmesh = quadmesh;
	int vcnt = surface_mesh.vsize();
	triplets.clear();
	corner_weights_sum.clear();
	corner_weights_sum.resize(vcnt);

	for (auto& vp : surface_mesh.allvertices()) {
		calc_corner_point_aux(vp.first);
	}

	// 建立稀疏系数矩阵
	Eigen::SparseMatrix<double> A;
	A.resize(vcnt, vcnt);
	A.setFromTriplets(triplets.begin(), triplets.end());

	// LU 预处理
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);

	Eigen::VectorXd bx, by, bz;
	bx.setZero(vcnt), by.setZero(vcnt), bz.setZero(vcnt);
	for (int vi = 0; vi < vcnt; ++vi) {
		VH vh(vi);
		auto& v = surface_mesh.vertices(vh);
		bx[vi] = v.x() * corner_weights_sum[vi];
		by[vi] = v.y() * corner_weights_sum[vi];
		bz[vi] = v.z() * corner_weights_sum[vi];
	}

	// 解方程
	Eigen::VectorXd xx, xy, xz;
	xx = solver.solve(bx);
	xy = solver.solve(by);
	xz = solver.solve(bz);

	for (auto& vp : auxmesh.allvertices()) {
		auto& v = vp.second;
		v.setPosition(xx[vp.first], xy[vp.first], xz[vp.first]);
	}
	
	//// 检查控制顶点是否插值
	//surface_mesh = auxmesh;
	//for (auto& vp : surface_mesh.allvertices()) {
	//	Vex v_new = calc_corner_point_aux(vp.first);
	//	double dist2 = (v_new - quadmesh.vertices(vp.first)).norm2();
	//	std::cout << "The norm2 is " << dist2 << " ." << std::endl;
	//}

}

Vex CAD_Reconstruction::calc_corner_point_aux(VH vh) {

	unordered_map<VH, double> weight;
	auto& v = surface_mesh.vertices(vh);
	const auto& adjvhs = surface_mesh.NeighborVh(vh);
	int n = adjvhs.size();// 顶点的度
	double weight_sum = 0.0;

	Vex v_new = v * n * n;// 计算新顶点
	weight[vh] += n * n;// 保存权值
	weight_sum += n * n;// 计算总权值
	triplets.emplace_back(vh, vh, n * n);

	for (auto& adjvh : adjvhs) {
		auto& adjv = surface_mesh.vertices(adjvh);
		v_new += adjv * 4.0;
		weight[adjvh] += 4.0;
		weight_sum += 4.0;
		triplets.emplace_back(vh, adjvh, 4);
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
		triplets.emplace_back(vh, vh_oppo, 1);
	}

	corner_weights_sum[vh] = weight_sum;
	v_new /= weight_sum;

	return v_new;

}

double CAD_Reconstruction::dist_point_line_segment(MeshKernel::iGameVertex& v, MeshKernel::iGameVertex& v0,
	MeshKernel::iGameVertex& v1, MeshKernel::iGameVertex& nearest_vertex) {
	auto vec1 = v - v0;
	auto vec2 = v1 - v0;
	double t = vec2 * vec2;
	auto min_v = v0;
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