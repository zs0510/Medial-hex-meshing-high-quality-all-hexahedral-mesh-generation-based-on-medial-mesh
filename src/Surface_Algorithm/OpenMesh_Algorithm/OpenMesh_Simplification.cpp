#include "OpenMesh_Simplification.h"

void iGame_OpenMesh_Simplification::simplification_qem(size_t target_vcnt) {

	auto vcnt = mesh.n_vertices();
	if (target_vcnt >= vcnt) return;

	mesh.add_property(error_mats);
	mesh.add_property(collapse_nodes);
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	init_error_mats();

	for (auto e_it = mesh.edges_sbegin(); e_it != mesh.edges_end(); ++e_it) {
		calc_collapse_cost(e_it);
	}
	std::cout << "Init collapse nodes success, heap size = " << collapse_heap.size() << std::endl;

	size_t iter_times = vcnt - target_vcnt;

	OM_VH vh_remove, vh_keep;
	OM_HEH heh;
	std::time_t time_beg = clock();
	while (iter_times > 0 && !collapse_heap.empty()) {

		CollapseNode* node = collapse_heap.top();
		collapse_heap.pop();
		while (!collapse_heap.empty() && node->valid == false) {
			node = collapse_heap.top();
			collapse_heap.pop();
		}
		//std::cout << "The selected cost is " << node->cost << std::endl;
		if (node->valid == false) break;

		heh = mesh.halfedge_handle(node->eh, 0);
		vh_remove = mesh.from_vertex_handle(heh);
		vh_keep = mesh.to_vertex_handle(heh);

		if (!mesh.is_valid_handle(vh_remove) || !mesh.is_valid_handle(vh_keep)) continue;

		if (mesh.valence(vh_keep) < mesh.valence(vh_remove)) {
			std::swap(vh_remove, vh_keep);
			heh = mesh.opposite_halfedge_handle(heh);
		}

		if (mesh.is_collapse_ok(heh)) {
			
			//std::cout << "HEH #" << heh << " is collapse OK." << std::endl;

			// 清除无效的收缩信息
			for (auto* _node : mesh.property(collapse_nodes, vh_remove)) {
				_node->valid = false;
			}
			mesh.property(collapse_nodes, vh_remove).clear();
			for (auto* _node : mesh.property(collapse_nodes, vh_keep)) {
				_node->valid = false;
			}
			mesh.property(collapse_nodes, vh_keep).clear();
			//std::cout << "Invalidate outdated collapse nodes success." << std::endl;

			// 累加误差矩阵
			mesh.property(error_mats, vh_keep) += mesh.property(error_mats, vh_remove);
			// 设置坍塌位置
			mesh.set_point(vh_keep, node->pos);
			
			/*std::cout << "VH #" << vh_keep << ": " << mesh.valence(vh_keep) << ", VH #" << vh_remove << ": " << mesh.valence(vh_remove) 
				<< "; heh to vertex: VH #" << mesh.to_vertex_handle(heh) << std::endl;*/

			// 边收缩
			mesh.collapse(heh);

			//assert(mesh.is_valid_handle(vh_keep) && "vh_keep should be valid!");
			//std::cout << "Collapse success." << std::endl;

			for (auto ve_it = mesh.ve_begin(vh_keep); ve_it.is_valid(); ++ve_it) {
				calc_collapse_cost(*ve_it);
				//std::cout << "Calculate collapse nodes success." << std::endl;
			}

			iter_times--;

		}
	}
	std::time_t time_end = clock();

	mesh.delete_isolated_vertices();
	mesh.garbage_collection();

	mesh.remove_property(error_mats);
	mesh.remove_property(collapse_nodes);
	mesh.release_face_status();
	mesh.release_edge_status();
	mesh.release_vertex_status();

	std::cout << "Input vertices size = " << vcnt << ", output vertices size = " << mesh.n_vertices() << 
		". cost: " << int(time_end - time_beg) << std::endl;


}

void iGame_OpenMesh_Simplification::init_error_mats() {

	mesh.request_face_normals();
	mesh.update_face_normals();

#pragma omp parallel for
	for (int i = 0; i < mesh.n_vertices(); ++i) {
		auto vh = mesh.vertex_handle(i);
		auto& mat = mesh.property(error_mats, vh);
		mat = Eigen::Matrix4d::Zero();
		double d;
		Eigen::Vector4d p;
		for (auto vf_it = mesh.vf_begin(vh); vf_it.is_valid(); ++vf_it) {
			auto nor = mesh.normal(*vf_it).data();
			d = mesh.normal(*vf_it).dot(mesh.point(vh));
			p = Eigen::Vector4d(nor[0], nor[1], nor[2], -d);
			mat += (p * p.transpose());
		}
	}

	/*for (auto v_it = mesh.vertices_sbegin(); v_it != mesh.vertices_end(); ++v_it) {
		auto& vh = *v_it;
		
	}*/

}

void iGame_OpenMesh_Simplification::calc_collapse_cost(OM_EH eh) {

	if (!mesh.is_valid_handle(eh)) return;

	double cost = (double)std::numeric_limits<double>::max();
	OM_HEH heh = mesh.halfedge_handle(eh, 0);
	const auto& vh0 = mesh.from_vertex_handle(heh);
	const auto& vh1 = mesh.to_vertex_handle(heh);
	Eigen::Vector4d X;// 最佳新顶点位置
	Eigen::Matrix4d Q = mesh.property(error_mats, vh0) + mesh.property(error_mats, vh1);
	Eigen::Matrix4d A;
	A << Q(0, 0), Q(0, 1), Q(0, 2), Q(0, 3),
		Q(0, 1), Q(1, 1), Q(1, 2), Q(1, 3),
		Q(0, 2), Q(1, 2), Q(2, 2), Q(2, 3),
		0.f, 0.f, 0.f, 1.f;

	if (A.determinant()) {
		// A 可逆，可通过解方程获取最佳位置
		Eigen::Vector4d b(0.f, 0.f, 0.f, 1.f);
		X = A.colPivHouseholderQr().solve(b);
		cost = X.transpose() * Q * X;
	} else {
		// A 不可逆，在 线段 v1-v2 上寻找局部最佳位置
		auto& v0 = mesh.point(vh0);
		auto& v1 = mesh.point(vh1);
		for (double i = 0.f; i <= 1.f; i += 0.05f) {
			auto tmp_pos = (v0 * i + v1 * (1 - i)).data();
			Eigen::Vector4d vec4(tmp_pos[0], tmp_pos[1], tmp_pos[2], 1.f);
			double tmp_cost = vec4.transpose() * Q * vec4;
			if (tmp_cost < cost) {
				cost = tmp_cost;
				X = vec4;
			}
		}
	}

	CollapseNode* node = new CollapseNode(cost, OM_Vex(X[0], X[1], X[2]), eh);

	collapse_heap.push(node);
	mesh.property(collapse_nodes, vh0).push_back(node);
	mesh.property(collapse_nodes, vh1).push_back(node);

}