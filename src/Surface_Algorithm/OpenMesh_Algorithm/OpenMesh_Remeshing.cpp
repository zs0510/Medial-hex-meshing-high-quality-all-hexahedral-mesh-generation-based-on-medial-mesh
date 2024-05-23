#include "OpenMesh_Remeshing.h"

using namespace std;

void iGame_OpenMesh_Remeshing::init_aabb_tree() {

	std::vector<Vector3f> vertices;
	vertices.reserve(mesh.n_faces() * 3);
	for (auto f_it = mesh.faces_sbegin(); f_it != mesh.faces_end(); ++f_it) {
		auto fvs = f_it->vertices();
		for (auto fv_it = fvs.begin(); fv_it != fvs.end(); ++fv_it) {
			auto pos = mesh.point(*fv_it).data();
			vertices.push_back(Vector3f(pos[0], pos[1], pos[2]));
		}
	}
	ab_tree = new AABB_Tree(vertices);

}

void iGame_OpenMesh_Remeshing::init_target_length(double ratio) {

	double sum_area = 0.f, count = 0;
	for (auto f_it = mesh.faces_sbegin(); f_it != mesh.faces_end(); ++f_it) {
		sum_area += mesh.calc_face_area(*f_it);
		count++;
	}
	double target_length = std::sqrt(sum_area / (1.5 * mesh.n_vertices()));// target number of vertices of output mesh
	target_length *= (2 / pow(3, 0.25f));
	target_length *= ratio;
	lower_length = target_length * 4 / 5;
	upper_length = target_length * 4 / 3;
	std::cout << "Target length: " << target_length << std::endl;
}

void iGame_OpenMesh_Remeshing::remeshing_uniform(double ratio, size_t iter_times) {

	init_target_length(ratio);
	std::time_t before, after;
	mesh.request_face_status();
	mesh.request_edge_status();
	mesh.request_vertex_status();
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	std::cout << "Remeshing begin.\n";
	for (size_t it = 0; it < iter_times; ++it) {

		before = clock();
		split_long_edge();
		after = clock();
		printf("iter %d: split long edge success, cost: %dms, edges size: %d\n", int(it), int(after - before), (int)mesh.n_edges());

		before = clock();
		collapse_short_edge();
		after = clock();
		printf("iter %d: collapse short edge success, cost: %dms, edges size: %d\n", int(it), int(after - before), (int)mesh.n_edges());

		before = clock();
		equalize_valence();
		after = clock();
		printf("iter %d: equalize valence success, cost: %dms, edges size: %d\n", int(it), int(after - before), (int)mesh.n_edges());


		before = clock();
		tangential_relaxation(3);
		after = clock();
		printf("iter %d: tangential relaxation success, cost: %dms\n", int(it), int(after - before));

	}

}

void iGame_OpenMesh_Remeshing::split_long_edge() {

	int long_edge_count = 0;
	auto edges = mesh.edges();
	OM_VH vh, vh1, vhleft, vhright;
	OM_Vex v_new;
	double len;
	OM_HEH heh0, heh1;
	OM_FH fh0, fh1;
	
	
	for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it) {
		if (!mesh.is_valid_handle(*e_it)) continue;
		auto& eh = *e_it;
		if (mesh.is_boundary(eh)) continue;
		len = mesh.calc_edge_length(eh);
		if (len > upper_length) {
			v_new = mesh.calc_edge_midpoint(eh);
			vh = mesh.add_vertex(v_new);
			heh0 = e_it->h0();
			heh1 = e_it->h1();
			// 硬核边分裂，速度较慢……
			std::vector<std::vector<OM_VH>> old_vhs;
			if (!mesh.is_boundary(heh0)) {
				fh0 = mesh.face_handle(heh0);
				old_vhs.push_back(get_vhs(fh0));
			}
			if (!mesh.is_boundary(heh1)) {
				fh1 = mesh.face_handle(heh1);
				old_vhs.push_back(get_vhs(fh1));
			}
			mesh.delete_edge(eh);// 必须在添加面之前删除旧的边
			for (auto& _vhs : old_vhs) {
				for (int i = 0; i < 3; ++i) {
					if (_vhs[i] == e_it->v0() || _vhs[i] == e_it->v1()) {
						auto vhs_ = _vhs;
						vhs_[i] = vh;
						mesh.add_face(vhs_);
					}
				}
			}
			long_edge_count++;
		}
	}

	mesh.garbage_collection();

	std::cout << "Split " << long_edge_count << " long edges.\n";

}

void iGame_OpenMesh_Remeshing::collapse_short_edge() {

	int short_edge_count = 0;
	auto edges = mesh.edges();
	OM_VH vh;
	double len;
	OM_HEH heh;
	

	for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it) {
		auto& eh = *e_it;
		if (!mesh.is_valid_handle(eh) || mesh.is_boundary(eh)) continue;
		len = mesh.calc_edge_length(eh);
		if (len < lower_length) {
			auto vh1 = e_it->v0();
			auto vh2 = e_it->v1();
			// 边收缩，将 from 收缩到 to，from 度数小可加快收缩速度
			heh = (mesh.valence(vh1) < mesh.valence(vh2)) ? e_it->h0() : e_it->h1();
			if (mesh.is_collapse_ok(heh)) mesh.collapse(heh);
			short_edge_count++;
		}
	}
	int ios_vcnt = mesh.delete_isolated_vertices();
	mesh.garbage_collection();

	std::cout << "Collapse " << short_edge_count << " short edges. Delete " << ios_vcnt << " isolated vertices.\n";

}

void iGame_OpenMesh_Remeshing::equalize_valence() {

	int flip_cnt = 0;
	auto edges = mesh.edges();
	OM_VH vh0, vh1, vh2, vh3;
	OM_HEH heh0, heh1;
	int diff0, diff1, diff2, diff3, diff_ori, diff_flip;
	for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it) {
		auto& eh = *e_it;
		if (mesh.is_valid_handle(eh) && mesh.is_flip_ok(eh)) {
			vh0 = e_it->v0();
			vh1 = e_it->v1();
			heh0 = mesh.next_halfedge_handle(e_it->h0());
			heh1 = mesh.next_halfedge_handle(e_it->h1());
			vh2 = mesh.to_vertex_handle(heh0);
			vh3 = mesh.to_vertex_handle(heh1);
			diff0 = mesh.valence(vh0) - 6;
			diff1 = mesh.valence(vh1) - 6;
			diff2 = mesh.valence(vh2) - 6;
			diff3 = mesh.valence(vh3) - 6;
			diff_ori = diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
			diff0--, diff1--;
			diff2++, diff3++;
			diff_flip = diff0 * diff0 + diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
			if (diff_flip < diff_ori) {
				mesh.flip(eh);
				flip_cnt++;
			}
		}
	}
	mesh.garbage_collection();
	std::cout << "Flip " << flip_cnt << " edges.\n";
}

void iGame_OpenMesh_Remeshing::tangential_relaxation(size_t iter_times) {

	

	OpenMesh::FPropHandleT<iGame_OpenMesh_TriMesh::Point> centroids;
	mesh.add_property(centroids);
	OpenMesh::VPropHandleT<iGame_OpenMesh_TriMesh::Point> avgs;
	mesh.add_property(avgs);

	for (size_t it = 0; it < iter_times; ++it) {

		mesh.update_face_normals();
		mesh.update_vertex_normals();

		for (auto f_it = mesh.faces_sbegin(); f_it != mesh.faces_end(); ++f_it) {
			// 避免重复计算
			mesh.property(centroids, *f_it) = mesh.calc_face_centroid(*f_it);
		}


		int cnt;
		for (auto v_it = mesh.vertices_sbegin(); v_it != mesh.vertices_end(); ++v_it) {
			auto& vh = v_it;
			if (mesh.is_boundary(vh)) continue;
			mesh.property(avgs, vh).vectorize(0.f);
			cnt = 0;
			for (auto vf_it = mesh.vf_begin(vh); vf_it.is_valid(); ++vf_it) {
				mesh.property(avgs, vh) += mesh.property(centroids, *vf_it);
				cnt++;
			}
			if (cnt) {
				mesh.property(avgs, vh) /= cnt;
			}
		}

#pragma omp parallel for
		for (int i = 0; i < mesh.n_vertices(); ++i) {
			auto vh = mesh.vertex_handle(i);
			Eigen::Vector3d nor(mesh.normal(vh).data());
			Eigen::Vector3d pos(mesh.point(vh).data());
			Eigen::Vector3d avg(mesh.property(avgs, vh).data());
			Eigen::Vector3d u = avg + nor * (nor.transpose() * (pos - avg));
			Vector3f v_pos(u[0], u[1], u[2]), nv_pos;
			ab_tree->findNearstPoint(v_pos, nv_pos);
			mesh.set_point(vh, OM_Vex(nv_pos[0], nv_pos[1], nv_pos[2]));
		}

	}

	mesh.remove_property(centroids);
	mesh.remove_property(avgs);

}

std::vector<OM_VH> iGame_OpenMesh_Remeshing::get_vhs(OM_FH fh) {
	if (!mesh.is_valid_handle(fh)) {
		return {};
	}
	std::vector<OM_VH> vhs;
	for (auto fv_it = mesh.fv_begin(fh); fv_it.is_valid(); ++fv_it) {
		vhs.push_back(*fv_it);
	}
	return vhs;
}
