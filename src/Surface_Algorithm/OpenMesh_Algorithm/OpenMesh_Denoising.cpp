#include "OpenMesh_Denoising.h"

void iGame_OpenMesh_Denoising::initMeshData() {

	// 计算法向量
	mesh.request_face_normals();
	mesh.update_face_normals();
	mesh.add_property(normals);
	mesh.add_property(new_normals);
#pragma omp parallel for
	for (int i = 0; i < mesh.n_faces(); ++i) {
		auto fh = mesh.face_handle(i);
		mesh.property(normals, fh) = mesh.normal(fh);
	}

	// 计算质心
	mesh.add_property(centroids);
#pragma omp parallel for
	for (int i = 0; i < mesh.n_faces(); ++i) {
		auto fh = mesh.face_handle(i);
		mesh.property(centroids, fh) = mesh.calc_face_centroid(fh);
	}
	
	// 计算面积
	mesh.add_property(areas);
#pragma omp parallel for
	for (int i = 0; i < mesh.n_faces(); ++i) {
		auto fh = mesh.face_handle(i);
		mesh.property(areas, fh) = mesh.calc_face_area(fh);
	}

	
	sigma_centroid = 0;
	int cnt = 0;
	for (auto f_it = mesh.faces_sbegin(); f_it != mesh.faces_end(); ++f_it) {
		for (auto ff_it = mesh.ff_begin(*f_it); ff_it.is_valid(); ++ff_it) {
			sigma_centroid += (mesh.property(centroids, *f_it) - mesh.property(centroids, *ff_it)).norm();
			cnt++;
		}
	}
	sigma_centroid /= cnt;

	std::cout << "init mesh data success. sigma centroid = " << sigma_centroid << ", sigma normal = " << sigma_normal << std::endl;

}

void iGame_OpenMesh_Denoising::denoising_bnf() {

	
	BilateralNormalFiltering();

	VertexUpdating();
	
	// 记得回收存储空间
	mesh.remove_property(centroids);
	mesh.remove_property(areas);
	mesh.remove_property(normals);
	mesh.remove_property(new_normals);

}

void iGame_OpenMesh_Denoising::BilateralNormalFiltering() {

	double double_sigma_centoid_square = 2 * sigma_centroid * sigma_centroid;
	double double_sigma_normal_square = 2 * sigma_normal * sigma_normal;

	for (size_t it = 0; it < iter_normal; ++it) {

#pragma omp parallel for
		for (int i = 0; i < mesh.n_faces(); ++i) {
			auto fh = mesh.face_handle(i);
			std::set<OM_FH> adjfhs;
			// adjfhs are hanldes of faces who share common vertex with fh.
			for (auto fv_it = mesh.fv_begin(fh); fv_it.is_valid(); ++fv_it) {
				for (auto vf_it = mesh.vf_begin(*fv_it); vf_it.is_valid(); ++vf_it) {
					if (*vf_it == fh) continue;
					adjfhs.insert(*vf_it);
				}
			}
			double weight_sum = 0.f;
			mesh.property(new_normals, fh).vectorize(0.f);
			for (auto& adjfh : adjfhs) {
				double delta_centroid = (mesh.property(centroids, fh) - mesh.property(centroids, adjfh)).norm();
				double delta_normal = (mesh.property(normals, fh) - mesh.property(normals, adjfh)).norm();
				double weight_centroid = std::exp(-(delta_centroid * delta_centroid) / double_sigma_centoid_square);
				double weight_normal = std::exp(-(delta_normal * delta_normal) / double_sigma_normal_square);
				double weight_adjfh = weight_centroid * weight_normal * mesh.property(areas, adjfh);
				mesh.property(new_normals, fh) += weight_adjfh * mesh.property(normals, adjfh);
				weight_sum += weight_adjfh;
			}
			mesh.property(new_normals, fh) /= weight_sum;
			mesh.property(new_normals, fh).normalize_cond();
		}

#pragma omp parallel for
		for (int i = 0; i < mesh.n_faces(); ++i) {
			auto fh = mesh.face_handle(i);
			mesh.property(normals, fh) = mesh.property(new_normals, fh);
		}

	}

}

void iGame_OpenMesh_Denoising::VertexUpdating() {

	for (size_t i = 0; i < iter_vex; ++i) {

#pragma omp parallel for
		for (int i = 0; i < mesh.n_vertices(); ++i) {
			auto vh = mesh.vertex_handle(i);
			OM_Vec move(0, 0, 0);
			size_t count = 0;
			for (auto vf_it = mesh.vf_begin(vh); vf_it.is_valid(); ++vf_it) {
				auto nor = mesh.property(normals, *vf_it).data();
				Eigen::Vector3d N(nor[0], nor[1], nor[2]);
				auto dif = (mesh.property(centroids, *vf_it) - mesh.point(vh)).data();
				Eigen::Vector3d D(dif[0], dif[1], dif[2]);
				double weight = N.transpose() * D;
				move += weight * mesh.property(normals, *vf_it);
				count++;
			}
			move /= count;
			mesh.set_point(vh, mesh.point(vh) + move);
		}

	}

}