#include "OpenMesh_Smoothing.h"

void iGame_OpenMesh_Smoother::smoothing(size_t iter_times) {

	OpenMesh::VPropHandleT<iGame_OpenMesh_TriMesh::Point> centroids;
	mesh.add_property(centroids);
	
	for (size_t it = 0; it < iter_times; ++it) {

#pragma omp parallel for
		for (int i = 0; i < mesh.n_vertices(); ++i) {
			auto vh = mesh.vertex_handle(i);
			mesh.property(centroids, vh).vectorize(0.f);
			for (auto vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it) {
				mesh.property(centroids, vh) += mesh.point(*vv_it);
			}
			auto val = mesh.valence(vh);
			if (val > 0) {
				mesh.property(centroids, vh) /= val;
			}
		}

#pragma omp parallel for
		for (int i = 0; i < mesh.n_vertices(); ++i) {
			auto vh = mesh.vertex_handle(i);
			if (!mesh.is_boundary(vh)) {
				mesh.set_point(vh, mesh.property(centroids, vh));
			}
		}

	}

	mesh.remove_property(centroids);

}