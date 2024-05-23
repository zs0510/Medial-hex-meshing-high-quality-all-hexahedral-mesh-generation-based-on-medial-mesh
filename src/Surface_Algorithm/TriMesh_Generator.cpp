#include "TriMesh_Generator.h"

void TriMesh_Generator::get_triangle_mesh(MeshKernel::VolumeMesh& input_mesh, MeshKernel::SurfaceMesh& trimesh) {
	
	
	// 1. 添加点
	for (auto& vp : input_mesh.allvertices()) {
		if (input_mesh.isOnBoundary(vp.first)) {
			oldvh_newvh[vp.first] = trimesh.AddVertex(vp.second);
		}
	}
	// 2. 添加面
	for (auto& fp : input_mesh.allfaces()) {
		if (!input_mesh.isOnBoundary(fp.first)) continue;
		auto ch = *(input_mesh.NeighborCh(fp.first).begin());
		auto cell_center = input_mesh.getCellCenter(ch);
		auto face_center = input_mesh.getFaceCenter(fp.first);
		Vec out_dir = (face_center - cell_center).normalized();
		auto vhs = fp.second.getVertexHandle();
		for (int i = 2; i < vhs.size(); ++i) {
			bool is_valid = true;
			std::vector<VH> tri_vhs = { vhs[0], vhs[i-1], vhs[i] };
			for (auto& vh : tri_vhs) {
				if (!oldvh_newvh.count(vh)) {
					std::cerr << "Error: Face is not on the boundary. Donot need to save to triangle mesh.\n";
					is_valid = false;
					break;
				}
				vh = oldvh_newvh[vh];
			}
			if (is_valid) {
				auto& v0 = trimesh.vertices(tri_vhs[0]);
				auto& v1 = trimesh.vertices(tri_vhs[1]);
				auto& v2 = trimesh.vertices(tri_vhs[2]);
				Vec vec01 = v1 - v0;
				Vec vec02 = v2 - v0;
				Vec tri_normal = (vec01.cross(vec02)).normalized();
				double cosine = tri_normal.dot(out_dir);
				if (cosine < 0) {
					std::swap(tri_vhs[0], tri_vhs[2]);
				}
				FH newfh = trimesh.AddFace(tri_vhs);
				newfh_oldfh[newfh] = fp.first;
			}
		}
	}

}

Vex TriMesh_Generator::get_center(const std::vector<Vex>& vertices) {
	Vex avg(0, 0, 0);
	int cnt = vertices.size();
	if (cnt == 0) return avg;
	for (auto& v : vertices) {
		avg += v;
	}
	avg /= cnt;
	return avg;
}