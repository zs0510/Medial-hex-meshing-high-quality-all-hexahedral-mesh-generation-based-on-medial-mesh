#include "Project_to_Surface.h"

void Project_to_Surface::project_to_nearest_point() {

	// 1. 初始化 AABB 树
	std::vector<Vector3f> vertices;
	for (auto& fp : refmesh.allfaces()) {
		auto vhs = fp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& v = refmesh.vertices(vh);
			vertices.push_back(Vector3f(v.x(), v.y(), v.z()));
		}
	}
	ab_tree = new AABB_Tree(vertices);


	project_exterior();

	HexMesh_Smoothing app(hexmesh);
	app.Laplacian_Level(5, false);


}

void Project_to_Surface::project_exterior() {

	// 2. 投影六面体的每个表面顶点
	int vcnt = hexmesh.vsize();

#pragma omp parallel for
	for (int vid = 0; vid < vcnt; ++vid) {
		VH vh(vid);
		if (hexmesh.isOnBoundary(vh)) {
			auto& v = hexmesh.vertices(vh);
			Vector3f pos(v.x(), v.y(), v.z());
			Vector3f nearest_pos;
			ab_tree->findNearstPoint(pos, nearest_pos);
			v.setPosition(Vector3d(nearest_pos[0], nearest_pos[1], nearest_pos[2]));
		}
	}

}