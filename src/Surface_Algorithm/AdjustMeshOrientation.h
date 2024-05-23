#pragma once
#include <queue>
#include "Kernel/Mesh.h"

class MeshAdjuster {
public:
	MeshAdjuster() {}
	void execute(MeshKernel::SurfaceMesh& _mesh, FH right_fh);

private:


};

void MeshAdjuster::execute(MeshKernel::SurfaceMesh& mesh, FH right_fh) {// 请保证输入的是无 无效中间handle 的网格
	MeshKernel::SurfaceMesh newmesh;
	int vcnt = mesh.vsize(), fcnt = mesh.fsize();
	for (int i = 0; i < vcnt; ++i) {
		VH vh(i);
		if (!mesh.isValid(vh)) {
			std::cerr << "There exist invalid vh.\n";
			return;
		}
		auto& v = mesh.vertices(vh);
		VH newvh = newmesh.AddVertex(v);
		if (newvh != vh) {
			std::cerr << "Invalid newvh.\n";
			return;
		}
	}
	if (mesh.vsize() != newmesh.vsize()) {
		std::cerr << "Error: new mesh has wrong vertices number.\n";
		return;
	}

	std::unordered_map<FH, FH> old_to_new;
	std::queue<FH> que;// 旧网格的fh
	const auto& tmp_vhs = mesh.faces(right_fh).getVertexHandle();
	que.push(right_fh);
	old_to_new[right_fh] = newmesh.AddFace(tmp_vhs);
	
	while (!que.empty()) {
		FH oldfh = que.front();
		que.pop();
		FH newfh = old_to_new[oldfh];
		const auto& vhs_base = newmesh.faces(newfh).getVertexHandle();
		for (auto& adjfh : mesh.NeighborFh(oldfh)) {
			if (old_to_new.count(adjfh)) continue;
			que.push(adjfh);
			auto vhs_old = mesh.faces(adjfh).getVertexHandle();
			int prevh_old = vhs_old[0];
			int nextvh_old = vhs_old[1];
			if (prevh_old == vhs_base[0]) {
				if (nextvh_old == vhs_base[1]) {
					std::swap(vhs_old[0], vhs_old[2]);
					std::cout << "We reorder a face vhs.\n";
				}
			} else if (prevh_old == vhs_base[1]) {
				if (nextvh_old == vhs_base[2]) {
					std::swap(vhs_old[0], vhs_old[2]);
					std::cout << "We reorder a face vhs.\n";
				}
			} else if (prevh_old == vhs_base[2]) {
				if (nextvh_old == vhs_base[0]) {
					std::swap(vhs_old[0], vhs_old[2]);
					std::cout << "We reorder a face vhs.\n";
				}
			}
			old_to_new[adjfh] = newmesh.AddFace(vhs_old);

		}


	}
	mesh = newmesh;
	if (mesh.fsize() != newmesh.fsize()) {
		std::cerr << "Error: new mesh has wrong faces number.\n";
		std::cout << "Before: fcnt = " << mesh.fsize() << "; After: fcnt = " << newmesh.fsize() << "\n";
		return;
	}

	

}
