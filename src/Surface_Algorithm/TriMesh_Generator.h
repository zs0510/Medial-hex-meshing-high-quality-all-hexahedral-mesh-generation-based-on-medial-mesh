#pragma once
#include "Kernel/Mesh.h"

class TriMesh_Generator {
public:
	TriMesh_Generator(){};
	void get_triangle_mesh(MeshKernel::VolumeMesh& volume_mesh, MeshKernel::SurfaceMesh& trimesh);
	void get_oldvh_to_newvh(std::unordered_map<VH, VH>& res) { res = oldvh_newvh; }
	void get_newfh_to_oldfh(std::unordered_map<FH, FH>& res) { res = newfh_oldfh; }
private:
	std::unordered_map<VH, VH> oldvh_newvh;
	std::unordered_map<FH, FH> newfh_oldfh;
	Vex get_center(const std::vector<Vex>& vertices);

};