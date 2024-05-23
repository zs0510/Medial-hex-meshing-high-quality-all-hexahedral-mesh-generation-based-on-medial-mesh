#pragma once
#include "AABB_Tree.h"
#include "Kernel/Mesh.h"
#include "ColorBar.h"
#include "Surface_Algorithm/TriMesh_Generator.h"

namespace ErrorRender {

	void get_Hausdorff_distance(MeshKernel::SurfaceMesh& src_mesh, MeshKernel::VolumeMesh& ref_mesh, std::unordered_map<VH, Vector3d>& vh2color);
	void get_Hausdorff_distance(MeshKernel::SurfaceMesh& src_mesh, MeshKernel::SurfaceMesh& ref_mesh, std::unordered_map<VH, Vector3d>& vh2color);
	

	void get_Hausdorff_distance(MeshKernel::VolumeMesh& src_mesh, MeshKernel::SurfaceMesh& ref_mesh, std::unordered_map<VH, Vector3d>& vh2color);
	void get_Hausdorff_distance(MeshKernel::VolumeMesh& src_mesh, MeshKernel::VolumeMesh& ref_mesh, std::unordered_map<VH, Vector3d>& vh2color);

	// 无论输入是什么样的网格, 最后都转化成一个表面网格与一个 AABB 树的计算
	void get_Hausdorff_distance(MeshKernel::SurfaceMesh& src_mesh, AABB_Tree& aabb_tree, std::unordered_map<VH, Vector3d>& vh2color);

};