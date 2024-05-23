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

	// ����������ʲô��������, ���ת����һ������������һ�� AABB ���ļ���
	void get_Hausdorff_distance(MeshKernel::SurfaceMesh& src_mesh, AABB_Tree& aabb_tree, std::unordered_map<VH, Vector3d>& vh2color);

};