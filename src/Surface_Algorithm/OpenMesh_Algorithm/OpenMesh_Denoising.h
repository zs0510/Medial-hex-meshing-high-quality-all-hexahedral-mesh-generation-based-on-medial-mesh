#pragma once
#include <Eigen/Core>
#include <omp.h>
#include "OpenMesh_Typedef.h"

// Author: Zhang Sheng
// Reference: Bilateral Normal Filtering for Mesh Denoising - zheng2011

class iGame_OpenMesh_Denoising {
public:
	iGame_OpenMesh_Denoising(OM_TriMesh& _mesh, double _sigma_normal = 0.5f, size_t _iter_normal = 5, size_t _iter_vex = 15) : 
		mesh(_mesh), sigma_normal(_sigma_normal), iter_normal(_iter_normal), iter_vex(_iter_vex) {
		initMeshData();
	}
	void denoising_bnf();
private:

	OM_TriMesh& mesh;
	OpenMesh::FPropHandleT<OM_Pos> centroids;
	OpenMesh::FPropHandleT<double> areas;
	OpenMesh::FPropHandleT<OM_Vec> normals;
	OpenMesh::FPropHandleT<OM_Vec> new_normals;
	double sigma_normal = 0.5f, sigma_centroid = 0.f;
	size_t iter_normal = 5, iter_vex = 15;
	void initMeshData();
	void BilateralNormalFiltering();
	void VertexUpdating();

};

