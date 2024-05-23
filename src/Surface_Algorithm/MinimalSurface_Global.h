#pragma once
#include "Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

// Author: Zhang Sheng
// Time: 2021.10.19

// Usage:
// MinimalSurface_Global tmp(triMesh);
// tmp.Execute();

#ifndef  EPSILON
#define EPSILON 1E-6F
#endif // ! EPSILON

typedef Eigen::Triplet<double> Tri;


class MinimalSurface_Global {
public:
	MinimalSurface_Global(MeshKernel::SurfaceMesh& _mesh): mesh(_mesh) {}
	void Execute();
private:
	MeshKernel::SurfaceMesh& mesh;
	void initCotWeight();
};