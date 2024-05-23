#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <set>


class bff {
public:
	bff(MeshKernel::SurfaceMesh& _mesh) :
		mesh(_mesh) {
	};
	void Init();
	double CotanAngles(MeshKernel::iGameEdgeHandle eh);
	void Parameterization();
	void InitCurvature();
	std::vector<MeshKernel::iGameEdgeHandle> AroundEdge(MeshKernel::iGameVertexHandle vh);
	MeshKernel::iGameEdgeHandle NextEdge(MeshKernel::iGameVertexHandle vh);
	void Laplace();
	void Execute();
private:
	MeshKernel::SurfaceMesh& mesh;
	int n, bn, in;
	std::vector<MeshKernel::iGameVertexHandle> points;
	std::map<int, int> to_new, to_old;
	std::map<int, bool> vis;
	Eigen::MatrixXd A, Aii, Aib, Abi, Abb;
	Eigen::MatrixXd K, k;
};