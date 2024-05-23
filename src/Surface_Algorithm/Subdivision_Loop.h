#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>

class Subdivision_Loop {
public:
	Subdivision_Loop(MeshKernel::SurfaceMesh& _mesh, size_t _iter = 1): mesh(_mesh), iter(_iter){}
	void Execute();
private:
	MeshKernel::SurfaceMesh& mesh;
	size_t iter;
};