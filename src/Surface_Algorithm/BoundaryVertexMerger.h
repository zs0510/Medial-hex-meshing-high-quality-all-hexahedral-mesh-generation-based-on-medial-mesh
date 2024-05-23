#pragma once
#include "Kernel/Mesh.h"

class BoundaryVertexMerger {
public:

	BoundaryVertexMerger(MeshKernel::SurfaceMesh& _mesh): mesh(_mesh) {}

	void merge();

private:
	MeshKernel::SurfaceMesh& mesh;

};