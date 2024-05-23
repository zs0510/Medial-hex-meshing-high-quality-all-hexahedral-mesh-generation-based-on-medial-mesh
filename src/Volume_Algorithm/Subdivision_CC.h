#pragma once
#include "Kernel/Mesh.h"
#include <omp.h>

#ifndef IGMAE_KERNEL_SIMPLIFY
#define IGMAE_KERNEL_SIMPLIFY
typedef MeshKernel::iGameVertex Vex;
typedef MeshKernel::iGameVertex Vec;
typedef MeshKernel::iGameVertexHandle VH;
typedef MeshKernel::iGameEdgeHandle EH;
typedef MeshKernel::iGameFaceHandle FH;
typedef MeshKernel::iGameCellHandle CH;
#endif

class Subdivision_CC {
public:
	Subdivision_CC(MeshKernel::VolumeMesh& _hexMesh) : mesh(_hexMesh){

	}
	void subdivision_ordinary(uint16_t iter_times);

private:
	MeshKernel::VolumeMesh& mesh;
	void add_new_cell();
	
};
