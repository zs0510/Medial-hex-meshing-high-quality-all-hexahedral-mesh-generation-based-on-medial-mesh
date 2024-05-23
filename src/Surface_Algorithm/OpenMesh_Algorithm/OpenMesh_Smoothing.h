#pragma once
#include <omp.h>
#include "OpenMesh_Typedef.h"

class iGame_OpenMesh_Smoother {
private:
	iGame_OpenMesh_TriMesh& mesh;
public:
	iGame_OpenMesh_Smoother(iGame_OpenMesh_TriMesh& _mesh): mesh(_mesh) {
	
	}
	void smoothing(size_t iter_times);// Laplacian Smoothing
};