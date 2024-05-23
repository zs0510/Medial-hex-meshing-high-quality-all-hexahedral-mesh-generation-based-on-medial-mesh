#pragma once
#include "Kernel/Mesh.h"
#include <omp.h>

class Subdivision_Linear {
public:
	Subdivision_Linear(MeshKernel::VolumeMesh& _hexMesh) : mesh(_hexMesh) {

	}
	void subdivision(int iter_times, bool update_handll_flag = true);

private:
	MeshKernel::VolumeMesh& mesh;
	void add_new_cell();
	void add_new_cell_non_update_handle();

};