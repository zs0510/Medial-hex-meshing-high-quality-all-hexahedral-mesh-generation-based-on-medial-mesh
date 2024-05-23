#pragma once
#include <Eigen/Core>
#include <map>
#include <omp.h>
#include "OpenMesh_Typedef.h"

// Author: Zhang Sheng
// Reference: Interpolating Subdivision for Meshes with Arbitrary Topology

class iGame_OpenMesh_Subdivision {
public:
	iGame_OpenMesh_Subdivision(OM_TriMesh& _mesh): mesh(_mesh) {}

	void subdivision_butterfly(int iter = 1);// Interpolating Subdivision for Meshes with Arbitrary Topology

private:
	OM_TriMesh& mesh;

	std::map<int, std::vector<double>> weight_tables;
	//OpenMesh::EPropHandleT<OM_VH> edge_vh;
	OpenMesh::HPropHandleT<OM_VH> edge_vh;

	void init_weight_tables();// 初始化细分模板中各点的权值

	void generate_edge_vertex();

	void generate_new_faces();

};
