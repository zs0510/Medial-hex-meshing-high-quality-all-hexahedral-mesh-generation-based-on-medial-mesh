#pragma once
#include "Kernel/Mesh.h"
#include <algorithm>
#include <vector>
#include <numeric>

class Subdivision_Butterfly {
public:
	Subdivision_Butterfly(MeshKernel::SurfaceMesh& _mesh) : mesh(_mesh) { }

	void subdivision(int iter_times);

private:

	MeshKernel::SurfaceMesh& mesh;
	std::map<int, std::vector<double>> weight_tables;
	std::unordered_map<EH, VH> edge_vh;

	const double nine_sixteen = 9.0 / 16.0;
	const double one_sixteen = 1.0 / 16.0;
	const double five_twelve = 5.0 / 12.0;
	const double one_twelve = 1.0 / 12.0;
	const double three_eight = 3.0 / 8.0;
	const double one_eight = 1.0 / 8.0;

	void init_weight_tables();// 初始化细分模板中各点的权值

	void generate_edge_vertex();

	void generate_new_faces();

	void get_ordered_vvhs(VH center, VH start, std::vector<VH>& res);
	void get_vex_from_extraordinary_vertices(VH center, VH strat, Vex& res);
};