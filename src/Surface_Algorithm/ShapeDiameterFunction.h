#pragma once
#include <vector>
#include <omp.h>
#include <time.h>
#include "Kernel/Mesh.h"
#include "Tools/BVH.h"

using namespace std;

class ShapeDiameterFunction {
public:
	ShapeDiameterFunction(MeshKernel::SurfaceMesh& _mesh) : mesh(_mesh) {

		bvh_tree.buildBVH_Tree(mesh);

	}

	void get_faces_sdf(std::vector<double>& FaceDiameter, bool smooth_flag = true, bool normalize_flag = false);

private:
	MeshKernel::SurfaceMesh& mesh;
	BVH_Tree bvh_tree;
	std::vector<double> faces_sdf;
	std::vector<double> normalized_faces_sdf;
	std::vector<Vex> faces_center;
	double sigma_center = 0, sigma_sdf = 0;// for smoothing which use gaussian function
	double sdf_min = std::numeric_limits<double>::max(), sdf_max = 0;

	void calc_faces_sdf();
	void calc_sigma();
	void smooth_sdf(int iter_times = 5);
	void sample_on_circle(const Vex& center, const Vex& normal, std::vector<Vex>& sample_points);
	void update_sdf_max_min();
	void normalize_sdf();

};