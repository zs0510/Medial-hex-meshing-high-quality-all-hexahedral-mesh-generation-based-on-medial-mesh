#pragma once
#include <Eigen/Core>
#include "Tools/AABB_Tree.h"
#include <omp.h>
#include "OpenMesh_Typedef.h"

class iGame_OpenMesh_Remeshing {
private:
	iGame_OpenMesh_TriMesh& mesh;
public:
	iGame_OpenMesh_Remeshing(iGame_OpenMesh_TriMesh& _mesh) : mesh(_mesh) {
		init_aabb_tree();
		
	}
	void remeshing_uniform(double rate, size_t iter_times);// average smoothing

private:
	void init_aabb_tree();
	void init_target_length(double rate);
	AABB_Tree* ab_tree;
	double lower_length, upper_length;
	void split_long_edge();
	void collapse_short_edge();
	void equalize_valence();
	void tangential_relaxation(size_t);
	std::vector<OM_VH> get_vhs(OM_FH fh);
};