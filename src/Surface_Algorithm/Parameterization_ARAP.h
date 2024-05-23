#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <omp.h>
#include "Kernel/Mesh.h"
#include "LSCM_Parameterization.h"
#include "Tools/MeshMath.h"

using namespace std;

class Parameterization_ARAP {

public:
	Parameterization_ARAP(MeshKernel::SurfaceMesh& trimesh): mesh(trimesh) {
		vcnt = mesh.vsize();
		fcnt = mesh.fsize();
	}
	void execute(int iter_times = 15);

private:
	MeshKernel::SurfaceMesh& mesh;
	Eigen::MatrixX2d uv;
	vector<vector<double>> faces_cot;// vh1-vh2, vh2-vh0, vh0-vh1
	Eigen::MatrixXd local_coord;
	vector<Eigen::Triplet<double>> triplets;
	vector<Eigen::Matrix2d> Lts;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

	int vcnt, fcnt;

	void get_initial_parameterization();

	void calc_local_coord();

	void calc_all_cot();

	void build_coefficient_matrix();

	void calc_optimal_ration();

	void solve_linear_equations();

	//double getCot(FH fh, VH vh0, VH vh1);

};