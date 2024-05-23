#pragma once
#include "../Kernel/Mesh.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <unordered_map>
#include <omp.h>
//#pragma omp parallel for

using namespace std;
using namespace MeshKernel;

struct transform_node {
	iGameCellHandle ch;
	Eigen::Vector3d pos;
	transform_node(): ch(-1), pos(0, 0, 0){}
	transform_node(iGameCellHandle _ch, Eigen::Vector3d _pos): ch(_ch), pos(_pos){}
};

class HexMesh_Regularization {
public:

	HexMesh_Regularization(VolumeMesh& _mesh, bool _use_exterior_only ): mesh(_mesh), use_exterior_only(_use_exterior_only){

	}

	void regularization(uint32_t iter);

private:

	VolumeMesh& mesh;
	bool use_exterior_only = false;
	bool geometric_constraint_flag = false;
	bool check_weight_flag = false;
	bool check_coordinate_flag = false;
	double scaling_factor = 0.7f;// scaling factor, should > 0
	double weight_e = 5;// Local Regularization use, the weight of edge direction( >= 5 is ok)
	double weight_d = 5;// Global Optimization use, the weight of node displacement, default: 5
	double weight_g = 10;// Global Optimization use, the weight of geometric constraints, it is much larger than weight_d, default: 100
	double factor_e = 1.f / 6;
	double KEEPOPSITIVE = 1E-8F;
	double in_min_scaled_jacobian = 2.f, in_avg_scaled_jacobian = 0.f;
	double out_min_scaled_jacobian = 2.f, out_avg_scaled_jacobian = 0.f;

	unordered_map<int, vector<transform_node>> transforms;// save the temporary transformed element nodes
	unordered_map<int, double> weight_cell;// Global Optimization use, shoule be opsitive
	unordered_map<int, double> weight_edge;// Global Optimization use, shoule be opsitive
	unordered_map<int, Eigen::Vector3d> normal_vertex;// vh and its normal
	unordered_map<int, Eigen::Vector3d> ideal_position;// vh and its ideal position

	vector<vector<int>> tetFaces = {
		{1, 0, 2, 5}, {3, 2, 0, 7},
		{4, 0, 5, 7}, {6, 7, 5, 2}
	};

	void CalculateCellQuality(vector<iGameCellHandle>& chs, int mode);// mode 0 estemm initial mesh, mode 1 esteem result mesh
	void Preprocessing(vector<iGameEdgeHandle>& ehs, vector<iGameCellHandle>& chs);
	void LocalRegularization(const iGameCellHandle& ch);
	void GlobalOptimization();
	void ComputeTetrahedronLaplacianWeight(vector<iGameVertexHandle> vhs);//  the weights should be positive

	
	void GlobalOptimization_Exterior();
	void LaplacianSmoothing_Interior();

};
