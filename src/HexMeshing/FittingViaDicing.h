#pragma once
#include "Tools/AABB_Tree.h"
#include "Tools/KdTree.h"
#include "Surface_Algorithm/TriMesh_Generator.h"
#include "Volume_Algorithm/HexMesh_TopologicalOperators.h"
#include "VolumeFitting_CVIF.h"

using namespace std;

class FittingViaDicing {

public:
	struct HausdorffError {
		FH fh;// 六面体的面
		double dist;
		HausdorffError() { fh = FH(-1); dist = 0.0; }
		HausdorffError(FH _fh, double _dist) {
			fh = _fh;
			dist = _dist;
		}
	};
	struct HausdorffError_SmallHeapCMP {
		bool operator()(HausdorffError& he1, HausdorffError& he2) {
			return he1.dist < he2.dist;
		}
	};

public:
	FittingViaDicing(MeshKernel::VolumeMesh& hexmesh, MeshKernel::SurfaceMesh& refmesh, double Hausdorff_error_threshold = 0.03);
	
	void fitting();

private:
	MeshKernel::VolumeMesh& hexmesh;
	MeshKernel::SurfaceMesh& refmesh;
	double hausdorff_error_threshold = 0.03;
	double length_of_bbox;
	AABB_Tree* abtree_ref = nullptr;
	iGameKdTree* kdtree_hex = nullptr;
	vector<Vex> sample_points_ref;
	vector<pair<double, FH>> hausdorff_error_ref_to_hex;
	

	void init();

	void calc_hausdorff_dist_ref_to_hex();// 以六面体建立 KD 树
	
};