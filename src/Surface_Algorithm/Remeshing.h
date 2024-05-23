#pragma once
#include <Kernel/Mesh.h>
#include "Tools/AABB_Tree.h"
#include <Eigen/Core>
#include <omp.h>
#include <time.h>

#ifndef PI
#define PI 3.14159265358979
#endif // !PI

#ifndef DOUBLE_LIMITS
#define DOUBLE_LIMITS
#define Double_MAX std::numeric_limits<double>::max()
#define Double_MIN -std::numeric_limits<double>::max()
#endif // !DOUBLE_LIMITS

class Remeshing {
public:
	Remeshing(MeshKernel::SurfaceMesh& _mesh, bool _uniform_flag);
	void Execute();
private:
	MeshKernel::SurfaceMesh& mesh;
	const double rad2angle = 180 / PI;
	std::unordered_map<int, double> edge_norm2;
	double target_length, lower_length, upper_length;
	double lower_angle = 30, upper_angle = 90;
	double flip_ok_angle = 10;
	bool uniform_flag = false;
	bool flip_check_flag = true;
	bool convergence_flag = false;
	AABB_Tree* ab_tree;
	void initAABBTree();
	void initTargetEdgeLength();
	void initEdgesNorm2();
	void calEdgeNorm2(MeshKernel::iGameEdgeHandle eh);
	void removeLargeAngle(int count);
	void removeSmallAngle(int count);
	void equalizeValence();
	void tangentialRelaxation(int iter = 5);
	void projectToSurface();
	void splitLongEdge();
	void collapseShortEdge();
	bool is_flip_ok(MeshKernel::iGameEdgeHandle eh);
	void flip(MeshKernel::iGameEdgeHandle eh);
	double get_square_difference(MeshKernel::iGameEdgeHandle eh);

};

