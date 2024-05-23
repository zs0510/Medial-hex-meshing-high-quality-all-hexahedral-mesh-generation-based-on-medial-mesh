#pragma once
#include <Eigen/Dense>
#include <Eigen/Core>
#include "Kernel/Mesh.h"
#ifndef DOUBLE_LIMITS
#define DOUBLE_LIMITS
#define DOUBLE_MAX std::numeric_limits<double>::max()
#define DOUBLE_MIN -std::numeric_limits<double>::max()
#endif // !DOUBLE_LIMITS

// 带有法向量的 KD Tree
using namespace std;
using namespace MeshKernel;

class KDTree6d {
	
public:
	struct NearestNeighbor {
		double dist;
		vector<Vector6d> face;
		Vector6d nearest;
		int tests;
		FH fh = FH(-1);
	};

	KDTree6d(MeshKernel::SurfaceMesh trimesh, unsigned int max_faces = 10, unsigned int max_depth = 30);
	NearestNeighbor nearest(Vector6d& v);

	double diagonal_length, weight_normal = 0.2;// default weight_normal = 0.2, 原论文六面体网格较密, 我们生成的六面体网格较为稀疏, 权值或应不同
private:
	
	typedef vector<vector<Vector6d>> Triangles;
	typedef vector<FH> TrianglesFH;
	typedef Vector6d Vec6;// (v, w * l * n): v is 3d-position, w is weight of normal, l is bounding box diagonal length, n is normal
	struct Node {
		Node() :faces(nullptr), left_child(nullptr), right_child(nullptr) {};
		unsigned int axis;
		double split;
		Triangles* faces;
		TrianglesFH* fhs;
		Node* left_child;
		Node* right_child;
		~Node() {
			delete faces;
			delete left_child;
			delete right_child;
		}
	};

	SurfaceMesh mesh;
	double dist_point_triangle(Vector6d& v, vector<Vector6d>& face, Vector6d& nearest_vertex);
	
	// Recursive part of build()
	unsigned int build_recurse(Node* node, unsigned int max_handles, unsigned int depth);

	// Recursive part of nearest()
	void nearest_recurse(Node* node, Vector6d& v, NearestNeighbor& data);

	Node* root = nullptr;
	

};