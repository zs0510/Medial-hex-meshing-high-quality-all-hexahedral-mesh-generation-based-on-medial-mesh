#pragma once
#include <Eigen/Dense>
#include <Eigen/Core>
#include "Kernel/Mesh.h"

#ifndef DOUBLE_LIMITS
#define DOUBLE_LIMITS
#define DOUBLE_MAX std::numeric_limits<double>::max()
#define DOUBLE_MIN -std::numeric_limits<double>::max()
#endif // !DOUBLE_LIMITS

// Author: ZS
// Time: 2021.11.27
/*
	// 创建一个 KD树
	TriangleKdTree* kd_tree = new TriangleKdTree(*surface_mesh, 0);
	// 找到最近的面以及最近点
	auto face_nearest = kd_tree->nearest(v);
	Vex v_kd = face_nearest.nearest;
*/

namespace MeshKernel {

	class iGameKdTree {
	public:
		iGameKdTree(SurfaceMesh _mesh, unsigned int max_faces = 10, unsigned int max_depth = 30);
		
		struct NearestNeighbor {
			double dist;
			FH fh;
			MeshKernel::iGameVertex nearest;
			int tests;
		};
		NearestNeighbor nearest(MeshKernel::iGameVertex& v);

	private:

		bool is_medial_mesh = false;

		SurfaceMesh mesh;
		
		typedef std::vector<FH> Triangles;
		typedef Eigen::Vector3d Vec3;

		struct Node {
			Node() :faces(nullptr), left_child(nullptr), right_child(nullptr) {};
			unsigned int axis;
			double split;
			Triangles* faces;
			Node* left_child;
			Node* right_child;
			~Node() {
				delete faces;
				delete left_child;
				delete right_child;
			}
		};

		double dist_point_triangle(MeshKernel::iGameVertex& v, FH& face, MeshKernel::iGameVertex& nearest_vertex);
		double dist_point_line_segment(MeshKernel::iGameVertex& v, MeshKernel::iGameVertex& v0, MeshKernel::iGameVertex& v1, MeshKernel::iGameVertex& nearest_vertex);

		// Recursive part of build()
		unsigned int build_recurse(Node* node, unsigned int max_handles, unsigned int depth);

		// Recursive part of nearest()
		void nearest_recurse(Node* node, MeshKernel::iGameVertex& v, NearestNeighbor& data);

		Node* root = nullptr;

	};

};