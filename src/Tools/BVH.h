#pragma once
#include "Kernel/Mesh.h"
#include "TinyVector.h"
#include <queue>
#include "../Dependence/Geom/CPoint.h"
#include "Skeletal Mesh/SkeletalMesh.h"

// Author: Zhang Sheng
// Date: 2022.2.20

struct Ray;
struct Bound3;
struct Intersection;
struct BVH_Face;
struct BVH_Node;

struct Ray {
	Vector3d pos, dir, invDir;
	bool isNeg[3] = { false, false, false };
	Ray() {};
	Ray(Vector3d _pos, Vector3d _dir) : pos(_pos), dir(_dir) {
		dir.Normalize();
		if (dir[0] < 0.f) invDir[0] = 1.f / std::min(dir[0], (double)-1E-6F);
		else invDir[0] = 1.f / std::max(dir[0], (double)1E-6F);
		if (dir[1] < 0.f) invDir[1] = 1.f / std::min(dir[1], (double)-1E-6F);
		else invDir[1] = 1.f / std::max(dir[1], (double)1E-6F);
		if (dir[2] < 0.f) invDir[2] = 1.f / std::min(dir[2], (double)-1E-6F);
		else invDir[2] = 1.f / std::max(dir[2], (double)1E-6F);
		isNeg[0] = (invDir[0] < 0.f);
		isNeg[1] = (invDir[1] < 0.f);
		isNeg[2] = (invDir[2] < 0.f);
	};
};

struct Bound3 {// 三维包围盒
	Vector3d pos_min, pos_max, centroid;
	Bound3() {
		centroid = pos_min = pos_max = Vector3d(0.f, 0.f, 0.f);
	};
	Bound3(Vector3d _min, Vector3d _max) : pos_min(_min), pos_max(_max) {
		centroid = (pos_min + pos_max) / 2;
	};
};

struct Intersection {
	bool happened;// 是否相交
	double distance;
	Vector3d pos;
	int index, intersection_count;
	Intersection() {
		happened = false;
		index = -1;
		intersection_count = 0;
		distance = std::numeric_limits<double>::max();
		pos = Vector3d(0.f, 0.f, 0.f);
	};
};

struct BVH_Face {
	std::vector<Vector3d> vertices;
	int index;
	BVH_Face() {
		index = -1;
	};
};

struct BVH_Node {
	Bound3 box;
	BVH_Face face;// 此对象只在叶子结点中有效
	double area;
	BVH_Node* left;
	BVH_Node* right;
	BVH_Node() {
		area = 0.f;
		left = nullptr;
		right = nullptr;
	};
};



class BVH_Tree {
private:

	int leaves = 0, height = 0;
	BVH_Node* root;
	

	Bound3 Union(const Bound3&, const Bound3&);
	bool checkIntersect(const Ray&, const Bound3&);
	Intersection getIntersection(const Ray&, const std::vector<Vector3d>&);
	Intersection getIntersection(const Ray&, const BVH_Face&);
	Intersection getIntersection(const Ray&, BVH_Node*);

	BVH_Node* recursiveBuild(std::vector<BVH_Node*>& nodes, int _height);

public:
	//int intersection_count = 0;
	Intersection getIntersection(const Ray& ray) {// 检查相交，并返回交点。交点index = -1 表示未相交
		//intersection_count = 0;
		return getIntersection(ray, root);
	}

	void buildBVH_Tree(MeshKernel::VolumeMesh&);
	void buildBVH_Tree(MeshKernel::SurfaceMesh&);
	void buildBVH_Tree(MeshKernel::TetMesh&);
	void buildBVH_Tree(SkeletalMesh&);
	void buildBVH_Tree(std::vector<BVH_Face>&);

	void destoryBVH_Tree();
};

