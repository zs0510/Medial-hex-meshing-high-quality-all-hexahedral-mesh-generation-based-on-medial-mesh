#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <queue>
#include <omp.h>
#include <time.h>

// Author: Zhang Sheng
// Time: 2021.8.27

#ifndef DOUBLE_LIMITS
#define DOUBLE_LIMITS
#define Double_MAX std::numeric_limits<double>::max()
#define Double_MIN -std::numeric_limits<double>::max()
#endif // !DOUBLE_LIMITS

struct contract {
	double cost, x, y, z;// 收缩代价及收缩位置
	MeshKernel::iGameVertexHandle vh1, vh2;// 只要有一个为 -1 就表示该收缩信息无效
	contract(double c, double i, double j, double k, int v1_, int v2_) : cost(c), x(i), y(j), z(k), vh1(v1_), vh2(v2_) {};
	contract() : cost(Double_MAX), x(0), y(0), z(0), vh1(-1), vh2(-1) {};
};

struct cmp {
	bool operator()(contract* a, contract* b) {
		return a->cost > b->cost;
	}
};

class Simplification_QEM {
public:
	Simplification_QEM(MeshKernel::SurfaceMesh& _mesh, double _remainder = 0.5f) :// remainder is the num of remainder vertices
		mesh(_mesh), it_num(mesh.VertexSize() * (1 - _remainder)) {};
	Simplification_QEM(MeshKernel::SurfaceMesh& _mesh, int _remainder_vcnt) :// remainder_vcnt is the num of remainder vertices
		mesh(_mesh), it_num(mesh.VertexSize() - _remainder_vcnt) {
	};
	void Execute();
private:
	void initFaceNormal();// 计算每个平面的法向量
	void initVertexErrorMat();// 计算每个顶点的误差二次型
	void initContractCosts();
	void initContractCost(const EH&);// 计算每条边的收缩代价及最佳收缩位置
	MeshKernel::SurfaceMesh& mesh;
	int it_num;// the num of removed vertices
	std::vector<Eigen::Matrix4d> ErrorMats;// 每个顶点的误差二次型
	std::vector<Eigen::Vector3d> Normals;// 每个面的法向量
	std::priority_queue<contract*, std::vector<contract*>, cmp> Contracts;
	std::vector<std::vector<contract*> > adjContracts;
};

