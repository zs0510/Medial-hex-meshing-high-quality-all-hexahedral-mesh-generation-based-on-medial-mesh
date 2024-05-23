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
	double cost, x, y, z;// �������ۼ�����λ��
	MeshKernel::iGameVertexHandle vh1, vh2;// ֻҪ��һ��Ϊ -1 �ͱ�ʾ��������Ϣ��Ч
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
	void initFaceNormal();// ����ÿ��ƽ��ķ�����
	void initVertexErrorMat();// ����ÿ���������������
	void initContractCosts();
	void initContractCost(const EH&);// ����ÿ���ߵ��������ۼ��������λ��
	MeshKernel::SurfaceMesh& mesh;
	int it_num;// the num of removed vertices
	std::vector<Eigen::Matrix4d> ErrorMats;// ÿ���������������
	std::vector<Eigen::Vector3d> Normals;// ÿ����ķ�����
	std::priority_queue<contract*, std::vector<contract*>, cmp> Contracts;
	std::vector<std::vector<contract*> > adjContracts;
};

