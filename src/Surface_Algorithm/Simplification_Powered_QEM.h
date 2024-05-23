#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <queue>
#include <time.h>

// Author: Zhang Sheng
// Time: 2021.11.20

// compare with QEM
// this method can generate good triangles and preserve boundary
// but need more running time and weaker feature sensitive

#ifndef DOUBLE_LIMITS
#define DOUBLE_LIMITS
#define Double_MAX std::numeric_limits<double>::max()
#define Double_MIN -std::numeric_limits<double>::max()
#endif // !DOUBLE_LIMITS

#ifndef ASPECT_RATIO_THRESHOLD
#define ASPECT_RATIO_THRESHOLD 0.25f
#endif // !ASPECT_RATIO_THRESHOLD


using namespace std;

struct powered_contract {
	double cost, x, y, z;// �������ۼ�����λ��
	MeshKernel::iGameVertexHandle vh1, vh2;// ֻҪ��һ��Ϊ -1 �ͱ�ʾ��������Ϣ��Ч
	powered_contract(double c, double i, double j, double k, int v1_, int v2_) : cost(c), x(i), y(j), z(k), vh1(v1_), vh2(v2_) {};
	powered_contract() : cost(Double_MAX), x(0), y(0), z(0), vh1(-1), vh2(-1) {};
};

struct powered_cmp {
	bool operator()(powered_contract* a, powered_contract* b) {
		return a->cost > b->cost;
	}
};

class Simplification_Powered_QEM {
public:
	Simplification_Powered_QEM(MeshKernel::SurfaceMesh& _mesh, double _remainder = 0.5f) :// remainder is the num of remainder vertices
		mesh(_mesh), it_num(mesh.VertexSize() * (1 - _remainder)) {};
	Simplification_Powered_QEM(MeshKernel::SurfaceMesh& _mesh, int _remainder_vcnt) :// remainder_vcnt is the num of remainder vertices
		mesh(_mesh), it_num(mesh.VertexSize() - _remainder_vcnt) {
	};
	void Execute();
private:
	void initFaceNormal();// ����ÿ��ƽ��ķ�����
	void initVertexErrorMat();// ����ÿ���������������
	void initContractCosts();
	void initContractCost(MeshKernel::iGameEdge&);// ����ÿ���ߵ��������ۼ��������λ��

	double getAspectRatio(vector<Eigen::Vector3d>&);
	double getNormalizedRoundness(vector<Eigen::Vector3d>&);
	double getAspectRatio(MeshKernel::iGameFaceHandle&);

	Eigen::Vector3d getNormal(vector<Eigen::Vector3d>&);

	MeshKernel::SurfaceMesh& mesh;
	int it_num;// the num of removed vertices
	std::vector<Eigen::Matrix4d> ErrorMats;// ÿ���������������
	std::vector<Eigen::Vector3d> Normals;// ÿ����ķ�����
	std::priority_queue<powered_contract*, std::vector<powered_contract*>, powered_cmp> Contracts;
	std::vector<std::vector<powered_contract*> > adjContracts;
};

