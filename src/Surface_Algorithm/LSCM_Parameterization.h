#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>


class LSCM_Parameterization {
public:
	LSCM_Parameterization(MeshKernel::SurfaceMesh& _mesh) : mesh(_mesh) {}
	bool Execute(int fixedPoint1, int fixedPoint2);
    bool Execute();
private:
	MeshKernel::SurfaceMesh& mesh;
	//Least	Square Conformal Maps
	int v_count, f_count, F_P, S_P;
	std::vector<int> VertexMapping, antiVertexMapping;
	Eigen::SparseMatrix<double> realMf, imageMf, realMp, imageMp, A2, BM, B_Sparse;
	Eigen::MatrixXd NewP;
	Eigen::Vector4d U;
	std::vector< std::vector<Eigen::Vector2d> > edgeVectors;
	std::vector<double> area, New_U, New_V;
	//每个3D空间三角形对应的2D平面三角形的三条边，按面编号存储
	void CalculateEdgeVectors();
	//固定2个顶点
	bool FixTwoBoundryPoint(int fixedPoint1, int fixedPoint2);
	//将顶点序号与坐标相对应
	void MapBackUV();
	//矩阵转化为稀疏矩阵
	void buildSparseMatrix(Eigen::SparseMatrix<double>& A1_sparse, Eigen::MatrixXd A, int A_rows, int A_cols);
	//构建A和BM矩阵
	void build_A2_BM();

};