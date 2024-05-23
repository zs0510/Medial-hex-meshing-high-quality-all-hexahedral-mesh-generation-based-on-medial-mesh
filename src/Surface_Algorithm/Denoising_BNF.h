#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <omp.h>

// Author: Zhang Sheng
// Time: 2021.8.26
// 
//MeshKernel::SurfaceMesh mesh = io.ReadFile("bunny_random.obj");
//Denoising_BNF test(mesh, 0.5, 5, 15);
//test.Execute();
//io.WriteFile(mesh, "bunny_denoise.obj");


class Denoising_BNF {
public:
	Denoising_BNF(MeshKernel::SurfaceMesh& _mesh, double _SigmaNormal = 0.5f, int it_normal_filtering = 5, int it_vertex_updating = 15) :
		mesh(_mesh), SigmaNormal(_SigmaNormal), it_normal(it_normal_filtering), it_vertex(it_vertex_updating) {
	};
	void Execute();

private:

	MeshKernel::SurfaceMesh& mesh;

	// 算法参数
	double SigmaNormal, SigmaCenter = 0.f;
	int it_normal, it_vertex;

	// 三角形数据
	std::vector<Eigen::Vector3d> centers;
	std::vector<Eigen::Vector3d> normals;
	std::vector<double> areas;

	// 顶点数据
	std::vector<Eigen::Vector3d> positions;

	double weight = 0.001;// user-specified first-trem positive weight, set to 0.001 by default
	
	void initTrianglesData();
	void initSigmaCenter();
	void initPositions();
	void BilateralNormalFiltering();
	void VertexUpdate();

	void VertexUpdate_PreventFlip();
	

};




