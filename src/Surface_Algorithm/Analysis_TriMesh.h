#pragma once
#include <Kernel/Mesh.h>
#include <omp.h>

#ifndef D_PI
#define D_PI 6.2831852f
#endif // !D_PI

typedef MeshKernel::iGameVertex Vex;
typedef MeshKernel::iGameVertexHandle VH;
typedef MeshKernel::iGameEdgeHandle EH;
typedef MeshKernel::iGameFaceHandle FH;
typedef MeshKernel::iGameCellHandle CH;

class Analysis_TriMesh {
private:
	MeshKernel::SurfaceMesh& mesh;
public:
	Analysis_TriMesh(MeshKernel::SurfaceMesh& _mesh) : mesh(_mesh) {

	}
	/******************** 计算曲率 ********************/
	void calculateGaussCurvature(std::vector<double>& curvatures);
	Vex getCircleCenter(FH);
	/******************** 计算角度 ********************/
	void calculateTriangleQuality(std::vector<double>& face_quality);
};
