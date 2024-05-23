#pragma once
#include <Eigen/Core>
#include <queue>
#include <vector>
#include <omp.h>
#include "../Kernel/Mesh.h"
#include "Tools/MeshMath.h"

using namespace std;
using namespace MeshKernel;

class HexMesh_Smoothing {
public:

	HexMesh_Smoothing(MeshKernel::VolumeMesh& _mesh);
	void Laplacian_Smoothing(int max_iter = 1, std::vector<int> feature_cells = {}, std::vector<int> feature_faces = {}, 
		std::vector<int> feature_edges = {}, std::vector<int> feature_vertices = {});
	
	void Laplacian_Level(int iter_max = 5, bool boundary_smooth_flag = true);

private:
	MeshKernel::VolumeMesh& mesh;
	double threshold_J = 0.9;
	vector<iGameVertexHandle> interior_vhs;
	vector<iGameVertexHandle> exterior_vhs;
	double CalculateCellQuality(MeshKernel::iGameCellHandle ch);
	Vex center;
	void update_bbox_center() {
		mesh.initBBox();
		center = (mesh.BBoxMax + mesh.BBoxMin) * 0.5;
	}
	void move_to_original_center() {
		mesh.initBBox();
		Vex current_center = (mesh.BBoxMax + mesh.BBoxMin) * 0.5;
		Vec move = center - current_center;
		for (auto& vp : mesh.allvertices()) {
			auto& v = mesh.vertices(vp.first);
			v += move;
		}
	}

	void Laplacian_Smoothing_Diffusion(int max_iter = 1, std::vector<int> feature_cells = {}, std::vector<int> feature_faces = {},
		std::vector<int> feature_edges = {}, std::vector<int> feature_vertices = {});
	void Laplacian_Smoothing_HC(std::vector<int> feature_vhs, int iter_max = 5);

};

class SurfaceMesh_Smoothing {
public:
	SurfaceMesh_Smoothing(MeshKernel::SurfaceMesh& _mesh): mesh(_mesh) {}

	void Laplacian_smoothing(int iter_tims);

private:
	MeshKernel::SurfaceMesh& mesh;

};
