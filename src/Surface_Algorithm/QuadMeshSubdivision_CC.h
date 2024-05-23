#pragma once
#include <unordered_map>
#include "Kernel/Mesh.h"

using namespace std;

class QuadMeshSubdivision_CC {
public:
	QuadMeshSubdivision_CC(MeshKernel::SurfaceMesh& _control_mesh): control_mesh(_control_mesh) {}

	void subdivision(int num_of_iteration);
	void get_subdivision_mesh(MeshKernel::SurfaceMesh& res) {
		res = subdi_mesh;
	}
	void get_control_weights(unordered_map<VH, unordered_map<VH, double>>& ctrl_weights) {
		ctrl_weights = weights;
	}

private:
	MeshKernel::SurfaceMesh& control_mesh;
	MeshKernel::SurfaceMesh subdi_mesh;
	
	// 与现有顶点相关的各个控制顶点的权值
	// weights[vh_sub][vh_ctrl] 表示控制顶点 vh_ctrl 对 vh_sub 的权值
	unordered_map<VH, unordered_map<VH, double>> weights;

};