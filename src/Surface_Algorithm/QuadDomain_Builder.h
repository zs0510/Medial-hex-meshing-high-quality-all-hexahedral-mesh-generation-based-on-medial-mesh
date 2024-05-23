#pragma once
#include <vector>
#include "Kernel/Mesh.h"

using namespace std;

class QuadDomain_Builder {
public:
	QuadDomain_Builder(MeshKernel::SurfaceMesh& quad_mesh): mesh(quad_mesh){
		
		for (auto& vp : mesh.allvertices()) {
			auto& vh = vp.first;
			if (mesh.isOnBoundary(vh)) boundary_vhs[vh] = true;
			else if (mesh.valence(vh) != 4) {
				singular_vhs[vh] = true;
			}
		}

	}

	void get_quad_domains(vector<int>& ehs);
	void generate_quad_mesh(MeshKernel::SurfaceMesh& quad_mesh, vector<int>& vhs);

private:
	MeshKernel::SurfaceMesh& mesh;
	unordered_map<VH, bool> singular_vhs;
	unordered_map<VH, bool> boundary_vhs;
	unordered_map<VH, unordered_map<VH, bool>> singular_visited;
	vector<vector<EH>> quad_edges;
	vector<VH> quad_edges_srcvh;

	unordered_map<VH, bool> bdy_singular_vhs;
	

	void march_inr(EH, VH src, VH tgt, vector<EH>&, unordered_map<EH, bool>&);
	void march_bdy(EH, VH src, VH tgt, vector<EH>&, unordered_map<EH, bool>&);
	EH get_next_eh(EH, VH src, VH tgt);

	

};