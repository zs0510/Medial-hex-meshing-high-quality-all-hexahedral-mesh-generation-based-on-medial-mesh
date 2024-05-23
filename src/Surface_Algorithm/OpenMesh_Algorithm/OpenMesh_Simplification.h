#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <queue>
#include <omp.h>
#include <time.h>
#include "OpenMesh_Typedef.h"

// Author: Zhang Sheng
// Reference: A New Mesh Simplification Algorithm based on QEM

struct CollapseNode {
	double cost;
	OM_Vex pos;
	OM_EH eh;
	bool valid;
	CollapseNode(double _cost, OM_Vex _pos, OM_EH _eh): cost(_cost), pos(_pos), eh(_eh), valid(true) {};
};

struct CollapseNode_CMP {
	bool operator() (CollapseNode* n1, CollapseNode* n2) {
		return n1->cost > n2->cost;
	}
};

class iGame_OpenMesh_Simplification {

public:
	iGame_OpenMesh_Simplification(iGame_OpenMesh_TriMesh& _mesh) : mesh(_mesh) {
		
	}
	void simplification_qem(size_t target_vertices_size);

private:
	iGame_OpenMesh_TriMesh& mesh;
	OpenMesh::VPropHandleT<Eigen::Matrix4d> error_mats;
	OpenMesh::VPropHandleT<std::vector<CollapseNode*>> collapse_nodes;
	std::priority_queue<CollapseNode*, std::vector<CollapseNode*>, CollapseNode_CMP> collapse_heap;
	void init_error_mats();
	void calc_collapse_cost(OM_EH);
	
};
