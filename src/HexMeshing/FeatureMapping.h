#pragma once
#include "Kernel/Mesh.h"

// 建立参照网格与六面体网格之间的特征线约束

class FeatureMapping {

public:
	FeatureMapping(MeshKernel::SurfaceMesh& _refmesh, 
		MeshKernel::VolumeMesh& _hexmesh, 
		std::vector<int> _corners_ref, 
		std::vector<int> _corners_hex);

	

	void build_mapping(std::unordered_map<VH, int>& res_line_features_hex,
		std::vector<std::vector<EH>>& res_line_features_tri);

private:
	MeshKernel::SurfaceMesh& refmesh;
	MeshKernel::VolumeMesh& hexmesh;
	std::vector<int> corners_ref;
	std::vector<int> corners_hex;
	std::unordered_map<int, int> corners_cons_ref2hex;
	std::unordered_map<int, int> corners_cons_hex2ref;

	std::unordered_set<int> edge_features_ref;
	std::unordered_set<int> corner_features_ref;
	std::unordered_map<int, std::unordered_map<int, int>> cc_lines_mapping;// 两个角点之间的特征边链
	std::vector<std::vector<int>> edges_chain;
	std::unordered_map<int, int> eid2ecid;// 边ID 到 所属边链ID 的映射
	
	void feature_line_detect(double angle = 60);

	void build_graph();

};