#pragma once
#include "Kernel/Mesh.h"

// ������������������������֮���������Լ��

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
	std::unordered_map<int, std::unordered_map<int, int>> cc_lines_mapping;// �����ǵ�֮�����������
	std::vector<std::vector<int>> edges_chain;
	std::unordered_map<int, int> eid2ecid;// ��ID �� ��������ID ��ӳ��
	
	void feature_line_detect(double angle = 60);

	void build_graph();

};