#pragma once
#include "Subdivision_Skl.h"
#include "BranchMerger.h"

class CurveResample {
public:
	CurveResample(SkeletalMesh& _mesh): mesh(_mesh) {}
	void resample();

private:
	SkeletalMesh& mesh;

	std::unordered_map<EH, bool> curve_ehs;

	std::unordered_map<VH, bool> branch_vhs;
	std::unordered_map<VH, bool> joint_vhs;
	std::unordered_map<VH, bool> end_vhs;

	void curve_subdivision(EH eh, int iter = 2);
	void recusive_resample(int left, int right, std::vector<VH>& curve_vhs, std::vector<bool>& used);
	void init_nodes_map();
};