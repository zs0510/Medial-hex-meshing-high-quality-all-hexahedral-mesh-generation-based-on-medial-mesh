#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include "Skeletal Mesh/SkeletalMesh.h"
#include "Kernel/Mesh.h"

using namespace std;
using namespace MeshKernel;

class HexMeshing {

public:
	HexMeshing(SkeletalMesh* in_skel_mesh) {
		sklMesh = in_skel_mesh;
	}

	void hexmeshing_based_on_skeleton();
	void find_partition_boundary(std::vector<std::vector<EH>>& ehs, std::vector<FH>& fhs);

private:
	SkeletalMesh* sklMesh;
	VolumeMesh* hexMesh = nullptr;

	void classify_skl_vertices();
	void resolution_control();

	bool find_partition_boundary_recursive(EH eh, std::vector<EH>& loop, std::unordered_map<EH, bool>& stored, std::unordered_map<EH, bool>& is_on_boundary);
	

};