#include "HoleFilling.h"

void SkeletalMesh_HoleFilling::hole_filling(SkeletalMesh& mesh, std::vector<VH> vhs) {

	if (vhs.size() < 3) return;
	for (auto& vh : vhs) {
		if (!mesh.isValid(vh)) return;
	}
	if (vhs.size() == 3) {
		mesh.AddFace(vhs);
		return;
	}
	/*std::unordered_set<VH> vhs_set;
	vhs_set.insert(vhs.begin(), vhs.end());

	vhs.clear();*/


}