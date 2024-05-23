#include "HexMeshing_based_on_SkelMesh.h"

using namespace std;

void HexMeshing_based_on_SkelMesh::hexmeshing(SkeletalMesh* _skeletal_mesh, MeshKernel::HexMesh& hex_mesh) {

	sklMesh = _skeletal_mesh;
	if (!sklMesh) {
		std::cerr << "Skeletal mesh is null!\n";
		return;
	}

	hexmeshing_face();

	std::cout << "HexMeshing success.\n";

}

void HexMeshing_based_on_SkelMesh::hexmeshing_face() {

	detect_triangles();


}

void HexMeshing_based_on_SkelMesh::hexmeshing_curve() {

}

void HexMeshing_based_on_SkelMesh::detect_triangles() {

	

}