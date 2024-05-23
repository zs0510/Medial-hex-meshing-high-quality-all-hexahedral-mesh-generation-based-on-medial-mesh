#pragma once
#include <Kernel/Mesh.h>
#include <time.h>
#define random(a,b) (((double)rand()/RAND_MAX)*(b-a)+a)

typedef MeshKernel::iGameVertex Vex;
typedef MeshKernel::iGameVertex Vec;
typedef MeshKernel::iGameVertexHandle VH;
typedef MeshKernel::iGameEdgeHandle EH;
typedef MeshKernel::iGameFaceHandle FH;
typedef MeshKernel::iGameCellHandle CH;

class Add_Noise {
private:
	MeshKernel::VolumeMesh& mesh;
public:
	Add_Noise(MeshKernel::VolumeMesh& _mesh) : mesh(_mesh) {

	}
	void add_gauss_noise(double rate = 0.3f);
};

void Add_Noise::add_gauss_noise(double rate) {
	mesh.genAllEdgesLength();
	srand((unsigned int)time(0));
	for (auto& vp : mesh.allvertices()) {
		auto vh = vp.first;
		auto& v = mesh.vertices(vh);
		double e_x = random(0.f, 1.f);
		double e_y = random(0.f, 1.f);
		double e_z = random(0.f, 1.f);
		Vec dir(e_x, e_y, e_z);
		dir.normalize();
		//std::cout << "dir norm: " << dir.norm() << std::endl;
		double len_avg = 0.f;
		auto ehs = mesh.NeighborEh(vh);
		if (ehs.empty()) continue;
		for (auto& eh : ehs) {
			auto& edge = mesh.edges(eh);
			len_avg += edge.getLength();
		}
		len_avg /= ehs.size();
		v = v + dir * len_avg * rate;
	}
}
