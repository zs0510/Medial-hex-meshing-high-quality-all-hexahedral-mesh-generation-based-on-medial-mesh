#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <set>
#include <unordered_map>
#include <omp.h>
#include <Kernel/Mesh.h>

using namespace std;

class SkeletalMesh : public MeshKernel::SurfaceMesh {
public:
	
	double scale_factor = 1.f;
	std::unordered_map<MeshKernel::iGameVertexHandle, double> radius;

	SkeletalMesh() {};
	bool readMaFile(std::string filename);// .ma
	bool writeMaFile(std::string filename);// .ma
	bool readSkelFile(std::string filename);// .skel
	bool writeSkelFile(std::string filename);// .skel

	void scale2Uint();
	void scale2Origin();
	void updateAllHandles();
	
};

