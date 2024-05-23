#pragma once
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Kernel/Mesh.h>
#include <Tools/HexQuality.h>

namespace VolumeEditer {

	enum GraftingEnergyMode {
		MIN_SJ, GRAFTING
	};
	
	double get_min_sj(MeshKernel::VolumeMesh& hexmesh, std::vector<Eigen::Vector3d>);
	double get_grafting_energy(std::vector<Eigen::Vector3d> );

	void add_cell(MeshKernel::VolumeMesh& hexmesh, FH fh1, FH fh2);
	void add_cell_connected(MeshKernel::VolumeMesh& hexmesh, FH fh1, FH fh2);
	void add_cell_disconnected(MeshKernel::VolumeMesh& hexmesh, FH fh1, FH fh2);

}