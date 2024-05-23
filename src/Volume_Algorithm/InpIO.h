#pragma once
#include "Mesh.h"
#include <fstream>
#include <string>
#include <sstream>
#include <array>

namespace InpIO {

	MeshKernel::VolumeMesh get_volume_mesh_from_inp_file(std::string filename);

}