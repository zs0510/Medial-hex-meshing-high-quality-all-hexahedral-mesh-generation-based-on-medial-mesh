#pragma once
#include "Kernel/Mesh.h"

using namespace std;

namespace HexAddCell {

	void add_cell_2fh(MeshKernel::VolumeMesh& mesh, FH fh1, FH fh2);// 通过现存的两个面增加体

	void add_cell_3fh(MeshKernel::VolumeMesh& mesh, FH fh1, FH fh2, FH fh3);// 通过现存的三个面增加体, 这三个面各自均相连

	void add_cell_eh_fh(MeshKernel::VolumeMesh& mesh, EH eh, FH fh);// 通过现存的一边一面增加体

	void add_cell_fh(MeshKernel::VolumeMesh& mesh, FH fh);// 通过现存的一个面增加体

};
