#pragma once
#include "Kernel/Mesh.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/IO/File_medit.h>


#include <iostream>
#include <fstream>
#include <string>




class TetMesh_Remeshing {
	

public:
	TetMesh_Remeshing(MeshKernel::VolumeMesh& _mesh, std::string _filename): mesh(_mesh), filename(_filename) {
		
	};
	void tet_remeshing(double target_length_ratio = 0.5f);
	void tet_remeshing_via_offfile();

private:
	MeshKernel::VolumeMesh& mesh;
	std::string filename;

};

