#pragma once
#include "Kernel/Mesh.h"
#include "Kernel/IO.h"
#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/SurfaceRemeshing.h>
#include <omp.h>

class Remeshing_2 {

public:
	Remeshing_2(MeshKernel::SurfaceMesh& _mesh):mesh(_mesh) {

	}

	void execute_pmp(double edge_length_ratio = 0.5);
	void execute_igame(double edge_length_ratio = 0.5);

private:
	MeshKernel::SurfaceMesh& mesh;

	void split_long_edges(double target_len);
	void collapse_short_edges(double target_len);
	void equalize_valence();
	void relaxation(int iter = 5);

	void flip(EH);

};

void Remeshing_2::execute_pmp(double ratio) {

	std::string input_file = "C:/My Files/Graphics/model_data/obj_data/remeshing_2_in.obj";
	std::string output_file = "C:/My Files/Graphics/model_data/obj_data/remeshing_2_out.off";
	{
		MeshKernel::IO igame_io;
		igame_io.WriteObjFile(mesh, input_file);
	}
	
	pmp::SurfaceMesh pmpmesh;

	pmpmesh.read(input_file);

	printf("Remeshing_Boundary_Perserving: initial mesh:\t vertices: %d, edges: %d, faces: %d\n", pmpmesh.n_vertices(), pmpmesh.n_edges(), pmpmesh.n_faces());
	double avg_edge_length = 0;

	for (auto it = pmpmesh.edges_begin(); it != pmpmesh.edges_end(); ++it) {
		double len = pmpmesh.edge_length(*it);
		avg_edge_length += len;
	}
	avg_edge_length /= pmpmesh.edges_size();

	pmp::SurfaceRemeshing app(pmpmesh);
	//app.uniform_remeshing(avg_edge_length * 0.5);
	//app.adaptive_remeshing(0, avg_edge_length, avg_edge_length * 0.5);
	app.uniform_remeshing(avg_edge_length * ratio);

	pmpmesh.write(output_file);

	{
		MeshKernel::IO igame_io;
		mesh = igame_io.ReadOffFile(output_file);
	}

}