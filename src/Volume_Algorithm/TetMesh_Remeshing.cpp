#include "TetMesh_Remeshing.h"

void TetMesh_Remeshing::tet_remeshing_via_offfile() {
//	// Domain
//	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//	typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
//	typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
//#ifdef CGAL_CONCURRENT_MESH_3
//	typedef CGAL::Parallel_tag Concurrency_tag;
//#else
//	typedef CGAL::Sequential_tag Concurrency_tag;
//#endif
//	// Triangulation for Meshing
//	typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
//	typedef CGAL::Mesh_complex_3_in_triangulation_3<
//		Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index> C3t3;
//	// Criteria
//	typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
//	// Triangulation for Remeshing
//	typedef CGAL::Triangulation_3<typename Tr::Geom_traits,
//		typename Tr::Triangulation_data_structure> Triangulation_3;
//	// To avoid verbose function and named parameters call
//	using namespace CGAL::parameters;
//	
//		const char* fname = (argc > 1) ? argv[1] : "data/fandisk.off";
//		std::ifstream input(fname);
//		Polyhedron polyhedron;
//		input >> polyhedron;
//		if (input.fail()) {
//			std::cerr << "Error: Cannot read file " << fname << std::endl;
//			return EXIT_FAILURE;
//		}
//		if (!CGAL::is_triangle_mesh(polyhedron)) {
//			std::cerr << "Input geometry is not triangulated." << std::endl;
//			return EXIT_FAILURE;
//		}
//		// Create domain
//		Mesh_domain domain(polyhedron);
//		// Get sharp features
//		domain.detect_features();
//		// Mesh criteria
//		Mesh_criteria criteria(edge_size = 0.025,
//			facet_angle = 25, facet_size = 0.05, facet_distance = 0.005,
//			cell_radius_edge_ratio = 3, cell_size = 0.05);
//		// Mesh generation
//		C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
//		Triangulation_3 tr = CGAL::convert_to_triangulation_3(std::move(c3t3));
//		//note we use the move semantic, with std::move(c3t3),
//		//  to avoid a copy of the triangulation by the function
//		//  `CGAL::convert_to_triangulation_3()`
//		//  After the call to this function, c3t3 is an empty and valid C3t3.
//		//It is possible to use :  CGAL::convert_to_triangulation_3(c3t3),
//		//  Then the triangulation is copied and duplicated, and c3t3 remains as is.
//		const double target_edge_length = 0.1;//coarsen the mesh
//		CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length,
//			CGAL::parameters::number_of_iterations(3));
//		

}

void TetMesh_Remeshing::tet_remeshing(double ratio) {

	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;
	typedef Remeshing_triangulation::Point CG_Point;
	typedef Remeshing_triangulation::Edge CG_Edge;
	typedef Remeshing_triangulation::Triangle CG_Tri;
	typedef Remeshing_triangulation::Cell CG_Cell;
	typedef Remeshing_triangulation::Vertex_handle CG_VH;

	double target_edge_length = 0;
	mesh.genAllEdgesLength();
	for (auto& ep : mesh.alledges()) {
		auto& egde = mesh.edges(ep.first);
		target_edge_length += egde.getLength();
	}
	target_edge_length /= mesh.esize();
	target_edge_length *= ratio;


	//// 1. ≥ı ºªØ CGAL TetMesh
	//std::unordered_map<VH, CG_VH> igame_to_cgal;

	///*for (auto& vp : mesh.allvertices()) {
	//	igame_to_cgal[vp.first] = (tr.insert(CG_Point(1, 1, 1)));
	//}

	//for (auto& cp : mesh.allcells()) {
	//	const auto& vhs = cp.second.getVertexHandle();
	//	tr.cell
	//}*/

	filename = "C:/My Files/Graphics/model_data/mesh_data/~dlut/MMC3D11/before_remeshing.mesh";
	std::ofstream outputfile(filename, std::ios::out);
	outputfile << "MeshVersionFormatted 1" << std::endl;
	outputfile << "Dimension 3" << std::endl;
	outputfile << "Vertices\n" << mesh.vsize() << std::endl;
	for (int i = 0; i < mesh.vsize(); ++i) {
		VH vh(i);
		auto& v = mesh.vertices(vh);
		outputfile << v.x() << " " << v.y() << " " << v.z() << " 1" << std::endl;
	}

	Remeshing_triangulation tr;





	std::ifstream is(filename, std::ios_base::in);
	CGAL::read_MEDIT(is, tr);

	std::cout << "CGAL: Initialize tetrahedral mesh success. #V = " << tr.number_of_vertices() << ", #C = " << tr.number_of_cells() << std::endl;

	CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

	std::ofstream os("C:/My Files/Graphics/model_data/mesh_data/~dlut/MMC3D11/after_remeshing.mesh");
	CGAL::write_MEDIT(os, tr);

}