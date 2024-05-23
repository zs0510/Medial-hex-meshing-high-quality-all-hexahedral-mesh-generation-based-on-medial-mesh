#pragma once
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <fstream>
#include <string>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

#include "Tools/CircumCenter.h"
#include "Skeletal Mesh/SkeletalMesh.h"

using namespace std;

namespace VoronoiMesh {
	
	//typedef double simple_numbertype;
	//typedef CGAL::Cartesian<simple_numbertype> simple_kernel;
	//typedef CGAL::Polyhedron_3<simple_kernel, CGAL::Polyhedron_items_3> SimpleMesh;

	//typedef CGAL::Point_3< CGAL::Cartesian<simple_numbertype> > Point;

	class VertexInfo {
	public:
		int id;
		int tag;
		std::vector<int> tagvec;
		std::set<unsigned int> pole_bplist;
		double var_double;
	};

	class CellInfo {
	public:
		bool inside;
		int id;
		int tag;
		bool is_pole;
		std::set<unsigned int> pole_bplist;
		double dist_center_to_boundary; // approximate 
	};


	//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	//typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;
	//typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo, K> Cb;
	//typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
	//typedef CGAL::Delaunay_triangulation_3<K, Tds> Triangulation;


	///*class Triangulation_VoronoiMesh : public CGAL::Delaunay_triangulation_3<K, Tds> {
	//public:
	//	double TetCircumRadius(const Tetrahedron& tet);
	//	double TetLargestEdgeLength(const Tetrahedron& tet);
	//	double TetBoundingSphereRadius(const Tetrahedron& tet);
	//};*/

	//typedef Triangulation_VoronoiMesh::Vertex_handle Vertex_handle_t;
	//typedef Triangulation_VoronoiMesh::Cell_handle Cell_handle_t;

	//typedef Triangulation_VoronoiMesh::Vertex_iterator Vertex_iterator_t;
	//typedef Triangulation_VoronoiMesh::Edge_iterator Edge_iterator_t;
	//typedef Triangulation_VoronoiMesh::Facet_iterator Facet_iterator_t;
	//typedef Triangulation_VoronoiMesh::Cell_iterator Cell_iterator_t;

	//typedef Triangulation_VoronoiMesh::Finite_cells_iterator Finite_cells_iterator_t;
	//typedef Triangulation_VoronoiMesh::Finite_facets_iterator Finite_facets_iterator_t;
	//typedef Triangulation_VoronoiMesh::Finite_edges_iterator Finite_edges_iterator_t;
	//typedef Triangulation_VoronoiMesh::Finite_vertices_iterator Finite_vertices_iterator_t;

	//typedef Triangulation_VoronoiMesh::Facet_circulator Facet_circulator_t;
	//typedef Triangulation_VoronoiMesh::Cell_circulator Cell_circulator_t;

	//typedef Triangulation_VoronoiMesh::Point Point_t;


	////typedef CGAL::Polyhedron_3<K> Polyhedron;
	////typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;



	////typedef CGAL::Cartesian_d<double> CarK;
	////typedef CGAL::Min_sphere_annulus_d_traits_d<CarK> MSTraits;
	////typedef CGAL::Min_sphere_d<MSTraits> Min_sphere;

	////typedef CarK::Point_d MSPoint;

	void get_voronoi_mesh(vector<Eigen::Vector3d>& sample_points, SkeletalMesh& voronoi_mesh);

};





