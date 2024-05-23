#pragma once
#include <Kernel/Mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <iostream>

struct FaceInfo2 {
    FaceInfo2() {}
    int nesting_level;
    bool in_domain() {
        return nesting_level % 2 == 1;
    }
};


class Constrained_Delaunay_Triangulation_2 {

public:

    Constrained_Delaunay_Triangulation_2(){}

    void execute(std::vector<std::vector<std::pair<double, double>>>& polygons,
        MeshKernel::SurfaceMesh& mesh
        );// Note that the following code does not work if the boundaries of the polygons intersect.

    void execute(std::vector<std::vector<std::pair<double, double>>>& constrained_polygons, 
        std::vector<std::pair<double, double>>& interior_points_2, MeshKernel::SurfaceMesh& mesh);// 参数 1 是输入多边形，参数 2 是用来执行三角化的内部顶点，参数 3 是返回的三角化结果

private:

    typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
    typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
    typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>    Fbb;
    typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>        Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb, Fb>               TDS;
    typedef CGAL::Exact_predicates_tag                                Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
    typedef CDT::Point                                                Point;
    typedef CGAL::Polygon_2<K>                                        Polygon_2;
    typedef CDT::Face_handle                                          Face_handle;


    void mark_domains(CDT& ct, Face_handle start, int index, std::list<CDT::Edge>& border);
    void mark_domains(CDT& ct);


};
