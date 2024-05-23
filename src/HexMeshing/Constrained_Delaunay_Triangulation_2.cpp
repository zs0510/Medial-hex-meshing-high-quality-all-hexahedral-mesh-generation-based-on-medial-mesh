#include "Constrained_Delaunay_Triangulation_2.h"

void Constrained_Delaunay_Triangulation_2::execute(std::vector<std::vector<std::pair<double, double>>>& polygons,
    std::vector<std::pair<double, double>>& interior_points_2, MeshKernel::SurfaceMesh& mesh) {

    //construct some non-intersecting nested polygons
    int poly_cnt = polygons.size();
    std::vector<Polygon_2> polygons_2;
    for (auto& poly : polygons) {
        Polygon_2 tmp_poly;
        for (auto& point : poly) {
            tmp_poly.push_back(Point(point.first, point.second));
        }
        polygons_2.push_back(tmp_poly);
    }

    CDT cdt;

    for (auto& pos2 : interior_points_2) {// 可以在 bounding-box 内采样获得
        cdt.insert(Point(pos2.first, pos2.second));
    }

    //Insert the polygons into a constrained triangulation
    for (auto& poly : polygons_2) {
        cdt.insert_constraint(poly.vertices_begin(), poly.vertices_end(), true);
    }

    //Mark facets that are inside the domain bounded by the polygon
    mark_domains(cdt);

    int count = 0;
    for (Face_handle f : cdt.finite_face_handles()) {
        if (f->info().in_domain()) ++count;
    }

    MeshKernel::SurfaceMesh newmesh;

    std::unordered_map<CDT::Vertex_handle, VH> cgal_to_igame;
    for (CDT::Finite_vertices_iterator v_it = cdt.finite_vertices_begin(); v_it != cdt.finite_vertices_end(); v_it++) {
        Point p = v_it->point();
        cgal_to_igame[v_it] = newmesh.AddVertex(Vex(p.x(), p.y(), 0));
    }

    for (CDT::Finite_faces_iterator f_it = cdt.finite_faces_begin(); f_it != cdt.finite_faces_end(); f_it++) {
        if (!f_it->info().in_domain()) continue;
        std::vector<VH> newface;
        for (int i = 0; i < 3; ++i) {
            newface.push_back(cgal_to_igame[f_it->vertex(i)]);
        }
        newmesh.AddFace(newface);
    }

    mesh = newmesh;
    mesh.updateAllHandles();

    std::cout << "There are " << count << " facets in the domain." << std::endl;

}

void Constrained_Delaunay_Triangulation_2::execute(std::vector<std::vector<std::pair<double, double>>>& polygons,
    MeshKernel::SurfaceMesh& mesh) {

    //construct two non-intersecting nested polygons
    int poly_cnt = polygons.size();
    std::vector<Polygon_2> polygons_2;
    for (auto& poly : polygons) {
        Polygon_2 tmp_poly;
        for (auto& point : poly) {
            tmp_poly.push_back(Point(point.first, point.second));
        }
        polygons_2.push_back(tmp_poly);
    }
    
    CDT cdt;

    int cnt = 0;
    for (auto& vp : mesh.allvertices()) {
        if (!mesh.isOnBoundary(vp.first)) {
            if (cnt % 2) cdt.insert(Point(vp.second.x(), vp.second.y()));
            cnt++;
        }
    }

    //Insert the polygons into a constrained triangulation
    for (auto& poly : polygons_2) {
        cdt.insert_constraint(poly.vertices_begin(), poly.vertices_end(), true);
    }
    
    //Mark facets that are inside the domain bounded by the polygon
    mark_domains(cdt);

    int count = 0;
    for (Face_handle f : cdt.finite_face_handles()) {
        if (f->info().in_domain()) ++count;
    }

    MeshKernel::SurfaceMesh newmesh;

    std::unordered_map<CDT::Vertex_handle, VH> cgal_to_igame;
    for (CDT::Finite_vertices_iterator v_it = cdt.finite_vertices_begin(); v_it != cdt.finite_vertices_end(); v_it++) {
        Point p = v_it->point();
        cgal_to_igame[v_it] = newmesh.AddVertex(Vex(p.x(), p.y(), 0));
    }

    for (CDT::Finite_faces_iterator f_it = cdt.finite_faces_begin(); f_it != cdt.finite_faces_end(); f_it++) {
        if (!f_it->info().in_domain()) continue;
        std::vector<VH> newface;
        for (int i = 0; i < 3; ++i) {
           newface.push_back(cgal_to_igame[f_it->vertex(i)]);
        }
        newmesh.AddFace(newface);
    }

    mesh = newmesh;
    mesh.updateAllHandles();

    std::cout << "There are " << count << " facets in the domain." << std::endl;

}

void Constrained_Delaunay_Triangulation_2::mark_domains(CDT& ct,
    Face_handle start,
    int index,
    std::list<CDT::Edge>& border) {
    if (start->info().nesting_level != -1) {
        return;
    }
    std::list<Face_handle> queue;
    queue.push_back(start);
    while (!queue.empty()) {
        Face_handle fh = queue.front();
        queue.pop_front();
        if (fh->info().nesting_level == -1) {
            fh->info().nesting_level = index;
            for (int i = 0; i < 3; i++) {
                CDT::Edge e(fh, i);
                Face_handle n = fh->neighbor(i);
                if (n->info().nesting_level == -1) {
                    if (ct.is_constrained(e)) border.push_back(e);
                    else queue.push_back(n);
                }
            }
        }
    }
}
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void Constrained_Delaunay_Triangulation_2::mark_domains(CDT& cdt) {
    for (CDT::Face_handle f : cdt.all_face_handles()) {
        f->info().nesting_level = -1;
    }
    std::list<CDT::Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border);
    while (!border.empty()) {
        CDT::Edge e = border.front();
        border.pop_front();
        Face_handle n = e.first->neighbor(e.second);
        if (n->info().nesting_level == -1) {
            mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
        }
    }
}