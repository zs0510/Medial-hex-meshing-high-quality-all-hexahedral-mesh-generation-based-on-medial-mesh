#include "VoronoiMesh.h"

namespace VoronoiMesh {
	

	void get_voronoi_mesh(vector<Eigen::Vector3d>& sample_points, SkeletalMesh& voronoi_mesh) {
        // 有 Bug, 暂不能用
       
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K> Vb;
        typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo, K> Cb;
        typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
        typedef CGAL::Delaunay_triangulation_3<K, Tds> Triangulation_VoronoiMesh;
        typedef Triangulation_VoronoiMesh::Point Point_VoronoiMesh;

        voronoi_mesh.destory();

        Triangulation_VoronoiMesh* pt = new Triangulation_VoronoiMesh();

        Eigen::Vector3d BBoxMin(99999999, 99999999, 99999999), BBoxMax(-99999999, -99999999, -99999999);

        int vid = 0;
        for (auto& point : sample_points) {
            Triangulation_VoronoiMesh::Vertex_handle vh = pt->insert(Point_VoronoiMesh(point.x(), point.y(), point.z()));
            vh->info().id = vid;
            BBoxMin[0] = min(BBoxMin[0], point.x()); BBoxMin[1] = min(BBoxMin[1], point.y()); BBoxMin[2] = min(BBoxMin[2], point.z());
            BBoxMax[0] = max(BBoxMax[0], point.x()); BBoxMax[1] = max(BBoxMax[1], point.y()); BBoxMax[2] = max(BBoxMax[2], point.z());
        }
       /* for (auto vit = tmesh.vertices_begin(); vit != tmesh.vertices_end(); vit++, vid++) {
            auto& point = vit->point();
            Vertex_handle_t vh = pt->insert(Point_t(point.x(), point.y(), point.z()));
            vh->info().id = vid;
            BBoxMin[0] = min(BBoxMin[0], point.x()); BBoxMin[1] = min(BBoxMin[1], point.y()); BBoxMin[2] = min(BBoxMin[2], point.z());
            BBoxMax[0] = max(BBoxMax[0], point.x()); BBoxMax[1] = max(BBoxMax[1], point.y()); BBoxMax[2] = max(BBoxMax[2], point.z());
        }*/
        std::cout << "[Voronoi Mesh] Input_Points size = " << sample_points.size() << "\n";
        std::cout << "[Voronoi Mesh] Input_BBox: (" << BBoxMin[0] << "," << BBoxMin[1] << "," << BBoxMin[2] << ") --> (" <<
            BBoxMax[0] << "," << BBoxMax[1] << "," << BBoxMax[2] << ").\n";

        int fid = 0;
        for (Triangulation_VoronoiMesh::Finite_cells_iterator fci = pt->finite_cells_begin(); fci != pt->finite_cells_end(); fci++) {
            fci->info().id = fid++;
            vector<Eigen::Vector3d> points;
            for (int i = 0; i < 4; ++i) {
                auto p = fci->vertex(i)->point();
                points.emplace_back(Eigen::Vector3d(p.x(), p.y(), p.z()));
            }
            //auto center = CGAL::circumcenter(pt->tetrahedron(fci));
            //auto radius = pt->TetCircumRadius(pt->tetrahedron(fci));
            Eigen::Vector3d cent;
            double radius;
            CircumCenter::get_circum_center(points, cent, radius);
            fci->info().inside = true;
            if (cent.x() < BBoxMin.x() || cent.x() > BBoxMax.x()) {// 做最初步的判断
                fci->info().inside = false;
            } else if (cent.y() < BBoxMin.y() || cent.y() > BBoxMax.y()) {
                fci->info().inside = false;
            } else if (cent.z() < BBoxMin.z() || cent.z() > BBoxMax.z()) {
                fci->info().inside = false;
            }
        }

        int mas_vertex_count(0);
        int num_vor_v = 0;
        int num_vor_e = 0;
        int num_vor_f = 0;
        
        for (Triangulation_VoronoiMesh::Finite_cells_iterator fci = pt->finite_cells_begin(); fci != pt->finite_cells_end(); fci++) {
            if (fci->info().inside == false) {
                fci->info().tag = -1;
                continue;
            }
            fci->info().tag = mas_vertex_count++;

            vector<Eigen::Vector3d> points;
            for (int i = 0; i < 4; ++i) {
                auto p = fci->vertex(i)->point();
                points.emplace_back(Eigen::Vector3d(p.x(), p.y(), p.z()));
            }
            //auto center = CGAL::circumcenter(pt->tetrahedron(fci));
            //auto radius = pt->TetCircumRadius(pt->tetrahedron(fci));
            Eigen::Vector3d cent;
            double radius;
            CircumCenter::get_circum_center(points, cent, radius);
            voronoi_mesh.AddVertex(Vex(cent.x(), cent.y(), cent.z()));
            num_vor_v++;
        }

        for (Triangulation_VoronoiMesh::Finite_facets_iterator ffi = pt->finite_facets_begin(); ffi != pt->finite_facets_end(); ffi++) {
            //     Triangulation::Object o = pt->dual(*ffi);
            //     if (const Triangulation::Segment* s = CGAL::object_cast<Triangulation::Segment>(&o)) 		{
            //         if ((ffi->first->info().inside == false) || (pt->mirror_facet(*ffi).first->info().inside == false)) {
            //             continue;
            //         } 
            //         auto v1 = ffi->first->info().tag;
            //         auto v2 = pt->mirror_facet(*ffi).first->info().tag;
                     //num_vor_e++;// 曲线骨架的边
            //     }

            if ((ffi->first->info().inside == false) || (pt->mirror_facet(*ffi).first->info().inside == false)) {
                continue;
            }
            auto v1 = ffi->first->info().tag;
            auto v2 = pt->mirror_facet(*ffi).first->info().tag;
            num_vor_e++;// 曲线骨架的边
            voronoi_mesh.AddEdge(VH(v1), VH(v2));
        }

        for (Triangulation_VoronoiMesh::Finite_edges_iterator fei = pt->finite_edges_begin(); fei != pt->finite_edges_end(); fei++) {
            bool all_finite_inside = true;
            std::vector<Triangulation_VoronoiMesh::Cell_handle> vec_ch;
            Triangulation_VoronoiMesh::Cell_circulator cc = pt->incident_cells(*fei);
            do {
                if (pt->is_infinite(cc))
                    all_finite_inside = false;
                else if (cc->info().inside == false)
                    all_finite_inside = false;
                vec_ch.push_back(cc++);
            } while (cc != pt->incident_cells(*fei));

            if (!all_finite_inside) continue;

            for (unsigned k = 2; k < vec_ch.size() - 1; k++) {
                // 面片的边
                auto v1 = vec_ch[0]->info().tag;
                auto v2 = vec_ch[k]->info().tag;
            }

            for (unsigned k = 1; k < vec_ch.size() - 1; k++) {
                unsigned vid[3];
                vid[0] = vec_ch[0]->info().tag;
                vid[1] = vec_ch[k]->info().tag;
                vid[2] = vec_ch[k + 1]->info().tag;
                voronoi_mesh.AddFace({ VH(vid[0]), VH(vid[1]), VH(vid[2]) });
                num_vor_f++;
            }
        }
        
        //std::cout << "Medial Mesh: #V = " << num_vor_v << ", #E = " << num_vor_e << ", #F = " << num_vor_f << std::endl;

        std::cout << "[Voronoi Mesh]: #V = " << num_vor_v << ", #E = " << num_vor_e << ", #F = " << num_vor_f << std::endl;
	}

};