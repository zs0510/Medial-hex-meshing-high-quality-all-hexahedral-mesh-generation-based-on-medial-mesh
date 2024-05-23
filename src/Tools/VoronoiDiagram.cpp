#include "VoronoiDiagram.h"

void VoronoiDiagram::get_vertices_of_voronoi_diagram(const vector<Vector3d>& sample_points, vector<Vector3d>& vertices_of_voronoi_diagram) {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Regular_triangulation_3<K> Regular_triangulation;
    typedef K::Vector_3                                     Vector;
    typedef K::Point_3 Point;
    typedef K::Weighted_point_3 Wp;
    typedef K::Point_3                                     Point;
    typedef CGAL::Regular_triangulation_3<K>               Rt;

    vertices_of_voronoi_diagram.clear();

    vector<Regular_triangulation::Weighted_point> wpoints;
    vector<Point> points;
    vector<double> X, Y, Z;
    set<Point> Voronoi_vert;

    // �������в�����
    for (auto& p : sample_points) {
        double x = p[0];
        double y = p[1];
        double z = p[2];
        wpoints.push_back(Wp(Point(x, y, z), 0));
        points.push_back(Point(x, y, z));
        X.push_back(x);
        Y.push_back(y);
        Z.push_back(z);
    }

    Regular_triangulation rt(wpoints.begin(), wpoints.end());
    rt.is_valid();
    Regular_triangulation::Edge_iterator eit;
    Regular_triangulation::Point_iterator pit;

    // ���� Delaunay �����бߣ����� Delaunay ͼ�Ķ�żͼ���� Voronoi ͼ
    /*cout << "====all points in rt====" << endl;
    for (auto vit = rt.vertices_begin(); vit != rt.vertices_end(); ++vit) {
        cout << vit->point() << endl;

    }
    cout << "==== Voronoi vert ====" << endl;*/
    Rt::Finite_edges_iterator fei;
    for (fei = rt.finite_edges_begin(); fei != rt.finite_edges_end(); fei++) {
        // ����ÿ����
        auto face = rt.incident_facets(*fei);
        auto face_start = face;
        // �ҵ�������ڵ��棨��׼ȷ�ģ�����ÿ���ڽ��棩
        do {
            CGAL::Object  dualedge = rt.dual(*face);
            // �ҵ����Ӧ�ı�
            if (CGAL::object_cast<Rt::Segment>(&dualedge)) {
                Point source = CGAL::object_cast<Rt::Segment>(&dualedge)->source();
                Point target = CGAL::object_cast<Rt::Segment>(&dualedge)->target();
                Voronoi_vert.insert(source);
                Voronoi_vert.insert(target);
            } else if (CGAL::object_cast<Rt::Ray>(&dualedge)) {// ��������������ߣ����������
                Point target = CGAL::object_cast<Rt::Ray>(&dualedge)->source();
                Voronoi_vert.insert(target);
            }
            ++face;
        } while (face_start != face);
    }

    set<Point>::iterator iter = Voronoi_vert.begin();
    while (iter != Voronoi_vert.end()) {
        //cout << *iter << endl;
        Vector3d pos((*iter).x(), (*iter).y(), (*iter).z());
        vertices_of_voronoi_diagram.emplace_back(pos);
        ++iter;
    }

    cout << "Input sample points size = " << sample_points.size()
        << ", output voronoi diagram points size = " << vertices_of_voronoi_diagram.size() << std::endl;
}

void VoronoiDiagram::get_medial_mesh(const vector<Vector3d>& sample_points, SkeletalMesh& medial_mesh) {

    medial_mesh.destory();

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Regular_triangulation_3<K> Regular_triangulation;
    typedef K::Vector_3                                     Vector;
    typedef K::Point_3 Point;
    typedef K::Weighted_point_3 Wp;
    typedef K::Point_3                                     Point;
    typedef CGAL::Regular_triangulation_3<K>               Rt;

    vector<Regular_triangulation::Weighted_point> wpoints;
    vector<Point> points;
    vector<double> X, Y, Z;
    unordered_map<Point, VH> Voronoi_vert;

    // �������в�����
    for (auto& p : sample_points) {
        double x = p[0];
        double y = p[1];
        double z = p[2];
        wpoints.push_back(Wp(Point(x, y, z), 0));
        points.push_back(Point(x, y, z));
        X.push_back(x);
        Y.push_back(y);
        Z.push_back(z);
    }

    Regular_triangulation rt(wpoints.begin(), wpoints.end());
    rt.is_valid();
    Regular_triangulation::Edge_iterator eit;
    Regular_triangulation::Point_iterator pit;

    // ���� Delaunay �����бߣ����� Delaunay ͼ�Ķ�żͼ���� Voronoi ͼ
    /*cout << "====all points in rt====" << endl;
    for (auto vit = rt.vertices_begin(); vit != rt.vertices_end(); ++vit) {
        cout << vit->point() << endl;

    }
    cout << "==== Voronoi vert ====" << endl;*/
    Rt::Finite_edges_iterator fei;
    for (fei = rt.finite_edges_begin(); fei != rt.finite_edges_end(); fei++) {
        // ����ÿ����
        auto face = rt.incident_facets(*fei);
        auto face_start = face;
        // �ҵ�������ڵ��棨��׼ȷ�ģ�����ÿ���ڽ��棩
        do {
            CGAL::Object  dualedge = rt.dual(*face);
            // �ҵ����Ӧ�ı�
            if (CGAL::object_cast<Rt::Segment>(&dualedge)) {
                Point source = CGAL::object_cast<Rt::Segment>(&dualedge)->source();
                Point target = CGAL::object_cast<Rt::Segment>(&dualedge)->target();
                if (!Voronoi_vert.count(source)) {
                    Vex v(source.x(), source.y(), source.z());
                    Voronoi_vert[source] = medial_mesh.AddVertex(v);
                }
                if (!Voronoi_vert.count(target)) {
                    Vex v(target.x(), target.y(), target.z());
                    Voronoi_vert[target] = medial_mesh.AddVertex(v);
                }
                medial_mesh.AddEdge(Voronoi_vert[source], Voronoi_vert[target]);
            } else if (CGAL::object_cast<Rt::Ray>(&dualedge)) {// ��������������ߣ����������
                Point source = CGAL::object_cast<Rt::Ray>(&dualedge)->source();
                if (!Voronoi_vert.count(source)) {
                    Vex v(source.x(), source.y(), source.z());
                    Voronoi_vert[source] = medial_mesh.AddVertex(v);
                }
            }
            ++face;
        } while (face_start != face);
    }

    std::cout << "Compute medial mesh success. Vertices size = " << medial_mesh.vsize()
        << ", edges size = " << medial_mesh.esize() << std::endl;

}