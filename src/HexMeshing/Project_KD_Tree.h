#pragma once
#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/TriangleKdTree.h"
#include "pmp/algorithms/BarycentricCoordinates.h"
#include "Kernel/Mesh.h"
#include "Skeletal Mesh/SkeletalMesh.h"
#include "Tools/KdTree.h"

class Project_KD_Tree_2D_To_3D {

public:
	Project_KD_Tree_2D_To_3D(std::string filename, std::vector<std::pair<Vex, Vex>>& pos2_pos3) {

		mesh = new pmp::SurfaceMesh();
		mesh->read(filename);
		kd_tree = new pmp::TriangleKdTree(*mesh, 0);
		pro_pos3 = mesh->add_vertex_property<pmp::Point>("pos3");

		double dis_max = -1;

		for (auto vit = mesh->vertices_begin(); vit != mesh->vertices_end(); ++vit) {
			pmp::Point pos2_pmp = mesh->position(*vit);
			auto pos2_pmp_data = pos2_pmp.data();
			Vex pos2(pos2_pmp_data[0], pos2_pmp_data[1], pos2_pmp_data[2]);
			int nearest_id = -1;
			double nearest_dis = 999999;
			for (int i = 0; i < pos2_pos3.size(); ++i) {
				Vex& pos2_ref = pos2_pos3[i].first;
				double dis = (pos2 - pos2_ref).norm();
				if (dis < nearest_dis) {
					nearest_dis = dis;
					nearest_id = i;
				}
			}
			dis_max = std::max(nearest_dis, dis_max);
			if (nearest_id != -1) {
				Vex& pos3 = pos2_pos3[nearest_id].second;
				pro_pos3[*vit] = pmp::Point(pos3.x(), pos3.y(), pos3.z());
			}
		}
		std::cout << "Initialize Project_KD_Tree success. Distance(max) = " << dis_max << std::endl;
	}

	Vex get_project_pos3(Vex pos2);

private:

	pmp::TriangleKdTree* kd_tree = nullptr;
	pmp::SurfaceMesh* mesh = nullptr;
	pmp::VertexProperty<pmp::Point> pro_pos3;

};

Vex Project_KD_Tree_2D_To_3D::get_project_pos3(Vex pos2) {
	
	pmp::Point pos2_pmp(pos2.x(), pos2.y(), pos2.z());
	// find closest triangle of reference mesh
	pmp::TriangleKdTree::NearestNeighbor nn = kd_tree->nearest(pos2_pmp);
	const pmp::Point p = nn.nearest;
	const pmp::Face f = nn.face;

	// get face data
	pmp::SurfaceMesh::VertexAroundFaceCirculator fvIt = mesh->vertices(f);
	const pmp::Point p0 = mesh->position(*fvIt);
	const pmp::Point pro0 = pro_pos3[*fvIt];
	++fvIt;
	const pmp::Point p1 = mesh->position(*fvIt);
	const pmp::Point pro1 = pro_pos3[*fvIt];
	++fvIt;
	const pmp::Point p2 = mesh->position(*fvIt);
	const pmp::Point pro2 = pro_pos3[*fvIt];

	// get barycentric coordinates
	pmp::Point b = barycentric_coordinates(p, p0, p1, p2);

	// interpolate project result
	pmp::Point pro_pos3_pmp = pro0 * b[0] + pro1 * b[1] + pro2 * b[2];

	return Vex(pro_pos3_pmp[0], pro_pos3_pmp[1], pro_pos3_pmp[2]);

}

class Project_KD_Tree_SkeletalMesh {
public:
	Project_KD_Tree_SkeletalMesh(SkeletalMesh& skelmesh) {

		mesh = new pmp::SurfaceMesh();
		std::unordered_map<VH, pmp::Vertex> vh_v_mp;
		std::unordered_map<pmp::IndexType, std::unordered_set<pmp::IndexType>> to_vidx;// 记录 pmp 中每个顶点向外发射的边，以避免出现重复的半边
		std::unordered_map<pmp::IndexType, VH> v_vh_mp;
		// 三角形面片添加顶点的顺序只有两种
		for (auto& fp : skelmesh.allfaces()) {
			std::vector<pmp::Vertex> f_vs;
			for (auto& skl_vh : fp.second.getVertexHandle()) {
				if (!vh_v_mp.count(skl_vh)) {// 如果该顶点第一次被添加
					auto& v = skelmesh.vertices(skl_vh);
					pmp::Point pos(v.x(), v.y(), v.z());
					vh_v_mp[skl_vh] = mesh->add_vertex(pos);// sklmesh 到 pmp 的映射
					v_vh_mp[vh_v_mp[skl_vh].idx_] = skl_vh;// pmp 到 sklmesh 的映射
				}
				f_vs.push_back(vh_v_mp[skl_vh]);
			}
			for (int i = 0; i < 3; ++i) {// 检查有无重复边
				int j = (i + 1) % 3;
				auto vidx0 = f_vs[i].idx_;
				auto vidx1 = f_vs[j].idx_;
				if (to_vidx[vidx0].count(vidx1)) {// 重复边
					std::reverse(f_vs.begin(), f_vs.end());
					break;
				}
			}

			bool no_compliate_edge = true;
			for (int i = 0; i < 3; ++i) {// 再次检查有无重复边
				int j = (i + 1) % 3;
				auto vidx0 = f_vs[i].idx_;
				auto vidx1 = f_vs[j].idx_;
				if (to_vidx[vidx0].count(vidx1)) {// 重复边
					no_compliate_edge = false;
					std::cerr << "Find compliate edge. We choose ignore it.\n";
					break;
				}
			}

			if (no_compliate_edge) {
				for (int i = 0; i < 3; ++i) {
					int j = (i + 1) % 3;
					auto vidx0 = f_vs[i].idx_;
					auto vidx1 = f_vs[j].idx_;
					to_vidx[vidx0].insert(vidx1);
				}

				mesh->add_face(f_vs);
			}
			
		}

		kd_tree = new pmp::TriangleKdTree(*mesh, 0);
		pro_radius = mesh->add_vertex_property<double>("radius");
		
		for (auto vit = mesh->vertices_begin(); vit != mesh->vertices_end(); ++vit) {
			auto v = *vit;
			pro_radius[v] = skelmesh.radius[v_vh_mp[v.idx_]];
		}
		
	}

	double get_radius(Vex pos);
	Vex get_nearest(Vex pos);

private:
	pmp::TriangleKdTree* kd_tree = nullptr;
	pmp::SurfaceMesh* mesh = nullptr;
	pmp::VertexProperty<double> pro_radius;
	
};

double Project_KD_Tree_SkeletalMesh::get_radius(Vex pos) {

	pmp::Point pos_pmp(pos.x(), pos.y(), pos.z());
	// find closest triangle of reference mesh
	pmp::TriangleKdTree::NearestNeighbor nn = kd_tree->nearest(pos_pmp);
	const pmp::Point p = nn.nearest;
	const pmp::Face f = nn.face;

	// get face data
	pmp::SurfaceMesh::VertexAroundFaceCirculator fvIt = mesh->vertices(f);
	const pmp::Point p0 = mesh->position(*fvIt);
	double pro0 = pro_radius[*fvIt];
	++fvIt;
	const pmp::Point p1 = mesh->position(*fvIt);
	double pro1 = pro_radius[*fvIt];
	++fvIt;
	const pmp::Point p2 = mesh->position(*fvIt);
	double pro2 = pro_radius[*fvIt];

	// get barycentric coordinates
	pmp::Point b = barycentric_coordinates(p, p0, p1, p2);

	// interpolate project result
	double res = pro0 * b[0] + pro1 * b[1] + pro2 * b[2];

	return res;

}

Vex Project_KD_Tree_SkeletalMesh::get_nearest(Vex pos) {

	pmp::Point pos_pmp(pos.x(), pos.y(), pos.z());
	// find closest triangle of reference mesh
	pmp::TriangleKdTree::NearestNeighbor nn = kd_tree->nearest(pos_pmp);

	return Vex(nn.nearest[0], nn.nearest[1], nn.nearest[2]);
}

