#include "ErrorRender.h"

namespace ErrorRender {

	void get_Hausdorff_distance(MeshKernel::SurfaceMesh& src_mesh, MeshKernel::VolumeMesh& ref_mesh, std::unordered_map<VH, Vector3d>& vh2color) {

		vh2color.clear();
		// 1. build AABB tree
		std::vector<Vector3f> aabb_vertices;

		for (auto& fp : ref_mesh.allfaces()) {
			if (!ref_mesh.isOnBoundary(fp.first)) continue;
			const auto& vhs = fp.second.getVertexHandle();
			for (int i = 2; i < vhs.size(); ++i) {
				std::vector<VH> tri_vhs = { vhs[0], vhs[i - 1], vhs[i] };
				for (auto& vh : tri_vhs) {
					auto& v = ref_mesh.vertices(vh);
					aabb_vertices.push_back(Vector3f(v.x(), v.y(), v.z()));
				}
			}
		}

		AABB_Tree aabb_tree(aabb_vertices);

		get_Hausdorff_distance(src_mesh, aabb_tree, vh2color);

	}

	void get_Hausdorff_distance(MeshKernel::SurfaceMesh& src_mesh, MeshKernel::SurfaceMesh& ref_mesh, std::unordered_map<VH, Vector3d>& vh2color) {

		vh2color.clear();
		// 1. build AABB tree
		std::vector<Vector3f> aabb_vertices;

		for (auto& fp : ref_mesh.allfaces()) {
			const auto& vhs = fp.second.getVertexHandle();
			for (int i = 2; i < vhs.size(); ++i) {
				std::vector<VH> tri_vhs = { vhs[0], vhs[i - 1], vhs[i] };
				for (auto& vh : tri_vhs) {
					auto& v = ref_mesh.vertices(vh);
					aabb_vertices.push_back(Vector3f(v.x(), v.y(), v.z()));
				}
			}
		}

		AABB_Tree aabb_tree(aabb_vertices);
		std::cout << "[Hausdorff Diatance Calculater]: Build AABB Tree success.\n";
		get_Hausdorff_distance(src_mesh, aabb_tree, vh2color);

	}

	

	void get_Hausdorff_distance(MeshKernel::VolumeMesh& src_mesh, MeshKernel::SurfaceMesh& ref_mesh, std::unordered_map<VH, Vector3d>& vh2color) {

		TriMesh_Generator tg;
		MeshKernel::SurfaceMesh src_trimesh;
		std::unordered_map<VH, VH> oldvh2newvh;
		std::unordered_map<VH, Vector3d> newvh2color;
		tg.get_triangle_mesh(src_mesh, src_trimesh);
		tg.get_oldvh_to_newvh(oldvh2newvh);
		std::cout << "[Hausdorff Diatance Calculater]: Generate triangles mesh success.\n";
		get_Hausdorff_distance(src_trimesh, ref_mesh, newvh2color);

		for (auto& vp : src_mesh.allvertices()) {
			if (oldvh2newvh.count(vp.first)) {
				vh2color[vp.first] = newvh2color[oldvh2newvh[vp.first]];
			}
		}

	}

	void get_Hausdorff_distance(MeshKernel::VolumeMesh& src_mesh, MeshKernel::VolumeMesh& ref_mesh, std::unordered_map<VH, Vector3d>& vh2color) {

		TriMesh_Generator tg;
		MeshKernel::SurfaceMesh ref_trimesh;
		tg.get_triangle_mesh(ref_mesh, ref_trimesh);
		get_Hausdorff_distance(src_mesh, ref_trimesh, vh2color);

	}

	void get_Hausdorff_distance(MeshKernel::SurfaceMesh& src_mesh, AABB_Tree& aabb_tree, std::unordered_map<VH, Vector3d>& vh2color) {

		ColorBar color_bar;

		//// 2. get vertex's min distance
		//std::unordered_map<VH, double> vh2dis;
		//double dis_max = 0, dis_avg = 0;
		//for (auto& vp : src_mesh.allvertices()) {
		//	Vector3f pos(vp.second.x(), vp.second.y(), vp.second.z());
		//	Vector3f nearest_p;
		//	double dis = aabb_tree.findNearstPoint(pos, nearest_p);
		//	//dis_min = std::min(dis_min, dis);
		//	dis_max = std::max(dis_max, dis);
		//	dis_avg += dis;
		//	vh2dis[vp.first] = dis;
		//}
		//dis_avg /= src_mesh.vsize();

		//// 3. set vertex's color
		//for (auto& vp : src_mesh.allvertices()) {
		//	double r = vh2dis[vp.first] / dis_max;
		//	vh2color[vp.first] = color_bar.getColor(r);
		//}

		// 在面片上采样
		std::unordered_map<FH, double> fh2dis;
		double dis_max = 0, dis_avg = 0;
		src_mesh.initBBox();
		double diagonal_length = (src_mesh.BBoxMax - src_mesh.BBoxMin).norm();
		for (auto& fp : src_mesh.allfaces()) {
			Vex face_center = src_mesh.getFaceCenter(fp.first);
			Vector3f pos(face_center.x(), face_center.y(), face_center.z());
			Vector3f nearest_p;
			double dis = aabb_tree.findNearstPoint(pos, nearest_p) / diagonal_length;
			dis_max = std::max(dis_max, dis);
			dis_avg += dis;
			fh2dis[fp.first] = dis;
		}
		dis_avg /= src_mesh.vsize();
		bool user_define_max_error_color = true;
		double dis_max_color = user_define_max_error_color ? 0.005 : dis_max;// 控制colorbar, 绝对误差
		std::unordered_map<FH, Vector3d> fh2color;
		for (auto& fp : src_mesh.allfaces()) {
			double r = fh2dis[fp.first] / dis_max_color;
			fh2color[fp.first] = color_bar.getColor(r);
		}

		for (auto& vp : src_mesh.allvertices()) {
			Vector3d color(0, 0, 0);
			const auto& adjfhs = src_mesh.NeighborFh(vp.first);
			for (auto& adjfh : adjfhs) {
				color += fh2color[adjfh];
			}
			color /= adjfhs.size();
			vh2color[vp.first] = color;
		}

		
		printf("[HausDist / DiagLen (relative)]: max = %.4f, avg = %.4f\n",
			dis_max, dis_avg);

	}

};