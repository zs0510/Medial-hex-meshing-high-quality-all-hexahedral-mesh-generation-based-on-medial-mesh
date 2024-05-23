#pragma once
#include "pmp/SurfaceMesh.h"
#include "Kernel/Mesh.h"

#include "pmp/algorithms/SurfaceRemeshing.h"

namespace MeshConverter {

	pmp::SurfaceMesh remeshing_pmp(MeshKernel::SurfaceMesh& mesh, std::vector<int>& _selected_fhs, double uniform_len_ratio = 0.75, bool uniform_flag = true) {

		std::unordered_map<FH, bool> selected_fhs;
		for (auto& _fh : _selected_fhs) {
			selected_fhs[FH(_fh)] = true;
		}

		std::unordered_map<VH, bool> selected_vhs;
		for (auto& vp : mesh.allvertices()) {// 如果一个顶点周围所有的面都被选择，我们则认为这个顶点也被选择
			bool selected = true;
			for (auto& adjfh : mesh.NeighborFh(vp.first)) {
				if (!selected_fhs.count(adjfh)) {
					selected = false;
					break;
				}
			}
			if (selected) {
				selected_vhs[vp.first] = true;
			}
		}

		std::unordered_set<EH> selected_ehs;
		for (auto& fp : mesh.allfaces()) {
			auto& f = mesh.faces(fp.first);
			const auto& ehs = f.getEdgeHandle();
			for (auto& eh : ehs) {
				selected_ehs.insert(eh);
			}
		}


		pmp::SurfaceMesh pmp_mesh;
		pmp_mesh.add_vertex_property<bool>("v:selected", false);
		auto selected_pmp = pmp_mesh.get_vertex_property<bool>("v:selected");

		std::unordered_map<VH, pmp::Vertex> igame_to_pmp;
		for (auto& vp : mesh.allvertices()) {
			auto& vh = vp.first;
			auto& pos = vp.second;
			pmp::Point point(pos.x(), pos.y(), pos.z());
			auto vit = pmp_mesh.add_vertex(point);
			if (selected_vhs.count(vh)) {
				selected_pmp[vit] = true;// 标记为待remeshing的点
			}
			igame_to_pmp[vh] = vit;
		}

		for (auto& fp : mesh.allfaces()) {
			const auto& f_vhs = fp.second.getVertexHandle();
			std::vector<pmp::Vertex> vhs;
			for (auto& f_vh : f_vhs) {
				vhs.push_back(igame_to_pmp[f_vh]);
			}
			pmp_mesh.add_face(vhs);
		}

		printf("iGame Mesh convert to PMP Mesh success. vsize = %d, esize = %d, fsize = %d\n", 
			pmp_mesh.vertices_size(), pmp_mesh.edges_size(), pmp_mesh.faces_size());

		double len_min = 999999, len_max = -1, len_avg = 0;

		if (selected_ehs.empty()) {
			for (auto eit = pmp_mesh.edges_begin(); eit != pmp_mesh.edges_end(); eit++) {
				double len = pmp_mesh.edge_length(*eit);
				len_min = std::min(len_min, len);
				len_max = std::max(len_max, len);
				len_avg += len;
			}
			len_avg /= pmp_mesh.edges_size();
		} else {
			mesh.genAllEdgesLength();
			for (auto& eh : selected_ehs) {
				auto& e = mesh.edges(eh);
				double len = e.getLength();
				len_min = std::min(len_min, len);
				len_max = std::max(len_max, len);
				len_avg += len;
			}
			len_avg /= selected_ehs.size();
		}

		

		

		pmp::SurfaceRemeshing remeshing_app(pmp_mesh);

		if (uniform_flag) {
			remeshing_app.uniform_remeshing(len_avg * uniform_len_ratio);
			//pmp_mesh.write("C:/My Files/Graphics/model_data/uniform_remeshing_result.obj");
			std::cout << "We are running uniform remeshing.\n";
		} else {
			remeshing_app.adaptive_remeshing(len_min * 0.25, len_max * 4, len_min * 0.1);// 最小边长，最大边长，最大允许误差
			//pmp_mesh.write("C:/My Files/Graphics/model_data/adaptive_remeshing_result.obj");
			std::cout << "We are running adaptive remeshing.\n";
		}

		mesh.destory();
		pmp_mesh.garbage_collection();
		std::unordered_map<int, VH> pmp_to_igame;
		for (auto vit = pmp_mesh.vertices_begin(); vit != pmp_mesh.vertices_end(); ++vit) {
			int vidx = (*vit).idx();
			auto pos = pmp_mesh.position(*vit);
			Vex vex(pos.data()[0], pos.data()[1], pos.data()[2]);
			pmp_to_igame[vidx] = mesh.AddVertex(vex);
		}
		for (auto fit = pmp_mesh.faces_begin(); fit != pmp_mesh.faces_end(); ++fit) {
			pmp::Halfedge h = pmp_mesh.halfedge(*fit);
			pmp::Halfedge hend = h;
			std::vector<VH> vhs;
			do {
				int vidx = pmp_mesh.to_vertex(h).idx_;
				vhs.push_back(pmp_to_igame[vidx]);
				h = pmp_mesh.next_halfedge(h);
			} while (h != hend);

			mesh.AddFace(vhs);
		}

		_selected_fhs.clear();

		return pmp_mesh;

	}


};