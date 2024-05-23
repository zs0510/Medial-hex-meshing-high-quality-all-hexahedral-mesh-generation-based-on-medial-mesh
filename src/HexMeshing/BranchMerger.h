#pragma once
#include "Skeletal Mesh/SkeletalMesh.h"

namespace SkeltalMesh_Operation {

	void branch_merge(SkeletalMesh& mesh) {

		std::unordered_map<EH, bool> curve_ehs;
		std::unordered_map<VH, bool> branch_vhs;
		std::unordered_map<VH, bool> joint_vhs;
		std::unordered_map<VH, bool> end_vhs;

		for (auto& ep : mesh.alledges()) {
			const auto& adjfhs = mesh.NeighborFh(ep.first);
			if (adjfhs.empty()) {
				curve_ehs[ep.first] = true;
			}
		}

		for (auto& vp : mesh.allvertices()) {

			VH vh(vp.first);

			const auto& adjfhs = mesh.NeighborFh(vh);
			const auto& adjehs = mesh.NeighborEh(vh);

			bool neighbor_with_curve = false;
			for (auto& adjeh : adjehs) {
				if (curve_ehs.count(adjeh)) {
					neighbor_with_curve = true;
					break;
				}
			}

			if (!neighbor_with_curve) continue;// 纯纯的三角面片顶点, 无需处理

			if (adjfhs.empty()) {// 只处理线线分支

				if (adjehs.size() > 2) {
					branch_vhs[vh] = true;// 线线分支
				} else if (adjehs.size() == 2) {
					joint_vhs[vh] = true;// 线上结点
				} else if (adjehs.size() == 1) {
					end_vhs[vh] = true;// 线端结点
				}
			}
			//} else {
			//	branch_vhs[vh] = true;// 线面分支
			//}


		}

		mesh.genAllEdgesLength();
		double edge_len_avg = 0.0;
		for (auto& ep : mesh.alledges()) {
			edge_len_avg += ep.second.getLength();
		}
		edge_len_avg /= mesh.esize();
		double merge_threshold = edge_len_avg * 0.35;

		for (auto& vp : branch_vhs) {
			if (!mesh.isValid(vp.first)) continue;
			const auto& adjvhs = mesh.NeighborVh(vp.first);
			auto& v = mesh.vertices(vp.first);
			for (auto& adjvh : adjvhs) {
				if (branch_vhs.count(adjvh)) {
					auto& adjv = mesh.vertices(adjvh);
					double len = (v - adjv).norm();
					if (len < merge_threshold) {
						const auto& _adjvhs = mesh.NeighborVh(adjvh);
						for (auto& _adjvh : _adjvhs) {
							if (_adjvh == vp.first) continue;
							mesh.AddEdge(_adjvh, vp.first);
						}
						v = (v + adjv) * 0.5;
						mesh.radius[vp.first] = (mesh.radius[vp.first] + mesh.radius[adjvh]) * 0.5;
						mesh.DeleteVertex(adjvh);
					}
				}
			}
		}
		mesh.updateAllHandles();
	}

};