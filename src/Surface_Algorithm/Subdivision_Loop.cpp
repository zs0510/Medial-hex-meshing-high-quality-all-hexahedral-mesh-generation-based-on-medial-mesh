#include "Subdivision_Loop.h"

#ifndef PI
#define PI 3.14159265358979f
#endif // !PI

typedef MeshKernel::iGameVertexHandle VH;

void Subdivision_Loop::Execute() {

	while (iter--) {

		auto ori_faces = mesh.allfaces();
		auto ori_edges = mesh.alledges();
		auto ori_vertices = mesh.allvertices();
		std::vector<VH> e2v(mesh.EdgeSize());// origin edge to new vertex
		std::vector<Eigen::Vector3d> new_pos(mesh.VertexSize());

		// calculate origin vertex's new position
		for (auto& vp : ori_vertices) {
			auto vh = vp.first;
			auto v = vp.second;
			Eigen::Vector3d pos(v.x(), v.y(), v.z());
			auto adj_v = mesh.NeighborVh(vh);
			auto adj_sz = adj_v.size();
			if (adj_sz == 0) continue;
			if (mesh.isOnBoundary(vh)) {
				Eigen::Vector3d tmp = Eigen::Vector3d::Zero();
				for (auto adjvh : adj_v) {
					auto adjv = mesh.vertices(adjvh);
					tmp += (Eigen::Vector3d(adjv.x(), adjv.y(), adjv.z()));
				}
				tmp /= adj_sz;
				pos = pos * 0.75 + tmp * 0.25;
			} else {
				double a = 2 * std::pow((0.375 + 0.25 * std::cos(2.f * PI / adj_sz)), 2);// - 0.25f?
				double b = (1 - a) / adj_sz;
				pos *= a;
				for (auto adjvh : adj_v) {
					auto adjv = mesh.vertices(adjvh);
					pos += (Eigen::Vector3d(adjv.x(), adjv.y(), adjv.z()) * b);
				}
			}
			new_pos[vh] = pos;
		}

		// generate new vertex
		for (auto& ep : ori_edges) {
			auto eh = ep.first;
			auto e = ep.second;
			auto vh1 = e.vh1(), vh2 = e.vh2();
			auto v1 = mesh.vertices(vh1);
			auto v2 = mesh.vertices(vh2);
			VH vh(-1);
			Eigen::Vector3d pos(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
			if (mesh.isOnBoundary(eh)) {
				pos /= 2;
			} else {
				// 寻找边的两个对点
				MeshKernel::iGameVertexHandle vh3(-1), vh4(-1);
				for (auto fh : mesh.NeighborFh(eh)) {
					auto vhs = (mesh.faces(fh)).getVertexHandle();
					assert(vhs.size() > 2);
					int tmp = vhs[0] + vhs[1] + vhs[2] - vh1 - vh2;
					if (vh3 == -1) vh3 = MeshKernel::iGameVertexHandle(tmp);
					else vh4 = MeshKernel::iGameVertexHandle(tmp);
				}
				if (vh3 == -1 || vh4 == -1) continue;
				auto v3 = mesh.vertices(vh3);
				auto v4 = mesh.vertices(vh4);
				pos *= 0.375;
				Eigen::Vector3d tmp(v3.x() + v4.x(), v3.y() + v4.y(), v3.z() + v4.z());
				tmp *= 0.125;
				pos += tmp;

			}
			MeshKernel::iGameVertex vex(pos[0], pos[1], pos[2]);
			vh = mesh.AddVertex(vex);
			e2v[eh] = vh;
		}

		// update all faces
		for (auto& fp : ori_faces) {
			auto& face = mesh.faces(fp.first);
			auto ehs = face.getEdgeHandle();
			auto vhs = face.getVertexHandle();
			assert(ehs.size() == 3 && vhs.size() == 3);
			mesh.AddFace({ VH(e2v[ehs[0]]), VH(e2v[ehs[1]]), VH(e2v[ehs[2]]) });
			for (int i = 0; i < 3; ++i) {
				auto cur = vhs[i];
				auto pre = vhs[(i + 2) % 3];
				auto next = vhs[(i + 1) % 3];
				VH new_pre(-1), new_next(-1);
				for (auto eh : ehs) {
					auto e = mesh.edges(eh);
					int vex = e.vh1() + e.vh2() - cur;
					if (new_pre == -1 && vex == pre) {
						new_pre = e2v[eh];
					}
					if (new_next == -1 && vex == next) {
						new_next = e2v[eh];
					}
				}
				assert(new_pre != -1 && new_next != -1);
				mesh.AddFace({ new_pre, cur, new_next });
			}
		}

		// delete all origin edges
		for (auto& ep : ori_edges) {
			auto eh = ep.first;
			if (mesh.isValid(eh)) {
				mesh.DeleteEdge(eh);
			}
		}

		// update origin vertex position
		for (auto& vp : ori_vertices) {
			auto vh = vp.first;
			auto& v = mesh.vertices(vh);
			v.setPosition(new_pos[vh][0], new_pos[vh][1], new_pos[vh][2]);
		}

		mesh.updateAllHandles();
	}
}