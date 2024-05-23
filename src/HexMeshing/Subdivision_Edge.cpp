#include "Subdivision_Edge.h"

namespace Subdivision_Edge {

	void subdivision_edges(MeshKernel::VolumeMesh& hexmesh, EH input_eh, int num_of_segments) {

		if (!hexmesh.isValid(input_eh) || num_of_segments <= 1) return;

		std::unordered_map<VH, VH> left2right;// �з���ؼ�¼һ���ߵ���������handle
		std::unordered_set<VH> total_vhs;// ��¼���ж��� handle

		std::unordered_map<FH, bool> visited_fh;
		std::unordered_map<EH, bool> visited_eh;
		std::queue<EH> que;
		que.emplace(input_eh);
		VH vh_left = hexmesh.edges(input_eh).vh1(), vh_right = hexmesh.edges(input_eh).vh2();
		left2right[vh_left] = vh_right;
		visited_eh[input_eh] = true;
		total_vhs.insert(vh_left);
		total_vhs.insert(vh_right);

		while (!que.empty()) {// ʹ�� BFS, ��������Ҫϸ�ֵı��ҵ�
			auto cur_eh = que.front();
			que.pop();
			auto& edge = hexmesh.edges(cur_eh);
			VH vh1 = edge.vh1(), vh2 = edge.vh2();
			const auto& adjfhs = hexmesh.NeighborFh(cur_eh);
			for (auto& adjfh : adjfhs) {
				if (visited_fh.count(adjfh)) continue;
				visited_fh[adjfh] = true;
				auto& face = hexmesh.faces(adjfh);
				const auto& f_ehs = face.getEdgeHandle();
				for (auto& f_eh : f_ehs) {
					if (visited_eh.count(f_eh)) continue;
					auto& f_edge = hexmesh.edges(f_eh);
					if (f_edge.vh1() == vh1 || f_edge.vh2() == vh1
						|| f_edge.vh1() == vh2 || f_edge.vh2() == vh2) continue;// ��ͬһ�������Ҳ����ڵı�

					vh_left = f_edge.vh1(), vh_right = f_edge.vh2();
					bool vh_left_is_ok = false;
					for (auto& adjvh : hexmesh.NeighborVh(vh_left)) {
						if (left2right.count(adjvh)) {
							vh_left_is_ok = true;
							break;
						}
					}
					if (!vh_left_is_ok) std::swap(vh_left, vh_right);// ��������, ȷ��ϸ��ʱ��������Խ�
					que.emplace(f_eh);
					left2right[vh_left] = vh_right;
					visited_eh[f_eh] = true;
					total_vhs.insert(vh_left);
					total_vhs.insert(vh_right);
				}
			}
		}

		std::unordered_map<VH, std::vector<VH>> left2vhs;// ��¼ϸ������Ҫ�õ��Ķ��� handle
		for (auto& vp : left2right) {
			auto& v_left = hexmesh.vertices(vp.first);
			auto& v_right = hexmesh.vertices(vp.second);
			Vec vec = v_right - v_left;

			left2vhs[vp.first].emplace_back(vp.first);
			for (int i = 1; i < num_of_segments; ++i) {
				left2vhs[vp.first].emplace_back( hexmesh.AddVertex(v_left + vec * (i * 1.0 / num_of_segments)));
			}
			left2vhs[vp.first].emplace_back(vp.second);

		}


		std::unordered_map<CH, bool> visited_ch;
		for (auto& ep : visited_eh) {
			auto eh = ep.first;
			const auto& adjchs = hexmesh.NeighborCh(eh);
			for (auto& adjch : adjchs) {
				if (visited_ch.count(adjch)) continue;
				visited_ch[adjch] = true;

				bool is_involved = true;// ����Ƿ�����Ҫ��ϸ�ֵ���
				auto& cell = hexmesh.cells(adjch);
				const auto& c_vhs = cell.getVertexHandle();
				for (auto& c_ch : c_vhs) {
					if (!total_vhs.count(c_ch)) {
						is_involved = false;
						break;
					}
				}
				if (!is_involved) continue;

				const auto& c_fhs = cell.getFaceHandle();// �ҵ�����Ǹ���
				for (auto& c_fh : c_fhs) {
					auto& face = hexmesh.faces(c_fh);
					const auto& f_vhs = face.getVertexHandle();
					bool is_left_face = true;
					for (auto& f_vh : f_vhs) {
						if (!left2right.count(f_vh)) {
							is_left_face = false;
							break;
						}
					}
					if (is_left_face) {
						for (int i = 1; i <= num_of_segments; ++i) {
							int j = i - 1;
							hexmesh.AddCell( { left2vhs[f_vhs[0]][j], left2vhs[f_vhs[1]][j] ,left2vhs[f_vhs[2]][j] ,left2vhs[f_vhs[3]][j],
								left2vhs[f_vhs[0]][i], left2vhs[f_vhs[1]][i] ,left2vhs[f_vhs[2]][i] ,left2vhs[f_vhs[3]][i] });
						}
						break;// һ����ֻ��һ������, ֻ��ϸ��һ��
					}
				}


			}
		}

		for (auto& ep : visited_eh) {// ɾȥ��ϸ�ֵı�
			if (hexmesh.isValid(ep.first)) {
				hexmesh.DeleteEdge(ep.first);
			}
		}

		//hexmesh.updateAllHandles();

	}

	void subdivision_edges_weighted(MeshKernel::VolumeMesh& hexmesh, EH input_eh, double weight_vh1) {

		if (!hexmesh.isValid(input_eh) || weight_vh1 <= 0.0 || weight_vh1 >= 1.0) return;

		std::unordered_map<VH, VH> left2right;// �з���ؼ�¼һ���ߵ���������handle
		std::unordered_set<VH> total_vhs;// ��¼���ж��� handle

		std::unordered_map<FH, bool> visited_fh;
		std::unordered_map<EH, bool> visited_eh;
		std::queue<EH> que;
		que.emplace(input_eh);
		VH vh_left = hexmesh.edges(input_eh).vh1(), vh_right = hexmesh.edges(input_eh).vh2();
		left2right[vh_left] = vh_right;
		visited_eh[input_eh] = true;
		total_vhs.insert(vh_left);
		total_vhs.insert(vh_right);

		while (!que.empty()) {// ʹ�� BFS, ��������Ҫϸ�ֵı��ҵ�
			auto cur_eh = que.front();
			que.pop();
			auto& edge = hexmesh.edges(cur_eh);
			VH vh1 = edge.vh1(), vh2 = edge.vh2();
			const auto& adjfhs = hexmesh.NeighborFh(cur_eh);
			for (auto& adjfh : adjfhs) {
				if (visited_fh.count(adjfh)) continue;
				visited_fh[adjfh] = true;
				auto& face = hexmesh.faces(adjfh);
				const auto& f_ehs = face.getEdgeHandle();
				for (auto& f_eh : f_ehs) {
					if (visited_eh.count(f_eh)) continue;
					auto& f_edge = hexmesh.edges(f_eh);
					if (f_edge.vh1() == vh1 || f_edge.vh2() == vh1
						|| f_edge.vh1() == vh2 || f_edge.vh2() == vh2) continue;// ��ͬһ�������Ҳ����ڵı�

					vh_left = f_edge.vh1(), vh_right = f_edge.vh2();
					bool vh_left_is_ok = false;
					for (auto& adjvh : hexmesh.NeighborVh(vh_left)) {
						if (left2right.count(adjvh)) {
							vh_left_is_ok = true;
							break;
						}
					}
					if (!vh_left_is_ok) std::swap(vh_left, vh_right);// ��������, ȷ��ϸ��ʱ��������Խ�
					que.emplace(f_eh);
					left2right[vh_left] = vh_right;
					visited_eh[f_eh] = true;
					total_vhs.insert(vh_left);
					total_vhs.insert(vh_right);
				}
			}
		}
		double weight_vh2 = 1.0 - weight_vh1;
		std::unordered_map<VH, std::vector<VH>> left2vhs;// ��¼ϸ������Ҫ�õ��Ķ��� handle
		for (auto& vp : left2right) {
			auto& v_left = hexmesh.vertices(vp.first);
			auto& v_right = hexmesh.vertices(vp.second);
			Vex v_new = v_left * weight_vh1 + v_right * weight_vh2;
			left2vhs[vp.first].emplace_back(vp.first);
			left2vhs[vp.first].emplace_back(hexmesh.AddVertex(v_new));
			left2vhs[vp.first].emplace_back(vp.second);
		}

		std::unordered_map<CH, bool> visited_ch;
		for (auto& ep : visited_eh) {
			auto eh = ep.first;
			const auto& adjchs = hexmesh.NeighborCh(eh);
			for (auto& adjch : adjchs) {
				if (visited_ch.count(adjch)) continue;
				visited_ch[adjch] = true;

				bool is_involved = true;// ����Ƿ�����Ҫ��ϸ�ֵ���
				auto& cell = hexmesh.cells(adjch);
				const auto& c_vhs = cell.getVertexHandle();
				for (auto& c_ch : c_vhs) {
					if (!total_vhs.count(c_ch)) {
						is_involved = false;
						break;
					}
				}
				if (!is_involved) continue;

				const auto& c_fhs = cell.getFaceHandle();// �ҵ�����Ǹ���
				for (auto& c_fh : c_fhs) {
					auto& face = hexmesh.faces(c_fh);
					const auto& f_vhs = face.getVertexHandle();
					bool is_left_face = true;
					for (auto& f_vh : f_vhs) {
						if (!left2right.count(f_vh)) {
							is_left_face = false;
							break;
						}
					}
					if (is_left_face) {

						for (int i = 1; i <= 2; ++i) {
							int j = i - 1;
							hexmesh.AddCell({ left2vhs[f_vhs[0]][j], left2vhs[f_vhs[1]][j] ,left2vhs[f_vhs[2]][j] ,left2vhs[f_vhs[3]][j],
								left2vhs[f_vhs[0]][i], left2vhs[f_vhs[1]][i] ,left2vhs[f_vhs[2]][i] ,left2vhs[f_vhs[3]][i] });
						}
						break;// һ����ֻ��һ������, ֻ��ϸ��һ��
					}
				}


			}
		}

		for (auto& ep : visited_eh) {// ɾȥ��ϸ�ֵı�
			if (hexmesh.isValid(ep.first)) {
				hexmesh.DeleteEdge(ep.first);
			}
		}


	}

};