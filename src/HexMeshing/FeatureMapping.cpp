#include "FeatureMapping.h"

FeatureMapping::FeatureMapping(MeshKernel::SurfaceMesh& _refmesh, 
	MeshKernel::VolumeMesh& _hexmesh, 
	std::vector<int> _corners_ref,
	std::vector<int> _corners_hex): refmesh(_refmesh), hexmesh(_hexmesh) {

	this->corners_ref = _corners_ref;
	this->corners_hex = _corners_hex;

	for (int i = 0; i < corners_ref.size(); ++i) {
		corners_cons_hex2ref[corners_hex[i]] = corners_ref[i];
		corners_cons_ref2hex[corners_ref[i]] = corners_hex[i];
	}

}

void FeatureMapping::feature_line_detect(double angle) {

	double cosine_threshold = std::cos(angle / 180.0 * M_PI);
	std::cout << "cosine threshold = " << cosine_threshold << "\n";
	refmesh.genAllFacesNormal();
	for (auto& ep : refmesh.alledges()) {
		EH eh = ep.first;
		if (refmesh.isOnBoundary(eh)) {
			edge_features_ref.insert(eh);
		} else {
			auto adjfhs = refmesh.NeighborFh(eh);
			auto it = adjfhs.begin();
			FH fh1 = *it; ++it;
			FH fh2 = *it;
			Vex normal1 = refmesh.getFaceNormal(fh1);
			Vex normal2 = refmesh.getFaceNormal(fh2);
			double cosine = normal1.dot(normal2);
			if (cosine < cosine_threshold) {
				edge_features_ref.insert(eh);
			} else {
				//std::cout << "cosine = " << cosine << "\n";
			}
		}
	}
	std::cout << "[Feature Mapping]: " << edge_features_ref.size() << " edge features are found!!!\n";

}

void FeatureMapping::build_graph() {

	corner_features_ref.clear();
	corner_features_ref.insert(corners_ref.begin(), corners_ref.end());

	for (auto& eid : edge_features_ref) {

		if (eid2ecid.count(eid)) continue;// 已遍历过
		EH eh(eid);
		auto edge = refmesh.edges(eh);
		VH vh1 = edge.vh1();
		VH vh2 = edge.vh2();
		if (!corner_features_ref.count(vh1) && !corner_features_ref.count(vh2)) {
			continue;
		}
		VH vh_beg;
		if (corner_features_ref.count(vh1)) vh_beg = vh1;
		else vh_beg = vh2;

		std::unordered_set<int> visited_ehs;
		std::vector<int> chain;
		visited_ehs.insert(eh);
		chain.push_back(eh);
		VH vh_curr(vh1 + vh2 - vh_beg);
		while (!corner_features_ref.count(vh_curr)) {
			for (auto& adjeh : refmesh.NeighborEh(vh_curr)) {
				if (!edge_features_ref.count(adjeh) || visited_ehs.count(adjeh)) continue;
				auto adje = refmesh.edges(adjeh);
				VH vh_next(adje.vh1() + adje.vh2() - vh_curr);
				vh_curr = vh_next;
				visited_ehs.insert(adjeh);
				chain.push_back(adjeh);
				break;
			}
		}
		cc_lines_mapping[vh_beg][vh_curr] = cc_lines_mapping[vh_curr][vh_beg] = edges_chain.size();
		for (auto& eid : chain) {
			eid2ecid[eid] = edges_chain.size();
		}
		edges_chain.push_back(chain);

	}
	std::cout << "[Feature Mapping]: " << edges_chain.size() << " edges chain are found!!!\n";

}

void FeatureMapping::build_mapping(std::unordered_map<VH, int>& res_line_features_hex,
	std::vector<std::vector<EH>>& res_line_features_tri) {

	feature_line_detect(60);
	build_graph();

	auto shortest_path = [&](int src, int dest) {

		std::unordered_map<int, int> prev;
		std::queue<int> que;
		que.push(src);
		prev[src] = -1;
		while (!que.empty() && !prev.count(dest)) {
			auto vid = que.front();
			que.pop();
			VH vh(vid);
			for (auto& adjvh : hexmesh.NeighborVh(vh)) {
				if (prev.count(adjvh)) continue;
				prev[adjvh] = vh;
				que.push(adjvh);
			}
		}

		if (!prev.count(dest)) {
			std::cerr << "It is not connected between " << src << " <---> " << dest << std::endl;
			return std::vector<int>{};
		}

		std::vector<int> res;

		int vid = prev[dest];
		while (vid != src) {
			res.push_back(vid);
			vid = prev[vid];
		}

		return res;
	};

	res_line_features_tri.clear();
	for (auto& edge_c : edges_chain) {
		std::vector<EH> chain;
		for (auto& eid : edge_c) {
			chain.push_back(EH(eid));
		}
		res_line_features_tri.push_back(chain);
	}

	// 对于每个参照网格上的映射对, 找到六面体网格上对应的顶点链
	for (auto& ccp : cc_lines_mapping) {
		int vid1_ref = ccp.first;
		for (auto& vid_ecid : ccp.second) {
			int vid2_ref = vid_ecid.first;
			if (vid2_ref < vid1_ref) continue;// 防止重复遍历
			int vid1_hex = corners_cons_ref2hex[vid1_ref];
			int vid2_hex = corners_cons_ref2hex[vid2_ref];
			auto path = shortest_path(vid1_hex, vid2_hex);
			if (!path.empty()) {
				for (auto& vid_hex : path) {
					res_line_features_hex[VH(vid_hex)] = vid_ecid.second;
				}
			}
		}
	}

}