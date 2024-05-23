#include "FittingViaDicing.h"

FittingViaDicing::FittingViaDicing(MeshKernel::VolumeMesh& _hexmesh, 
	MeshKernel::SurfaceMesh& _refmesh, 
	double _Hausdorff_error_threshold):hexmesh(_hexmesh), refmesh(_refmesh) {

	hausdorff_error_threshold = _Hausdorff_error_threshold;
	hausdorff_error_threshold = 0.01;
	init();

}

void FittingViaDicing::init() {

	// 初始化 AABB树, 参照曲面上的采样点
	vector<Vector3f> vertices;
	for (auto& fp : refmesh.allfaces()) {
		const auto& vhs = fp.second.getVertexHandle();
		for (int i = 2; i < vhs.size(); ++i) {
			vector<VH> triangle = { vhs[0], vhs[i - 1], vhs[i] };
			for (auto& vh : triangle) {
				auto& v = refmesh.vertices(vh);
				Vector3f pos(v.x(), v.y(), v.z());
				vertices.emplace_back(pos);
			}
		}
		Vex center = refmesh.getFaceCenter(fp.first);
		sample_points_ref.emplace_back(center);
	}
	abtree_ref = new AABB_Tree(vertices);

	refmesh.initBBox();
	length_of_bbox = (refmesh.BBoxMax - refmesh.BBoxMin).norm();

}


void FittingViaDicing::fitting() {

	int iter = 0, iter_max = 5;
	

	while (iter < iter_max) {
		
		std::cout << "It " << iter++ << " : ";
		calc_hausdorff_dist_ref_to_hex();

		double dist_max = hausdorff_error_ref_to_hex.begin()->first;

		if (dist_max < hausdorff_error_threshold) {
			std::cout << "It " << iter << ": Hausdorff distance max = " << dist_max << ".It's adequate.\n";
			break;
		}

		unordered_set<FH> fhs_error;
		for (auto& dh : hausdorff_error_ref_to_hex) {
			if (dh.first < hausdorff_error_threshold) break;
			fhs_error.insert(dh.second);
		}
		
		for (auto& fh : fhs_error) {
			if (!hexmesh.isValid(fh)) continue;
			const auto& vhs = hexmesh.faces(fh).getVertexHandle();
			EH eh1 = hexmesh.getEdgeHandle(vhs[0], vhs[1]);
			EH eh2 = hexmesh.getEdgeHandle(vhs[1], vhs[2]);
			int cnt1 = INT_MAX, cnt2 = INT_MAX;
			if (HexMesh_TopologicalOperators::is_dicing_ok(hexmesh, eh1)) {
				unordered_map<VH, VH> left2right;
				HexMesh_TopologicalOperators::get_dicing_vhs(hexmesh, eh1, left2right);
				cnt1 = left2right.size();
			}
			if (HexMesh_TopologicalOperators::is_dicing_ok(hexmesh, eh2)) {
				unordered_map<VH, VH> left2right;
				HexMesh_TopologicalOperators::get_dicing_vhs(hexmesh, eh2, left2right);
				cnt2 = left2right.size();
			}
			if (cnt1 == INT_MAX && cnt2 == INT_MAX) continue;
			EH eh_dicing = (cnt1 < cnt2) ? eh1 : eh2;
			HexMesh_TopologicalOperators::dicing(hexmesh, eh_dicing, 2);
		}
		
		

		/*for (auto& vp : hexmesh.allvertices()) {
			if (vps_origin.count(vp.first)) continue;
			auto& v = hexmesh.vertices(vp.first);
			Vector3f pos_curr(v.x(), v.y(), v.z()), pos_nearest;
			abtree_ref->findNearstPoint(pos_curr, pos_nearest);
			v.setPosition(pos_nearest[0], pos_nearest[1], pos_nearest[2]);
		}*/
		hexmesh.updateAllHandles();

		VolumeFitting_CVIF app_cvif(hexmesh, refmesh);
		app_cvif.volume_fitting(1);

		//std::cout << "It " << iter++ << ": Hausdorff distance max = " << dist_max << ", the num of out-range faces = " << fhs_error.size() << std::endl;
		//break;
	}

	

	

}

void FittingViaDicing::calc_hausdorff_dist_ref_to_hex() {

	MeshKernel::SurfaceMesh trimesh_from_hexmesh;
	unordered_map<FH, FH> trifh_to_hexfh;
	TriMesh_Generator tg;
	tg.get_triangle_mesh(hexmesh, trimesh_from_hexmesh);
	tg.get_newfh_to_oldfh(trifh_to_hexfh);
	delete kdtree_hex;
	kdtree_hex = nullptr;
	kdtree_hex = new iGameKdTree(trimesh_from_hexmesh);
	hausdorff_error_ref_to_hex.clear();
	double dist_max = 0, dist_avg = 0;
	for (auto& pos : sample_points_ref) {
		auto nearest = kdtree_hex->nearest(pos);
		FH fh_hex = trifh_to_hexfh[nearest.fh];
		double dist = nearest.dist / length_of_bbox;
		dist_max = max(dist_max, dist);
		dist_avg += dist;
		hausdorff_error_ref_to_hex.emplace_back(dist, fh_hex);
	}
	dist_avg /= sample_points_ref.size();
	sort(hausdorff_error_ref_to_hex.begin(), hausdorff_error_ref_to_hex.end(), [&](auto& p1, auto& p2) {
		return p1.first > p2.first;
		});
	std::cout << "Hausdorff distance max = " << dist_max << ", avg = " << dist_avg << std::endl;
}