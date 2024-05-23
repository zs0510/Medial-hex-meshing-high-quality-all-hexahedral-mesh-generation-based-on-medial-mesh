#include "DomainProcessor.h"

using namespace std;

void DomainProcessor::detect_faces(SkeletalMesh& mesh, vector<vector<VH>>& domains) {

	unordered_map<VH, bool> visited;

	for (int i = 0; i < mesh.esize(); ++i) {
		EH eh(i);
		if (!mesh.isValid(eh)) continue;
		auto& edge = mesh.edges(eh);
		if (mesh.NeighborFh(eh).size() != 1) continue;// 非边界边直接返回
		auto vh1 = edge.vh1();
		auto vh2 = edge.vh2();
		if (visited.count(vh1) && visited.count(vh2)) continue;
		VH curvh = visited.count(vh1) ? vh2 : vh1;
		vector<VH> domain;
		bool is_valid = true;
		VH nextvh(-1), prevh(-1);
		do {
			visited[curvh] = true;
			for (auto& adjeh : mesh.NeighborEh(curvh)) {
				if (mesh.NeighborFh(adjeh).size() == 1) {
					auto& adje = mesh.edges(adjeh);
					VH adjvh(adje.vh1() + adje.vh2() - curvh);
					if (adjvh != prevh) {
						nextvh = adjvh;
						break;
					}
				}
			}
			if (nextvh == -1) {
				std::cerr << "Error: Can't find a valid domain.\n";
				is_valid = false;
				break;
			}
			domain.push_back(curvh);
			prevh = curvh;
			curvh = nextvh;
		} while (curvh != domain[0]);

		if (is_valid) {
			domains.push_back(domain);
		}
	}


	std::cout << "Find " << domains.size() << " domains.\n";

}

void DomainProcessor::get_input_domains(SkeletalMesh& mesh, vector<vector<VH>>& domains, vector<vector<Vex>>& expanded_domains) {
	expanded_domains.clear();
	// 将三角形面片的边界扩张
	for (auto& domain : domains) {
		vector<Vex> new_domain_points;
		
		for (int i = 0; i < domain.size(); ++i) {
			VH curvh = domain[i];
			VH prevh = domain[(i == 0) ? domain.size() - 1 : (i - 1)];
			VH nextvh = domain[(i + 1) % domain.size()];
			EH eh0 = mesh.getEdgeHandle(curvh, prevh);
			EH eh1 = mesh.getEdgeHandle(curvh, nextvh);
			if (mesh.NeighborFh(eh0).empty() || mesh.NeighborFh(eh1).empty()) {
				std::cerr << "\nError: Cann't find adjoining face...";
			}/* else {
				std::cout << "\n";
			}*/
			FH fh0 = (*(mesh.NeighborFh(eh0)).begin());
			FH fh1 = (*(mesh.NeighborFh(eh1)).begin());
			/*printf("curvh = %d, prevh = %d, nextvh = %d, eh0 = %d, eh1 = %d, fh0 = %d, fh1 = %d\n",
				curvh, prevh, nextvh, eh0, eh1, fh0, fh1);*/
			auto& curv = mesh.vertices(curvh);
			auto& prev = mesh.vertices(prevh);
			auto& nextv = mesh.vertices(nextvh);
			Vec vec0 = (curv - prev).normalized();
			Vec vec1 = (curv - nextv).normalized();
			Vec n_face0 = mesh.getFaceNormal(fh0);
			Vec n_face1 = mesh.getFaceNormal(fh1);
			Vec dir0 = (vec0 % n_face0).normalized();// dir 与 边向量、面法向量正交
			Vec dir1 = (vec1 % n_face1).normalized();
			Vex c_face0 = mesh.getFaceCenter(fh0);
			Vex c_face1 = mesh.getFaceCenter(fh1);
			Vex mid_e0 = mesh.getEdgeMidpoint(eh0);
			Vex mid_e1 = mesh.getEdgeMidpoint(eh1);
			Vec vec_mc0 = (c_face0 - mid_e0).normalized();// 朝内的方向
			Vec vec_mc1 = (c_face1 - mid_e1).normalized();
			/*printf("Normals' norm: %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n", 
				vec0.norm(), vec1.norm(), n_face0.norm(), n_face1.norm(), dir0.norm(), dir1.norm(), vec_mc0.norm(), vec_mc1.norm());*/
			if (dir0.dot(vec_mc0) > 0) dir0 *= -1;// dir 须朝外
			if (dir1.dot(vec_mc1) > 0) dir1 *= -1;
			Vec dir = (dir0 * vec1.norm() + dir1 * vec0.norm()).normalized();// 以距离作为权值, 离得越远, 权重越小
			Vex new_point = curv + dir * mesh.radius[curvh] * radius_scale;// 默认向外扩展

			new_domain_points.push_back(new_point);
		}
		expanded_domains.push_back(new_domain_points);
		/*double perimeter_ori = 0, perimeter_exp = 0;
		for (int i = 0; i < domain.size(); ++i) {
			int j = (i + 1) % domain.size();
			perimeter_ori += (mesh.vertices(domain[i]) - mesh.vertices(domain[j])).norm();
			perimeter_exp += (new_domain_points[i] - new_domain_points[j]).norm();
		}
		std::cout << "Bounadry expand success. Perimeter: " << perimeter_ori << " expand to " << perimeter_exp 
			<< ". Check size flag = " << int(domain.size() == new_domain_points.size()) << std::endl;*/
	}

}

void DomainProcessor::remeshing_3d(SkeletalMesh& mesh,
	vector<vector<VH>>& domains_vhs, 
	vector<vector<Vex>>& expanded_domains, 
	vector<MeshKernel::SurfaceMesh>& trimeshs) {

	int dcnt = domains_vhs.size();// domain count

	for (int d_idx = 0; d_idx < dcnt; ++d_idx) {

		MeshKernel::SurfaceMesh trimesh;
		std::unordered_map<VH, VH> sklvh_to_trivh;
		std::unordered_map<VH, VH> trivh_to_sklvh;
		auto& bd_sklvhs = domains_vhs[d_idx];// 边界顶点
		auto& sklvs = expanded_domains[d_idx];// 扩张后的边界顶点位置
		int vcnt = bd_sklvhs.size();
		//std::cout << "Project to 2D: vs.size = " << sklvs.size() << ", vhs.size = " << bd_sklvhs.size() << std::endl;
		for (int i = 0; i < vcnt; ++i) {// 所有边界顶点都存入扩张后的顶点位置
			VH trivh = trimesh.AddVertex(sklvs[i]);
			sklvh_to_trivh[bd_sklvhs[i]] = trivh;
			trivh_to_sklvh[trivh] = bd_sklvhs[i];
		}

		// 将内部顶点也放入 sklvh
		std::unordered_set<FH> sklfhs;// 这块曲面所有的 fh
		std::queue<FH> que_fh;

		for (auto& sklvh : bd_sklvhs) {
			for (auto& fh : mesh.NeighborFh(sklvh)) {
				if (sklfhs.count(fh)) continue;
				sklfhs.insert(fh);
				que_fh.emplace(fh);
			}
		}

		while (!que_fh.empty()) {
			FH fh = que_fh.front();
			que_fh.pop();
			for (auto& adjfh : mesh.NeighborFh(fh)) {
				if (sklfhs.count(adjfh)) continue;
				sklfhs.insert(adjfh);
				que_fh.emplace(adjfh);
			}

			// 保存所有顶点
			const auto& f_vhs = mesh.faces(fh).getVertexHandle();
			for (auto& sklvh : f_vhs) {
				if (sklvh_to_trivh.count(sklvh)) continue;
				VH trivh = trimesh.AddVertex(mesh.vertices(sklvh));
				sklvh_to_trivh[sklvh] = trivh;
				trivh_to_sklvh[trivh] = sklvh;
			}

		}

		unordered_map<VH, unordered_set<VH>> pre_vh;// 半边数据结构中的前一个顶点，防止出现复杂边
		unordered_map<FH, bool> added;
		queue<FH> add_cand_fhs;
		add_cand_fhs.emplace(*(sklfhs.begin()));
		added[*(sklfhs.begin())] = true;

		//std::cout << "Project to 2D 167: vs.size = " << sklvs.size() << ", vhs.size = " << bd_sklvhs.size() << std::endl;

		while (!add_cand_fhs.empty()) {
			auto sklfh = add_cand_fhs.front();
			add_cand_fhs.pop();
			auto& sklf = mesh.faces(sklfh);
			auto fvhs = sklf.getVertexHandle();
			bool add_flag = true;
			for (auto& fvh : fvhs) {
				if (!sklvh_to_trivh.count(fvh)) {
					add_flag = false;
					std::cerr << "DomainProcessor.cpp Error: get a invalid face.\n";
					break;
				}
				fvh = sklvh_to_trivh[fvh];
			}
			if (add_flag) {
				// 调整面的方向
				for (int i = 0; i < fvhs.size(); ++i) {
					VH prevh = fvhs[i == 0 ? (fvhs.size() - 1) : i - 1];
					VH curvh = fvhs[i];
					if (pre_vh[curvh].count(prevh)) {
						std::reverse(fvhs.begin(), fvhs.end());
						break;
					}
				}
				for (int i = 0; i < fvhs.size(); ++i) {
					VH prevh = fvhs[i == 0 ? (fvhs.size() - 1) : i - 1];
					VH curvh = fvhs[i];
					if (pre_vh[curvh].count(prevh)) {
						std::cerr << "\nDomainProcessor Error: find complicated edge!\n\n";
						break;
					}
					pre_vh[curvh].insert(prevh);
				}
				trimesh.AddFace(fvhs);
			}
			for (auto& adjfh : mesh.NeighborFh(sklfh)) {
				if (added.count(adjfh)) continue;
				add_cand_fhs.emplace(adjfh);
				added[adjfh] = true;
			}
		}
		
		MeshKernel::IO io_origin;
		io_origin.WriteOffFile(trimesh, "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(d_idx) + "_origin.off");
	
		trimeshs.emplace_back(trimesh);

		Remeshing_2 app_remeshing(trimesh);
		app_remeshing.execute_pmp(0.5);// 加密 3D 网格, 用于生成高质量的四边形网格

		MeshKernel::IO io_remeshing;
		io_remeshing.WriteOffFile(trimesh, "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(d_idx) + "_remeshing.off");
		
	}

}


bool DomainProcessor::is_manifold(SkeletalMesh& mesh) {
	// 检查有无邻接面个数大于3的边
	for (auto& ep : mesh.alledges()) {
		if (mesh.NeighborFh(ep.first).size() > 2) {
			return false;
		}
	}

	return true;
}

void DomainProcessor::expand_boundary(SkeletalMesh& mesh, vector<vector<Vex>>& expanded_bdy) {
	mesh.genAllFacesNormal();
	unordered_map<VH, Vex> vh2pos;
	unordered_map<VH, int> vh2cnt;
	for (auto& ep : mesh.alledges()) {
		if (mesh.NeighborFh(ep.first).size() == 1) {
			auto fh = *(mesh.NeighborFh(ep.first).begin());
			Vex normal = mesh.getFaceNormal(fh);
			Vex midpoint = mesh.getEdgeMidpoint(ep.first);
			Vex centroid = mesh.getFaceCenter(fh);
			Vex dir_in = (centroid - midpoint).normalized();
			Vex dir_edge = (mesh.vertices(ep.second.vh1()) - mesh.vertices(ep.second.vh2())).normalized();
			Vex dir_move = normal.cross(dir_edge);
			double len = (mesh.radius[ep.second.vh1()] + mesh.radius[ep.second.vh2()]) / 2.0;
			if (dir_move.dot(dir_in) > 0.f) {
				dir_move *= -1;
			}
			vector<VH> vhs = { ep.second.vh1(), ep.second.vh2() };
			for (auto& vh : vhs) {
				Vex v_new = mesh.vertices(vh) + dir_move * len;
				if (vh2cnt.count(vh)) {
					++vh2cnt[vh];
					vh2pos[vh] += v_new;
				} else {
					vh2cnt[vh] = 1;
					vh2pos[vh] = v_new;
				}
			}
		}
	}

	// 是否扩张边界?
	for (auto& vp : vh2pos) {
		Vex& v = mesh.vertices(vp.first);
		v = vp.second / vh2cnt[vp.first];
	}

}

void DomainProcessor::decomposition(SkeletalMesh& mesh) {
	// 将面片区域分解, 然后输出
	vector<MeshKernel::SurfaceMesh> trimeshes;
	unordered_set<FH> visited;
	for (auto& fp : mesh.allfaces()) {
		FH fh = fp.first;
		if (visited.count(fh)) continue;
		unordered_set<FH> domain_fhs;
		queue<FH> que;
		que.push(fh);
		domain_fhs.insert(fh);
		while (!que.empty()) {
			FH curr_fh = que.front();
			que.pop();
			auto face = mesh.faces(curr_fh);
			auto face_ehs = face.getEdgeHandle();
			for (auto& f_eh : face_ehs) {
				auto e_fhs = mesh.NeighborFh(f_eh);
				if (e_fhs.size() > 2) continue;// 非流形的边
				for (auto& e_fh : e_fhs) {
					if (domain_fhs.count(e_fh)) continue;
					que.push(e_fh);
					domain_fhs.insert(e_fh);
				}
			}
		}

		// 输出三角形网格
		MeshKernel::SurfaceMesh trimesh;
		unordered_map<VH, VH> sklvh2trivh;
		unordered_map<VH, unordered_set<VH>> pre_vh;// 半边数据结构中的前一个顶点，防止出现复杂边
		unordered_set<FH> added;
		que = queue<FH>();
		que.push(*domain_fhs.begin());
		added.insert(*domain_fhs.begin());
		while (!que.empty()) {
			auto fh = que.front();
			que.pop();
			auto f_vhs = mesh.faces(fh).getVertexHandle();
			for (auto& f_vh : f_vhs) {
				if (!sklvh2trivh.count(f_vh)) {
					sklvh2trivh[f_vh] = trimesh.AddVertex(mesh.vertices(f_vh));
				}
				f_vh = sklvh2trivh[f_vh];
			}
			// 调整面的方向
			for (int i = 0; i < f_vhs.size(); ++i) {
				VH prevh = f_vhs[i == 0 ? (f_vhs.size() - 1) : i - 1];
				VH curvh = f_vhs[i];
				if (pre_vh[curvh].count(prevh)) {
					std::reverse(f_vhs.begin(), f_vhs.end());
					break;
				}
			}
			for (int i = 0; i < f_vhs.size(); ++i) {
				VH prevh = f_vhs[i == 0 ? (f_vhs.size() - 1) : i - 1];
				VH curvh = f_vhs[i];
				if (pre_vh[curvh].count(prevh)) {
					std::cerr << "\nDomainProcessor Error: find complicated edge!\n\n";
					break;
				}
				pre_vh[curvh].insert(prevh);
			}
			trimesh.AddFace(f_vhs);
			for (auto& adjfh : mesh.NeighborFh(fh)) {
				if (domain_fhs.count(adjfh)) {
					if (added.count(adjfh)) continue;
					que.push(adjfh);
					added.insert(adjfh);
				}
			}
		}

		visited.insert(domain_fhs.begin(), domain_fhs.end());
		trimeshes.emplace_back(trimesh);

	}

	sort(trimeshes.begin(), trimeshes.end(), [&](MeshKernel::SurfaceMesh& trimesh1, MeshKernel::SurfaceMesh& trimesh2) {
		return trimesh1.vsize() > trimesh2.vsize();
		});

	for (uint i = 0; i < trimeshes.size(); ++i) {

		auto& trimesh = trimeshes[i];

		MeshKernel::IO io_origin;
		io_origin.WriteOffFile(trimesh, "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(i) + "_origin.off");

		Remeshing_2 app_remeshing(trimesh);
		app_remeshing.execute_pmp(0.5);// 加密 3D 网格, 用于生成高质量的四边形网格

		MeshKernel::IO io_remeshing;
		io_remeshing.WriteOffFile(trimesh, "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(i) + "_remeshing.off");

	}

	std::cout << "Find " << trimeshes.size() << " meshes.\n";

}