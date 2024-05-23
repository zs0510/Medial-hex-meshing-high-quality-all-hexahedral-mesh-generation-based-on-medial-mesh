#include "Simplification_Powered_QEM.h"

void Simplification_Powered_QEM::initFaceNormal() {// 计算每个平面的法向量
	auto& faces = mesh.allfaces();
	auto fcnt = faces.size();
	Normals.resize(fcnt, Eigen::Vector3d::Zero());
	for (auto f : faces) {
		auto vex = f.second.getVertexHandle();
		assert(vex.size() > 2);
		Eigen::Vector3d pos1(mesh.vertices(vex[0]).x(), mesh.vertices(vex[0]).y(), mesh.vertices(vex[0]).z());
		Eigen::Vector3d pos2(mesh.vertices(vex[1]).x(), mesh.vertices(vex[1]).y(), mesh.vertices(vex[1]).z());
		Eigen::Vector3d pos3(mesh.vertices(vex[2]).x(), mesh.vertices(vex[2]).y(), mesh.vertices(vex[2]).z());
		Normals[f.first] = (pos3 - pos2).cross(pos1 - pos2);
		Normals[f.first].normalize();
	}
	printf("init face normal success\n");
}

void Simplification_Powered_QEM::initVertexErrorMat() {
	auto& vertices = mesh.allvertices();
	auto vcnt = vertices.size();
	ErrorMats.resize(vcnt, Eigen::Matrix4d::Zero());
	for (auto v : vertices) {
		for (auto f : mesh.NeighborFh(v.first)) {
			double d = Normals[f][0] * v.second.x() + Normals[f][1] * v.second.y() + Normals[f][2] * v.second.z();
			Eigen::Vector4d p(Normals[f][0], Normals[f][1], Normals[f][2], -d);
			ErrorMats[v.first] += (p * p.transpose());
		}
		for (auto eh : mesh.NeighborEh(v.first)) {
			if (mesh.isOnBoundary(eh)) {
				auto uh = (mesh.edges(eh)).vh1() + (mesh.edges(eh)).vh2() - v.first;
				auto& u = mesh.vertices(MeshKernel::iGameVertexHandle(uh));
				Eigen::Vector3d vec_uv = Eigen::Vector3d(v.second.x() - u.x(), v.second.y() - u.y(), v.second.z() - u.z());
				auto fhs = mesh.NeighborFh(eh);
				if (fhs.empty()) break;
				MeshKernel::iGameFaceHandle fh;
				for (auto& tmp : fhs) {
					fh = tmp;
					break;
				}
				Eigen::Vector3d N = Normals[fh].cross(vec_uv);
				N.normalize();
				double d = N[0] * v.second.x() + N[1] * v.second.y() + N[2] * v.second.z();
				Eigen::Vector4d p(N[0], N[1], N[2], -d);
				ErrorMats[v.first] += (p * p.transpose());
			}
		}
	}
	
	printf("init vertex error matrix success\n");
}

void Simplification_Powered_QEM::initContractCost(MeshKernel::iGameEdge& e) {
	auto& eh = mesh.edgehandle(e);
	//if (mesh.isOnBoundary(eh)) return; // 边界边是否允许坍缩
	auto vh1 = e.vh1();
	auto vh2 = e.vh2();

	// intersection detection
	auto adj_vhs1 = mesh.NeighborFh(vh1);
	auto adj_vhs2 = mesh.NeighborFh(vh2);
	int adj_cnt = 0;
	for (auto& adj_vh1 : adj_vhs1) {
		if (find(adj_vhs2.begin(), adj_vhs2.end(), adj_vh1) != adj_vhs2.end()) {
			adj_cnt++;
			if (adj_cnt > 2) return;// dangerous!
		}
	}

	double cost = Double_MAX;
	Eigen::Vector4d X;// 最佳新顶点位置
	Eigen::Matrix4d Q = ErrorMats[vh1] + ErrorMats[vh2];
	Eigen::Matrix4d A;
	A << Q(0, 0), Q(0, 1), Q(0, 2), Q(0, 3),
		Q(0, 1), Q(1, 1), Q(1, 2), Q(1, 3),
		Q(0, 2), Q(1, 2), Q(2, 2), Q(2, 3),
		0.f, 0.f, 0.f, 1.f;

	if (A.determinant()) {
		// A 可逆，可通过解方程获取最佳位置
		Eigen::Vector4d b(0.f, 0.f, 0.f, 1.f);
		X = A.colPivHouseholderQr().solve(b);
		cost = X.transpose() * Q * X;
	} else {
		// A 不可逆，在 线段 v1-v2 上寻找局部最佳位置
		auto v1 = mesh.vertices(vh1);
		auto v2 = mesh.vertices(vh2);
		Eigen::Vector3d dir(v2.x() - v1.x(), v2.y() - v1.y(), v2.z() - v1.z());
		Eigen::Vector3d beg(v1.x(), v1.y(), v1.z());
		for (double i = 0.f; i <= 1.f; i += 0.05f) {
			Eigen::Vector3d cur = beg + dir * i;
			Eigen::Vector4d temp_pos(cur[0], cur[1], cur[2], 1);
			double temp_cost = temp_pos.transpose() * Q * temp_pos;
			if (temp_cost < cost) {
				cost = temp_cost;
				X = temp_pos;
			}
		}
	}

	// calculate triangle's aspect ratio
	double prio_ar, aspect_ratio_ori = 999.f, aspect_ratio_new = aspect_ratio_ori;

	std::vector<std::vector<MeshKernel::iGameVertexHandle>> vhs_vh1;
	for (auto fh : mesh.NeighborFh(vh1)) {
		aspect_ratio_ori = std::min(aspect_ratio_ori, getAspectRatio(fh));
		auto vex = (mesh.faces(fh)).getVertexHandle();
		if (find(vex.begin(), vex.end(), vh2) == vex.end()) {// 不含 vh2
			vhs_vh1.push_back(vex);
			vector<Eigen::Vector3d> pos(3);
			for (int i = 0; i < 2; ++i) {
				auto& v = mesh.vertices(vex[i]);
				pos[i] = Eigen::Vector3d(v.x(), v.y(), v.z());
			}
			Eigen::Vector3d N_ori = getNormal(pos);
			for (int i = 0; i < 2; ++i) {
				if (vex[i] == vh1) {
					pos[i] = Eigen::Vector3d(X[0], X[1], X[2]);
				}
			}
			Eigen::Vector3d N_new = getNormal(pos);
			double cos_alpha = N_ori.dot(N_new);
			if (cos_alpha < 0.f) return;
		}
	}
	for (auto fh : mesh.NeighborFh(vh2)) {
		aspect_ratio_ori = std::min(aspect_ratio_ori, getAspectRatio(fh));
		auto vex = (mesh.faces(fh)).getVertexHandle();
		if (find(vex.begin(), vex.end(), vh1) == vex.end()) {// 不含 vh1
			vector<Eigen::Vector3d> pos(3);
			for (int i = 0; i < 2; ++i) {
				auto& v = mesh.vertices(vex[i]);
				pos[i] = Eigen::Vector3d(v.x(), v.y(), v.z());
			}
			Eigen::Vector3d N_ori = getNormal(pos);
			for (int i = 0; i < 2; ++i) {
				if (vex[i] == vh2) {
					pos[i] = Eigen::Vector3d(X[0], X[1], X[2]);
				}
			}
			Eigen::Vector3d N_new = getNormal(pos);
			double cos_alpha = N_ori.dot(N_new);
			if (cos_alpha < 0.f) return;
		}
	}

	std::vector<std::vector<MeshKernel::iGameVertexHandle>> new_triangles;
	for (auto& fh : mesh.NeighborFh(vh2)) {
		auto vex = (mesh.faces(fh)).getVertexHandle();
		assert(vex.size() > 2);
		if (std::find(vex.begin(), vex.end(), vh1) != vex.end()) continue;
		for (int i = 0; i < vex.size(); ++i) {
			if (vex[i] == vh2) {
				vex[i] = vh1;
				new_triangles.push_back(vex);
				break;
			}
		}
	}
	
	new_triangles.insert(new_triangles.end(), vhs_vh1.begin(), vhs_vh1.end());
	vector<vector<Eigen::Vector3d>> positions(new_triangles.size(), vector<Eigen::Vector3d>(3));
	for (int i = 0; i < new_triangles.size(); ++i) {
		auto& vhs = new_triangles[i];
		for (int j = 0; j < 3; ++j) {
			auto& v = mesh.vertices(vhs[j]);
			if (vhs[j] != vh1) {
				positions[i][j] = Eigen::Vector3d(v.x(), v.y(), v.z());
			} else {
				positions[i][j] = Eigen::Vector3d(X[0], X[1], X[2]);
			}
		}
		aspect_ratio_new = min(aspect_ratio_new, getAspectRatio(positions[i]));
	}

	if (aspect_ratio_new > aspect_ratio_ori) prio_ar = 1 - aspect_ratio_new;
	else if (aspect_ratio_new > ASPECT_RATIO_THRESHOLD) prio_ar = 2 - aspect_ratio_new;
	else return;// 不合法

	double prio_nrd, min_nrd = 99999.f;

	for (auto& tri : positions) {
		min_nrd = min(min_nrd, getNormalizedRoundness(tri));
	}
	prio_nrd = 1 - min_nrd;
	//printf("prio_ar: %.4f, prio_nrd: %.4f\n", prio_ar, prio_nrd);

	cost /= (prio_ar * prio_nrd);

	powered_contract* tmp = new powered_contract(cost, X[0], X[1], X[2], vh1, vh2);
	adjContracts[vh1].push_back(tmp);
	adjContracts[vh2].push_back(tmp);
	Contracts.push(tmp);
}

void Simplification_Powered_QEM::initContractCosts() {
	auto& edges = mesh.alledges();
	for (auto e : edges) {
		initContractCost(e.second);
	}
	printf("init contract cost success\n");
}

Eigen::Vector3d Simplification_Powered_QEM::getNormal(vector<Eigen::Vector3d>& vec) {
	Eigen::Vector3d vec01 = vec[1] - vec[0];
	Eigen::Vector3d vec02 = vec[2] - vec[0];
	Eigen::Vector3d N = vec01.cross(vec02);
	N.normalize();
	return N;
}

double Simplification_Powered_QEM::getAspectRatio(vector<Eigen::Vector3d>& vex) {
	Eigen::Vector3d vec1 = vex[1] - vex[0];
	Eigen::Vector3d vec2 = vex[2] - vex[0];
	Eigen::Vector3d vec3 = vex[2] - vex[1];
	double len1 = vec1.norm(), len2 = vec2.norm(), len3 = vec3.norm();
	double double_area = (vec1.cross(vec2)).norm();
	double max_len = max(len1, max(len2, len3));
	double aspect_ratio = double_area / (max_len * max_len);
	return aspect_ratio;
}

double Simplification_Powered_QEM::getNormalizedRoundness(vector<Eigen::Vector3d>& vex) {
	Eigen::Vector3d vec1 = vex[1] - vex[0];
	Eigen::Vector3d vec2 = vex[2] - vex[0];
	Eigen::Vector3d vec3 = vex[2] - vex[1];
	double area = (vec1.cross(vec2)).norm() * 0.5f;
	double nrd = 1.1547f;
	double num1 = vec1.norm(), num2 = vec2.norm(), num3 = vec3.norm();
	double numerator = area * min(num1, min(num2, num3));
	double denominator = num1 * num2 * num3;
	nrd *= (numerator / denominator);
	return nrd;
}

double Simplification_Powered_QEM::getAspectRatio(MeshKernel::iGameFaceHandle& fh) {
	auto vhs = mesh.faces(fh).getVertexHandle();
	auto& v1 = mesh.vertices(vhs[0]);
	auto& v2 = mesh.vertices(vhs[1]);
	auto& v3 = mesh.vertices(vhs[2]);
	Eigen::Vector3d vec1 = Eigen::Vector3d(v2.x() - v1.x(), v2.y() - v1.y(), v2.z() - v1.z());
	Eigen::Vector3d vec2 = Eigen::Vector3d(v3.x() - v1.x(), v3.y() - v1.y(), v3.z() - v1.z());
	Eigen::Vector3d vec3 = Eigen::Vector3d(v2.x() - v3.x(), v2.y() - v3.y(), v2.z() - v3.z());
	double len1 = vec1.norm(), len2 = vec2.norm(), len3 = vec3.norm();
	double double_area = (vec1.cross(vec2)).norm();
	double max_len = max(len1, max(len2, len3));
	double aspect_ratio = double_area / (max_len * max_len);
	return aspect_ratio;
}

void Simplification_Powered_QEM::Execute() {// remainder is the num of remainder vertices
	time_t start = clock();
	int input_vcnt = mesh.allvertices().size();
	adjContracts.resize(input_vcnt);
	initFaceNormal();
	initVertexErrorMat();
	initContractCosts();
	while (it_num--) {

		powered_contract* min_contract;
		MeshKernel::iGameVertexHandle vh_keep(-1), vh_remove(-1);
		while (vh_keep == -1 || vh_remove == -1) {
			if (Contracts.empty()) break;
			min_contract = Contracts.top();
			Contracts.pop();
			vh_keep = min_contract->vh1;// 注意！
			vh_remove = min_contract->vh2;
		}
		if (vh_keep == -1 || vh_remove == -1) break;
		//printf("it%d: selected contract cost: %.6f\n", it_num, min_contract->cost);
		if (mesh.NeighborVh(vh_keep).size() < mesh.NeighborVh(vh_remove).size()) {
			std::swap(vh_keep, vh_remove);
		}
		auto& vex_keep = mesh.vertices(vh_keep);

		// 更新 vertex_keep 信息
		vex_keep.setPosition(min_contract->x, min_contract->y, min_contract->z);
		ErrorMats[vh_keep] += ErrorMats[vh_remove];

		std::vector<std::vector<MeshKernel::iGameVertexHandle>> new_triangles;
		for (auto& fh : mesh.NeighborFh(vh_remove)) {
			auto vex = (mesh.faces(fh)).getVertexHandle();
			assert(vex.size() > 2);
			if (std::find(vex.begin(), vex.end(), vh_keep) != vex.end()) continue;
			for (int i = 0; i < vex.size(); ++i) {
				if (vex[i] == vh_remove) {
					vex[i] = vh_keep;
					new_triangles.push_back(vex);
					break;
				}
			}
		}

		// 更新失效收缩
		for (auto* ctr : adjContracts[vh_remove]) {
			ctr->vh1 = ctr->vh2 = MeshKernel::iGameVertexHandle(-1);
		}
		adjContracts[vh_remove].clear();
		mesh.DeleteVertex(vh_remove);

		// 添加三角形
		for (auto& tri : new_triangles) {
			mesh.AddFace(tri);
		}

		// 更新失效收缩
		for (auto* ctr : adjContracts[vh_keep]) {
			ctr->vh1 = ctr->vh2 = MeshKernel::iGameVertexHandle(-1);
		}
		adjContracts[vh_keep].clear();

		// 添加新收缩
		for (auto& e : mesh.NeighborEh(vh_keep)) {
			initContractCost(mesh.edges(e));
		}
	}
	time_t end = clock();
	int output_vcnt = mesh.allvertices().size();
	printf("simplification success vcnt: %d -> %d\n", input_vcnt, output_vcnt);
	mesh.updateAllHandles();
	printf("running time: %dms\n", (int)(end - start));
}