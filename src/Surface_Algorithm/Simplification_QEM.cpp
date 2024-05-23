#include "Simplification_QEM.h"

void Simplification_QEM::initFaceNormal() {// 计算每个平面的法向量
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

void Simplification_QEM::initVertexErrorMat() {
	auto& vertices = mesh.allvertices();
	auto vcnt = vertices.size();
	ErrorMats.resize(vcnt, Eigen::Matrix4d::Zero());
	for (auto v : vertices) {
		for (auto f : mesh.NeighborFh(v.first)) {
			double d = Normals[f][0] * v.second.x() + Normals[f][1] * v.second.y() + Normals[f][2] * v.second.z();
			Eigen::Vector4d p(Normals[f][0], Normals[f][1], Normals[f][2], -d);
			ErrorMats[v.first] += (p * p.transpose());
		}
	}
	printf("init vertex error matrix success\n");
}

void Simplification_QEM::initContractCost(const EH& eh) {
	//auto& eh = mesh.edgehandle(e);
	auto& e = mesh.edges(eh);
	auto vh1 = e.vh1();
	auto vh2 = e.vh2();
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
	contract* tmp = new contract(cost, X[0], X[1], X[2], vh1, vh2);
	adjContracts[vh1].push_back(tmp);
	adjContracts[vh2].push_back(tmp);
	Contracts.push(tmp);
}

void Simplification_QEM::initContractCosts() {
	auto& edges = mesh.alledges();
	for (auto& e : edges) {
		initContractCost(e.first);
	}
	printf("init contract cost success\n");
}

void Simplification_QEM::Execute() {// remainder is the num of remainder vertices
	time_t start = clock();
	int input_vcnt = mesh.allvertices().size();
	adjContracts.resize(input_vcnt);
	initFaceNormal();
	initVertexErrorMat();
	initContractCosts();
	while (it_num--) {

		contract* min_contract;
		MeshKernel::iGameVertexHandle vh_keep(-1), vh_remove(-1);
		while (vh_keep == -1 || vh_remove == -1) {
			if (Contracts.empty()) break;// code - 1
			min_contract = Contracts.top();
			Contracts.pop();
			vh_keep = min_contract->vh1;
			vh_remove = min_contract->vh2;
		}
		if (vh_keep == -1 || vh_remove == -1) break;// code - 2
		//printf("it%d: selected contract cost: %.10f\n", it_num, min_contract->cost);
		if (mesh.NeighborVh(vh_keep).size() < mesh.NeighborVh(vh_remove).size()) {
			std::swap(vh_keep, vh_remove);
		}
		auto& vex_keep = mesh.vertices(vh_keep);

		// 更新 vertex_keep 信息
		vex_keep.setPosition(min_contract->x, min_contract->y, min_contract->z);
		ErrorMats[vh_keep] += ErrorMats[vh_remove];


		// 如何避免面翻转?
		// 1. 参数化，根据 uv 坐标逆时针来确定顺序
		// 2. 计算 三角形法向量 与 向量(顶点->质心) 的夹角
		// 3. 直接用 vex_keep 代替 vex_remove

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
		for (auto& eh : mesh.NeighborEh(vh_keep)) {
			initContractCost(eh);
		}
	}
	time_t end = clock();
	int output_vcnt = mesh.allvertices().size();
	printf("simplification success vcnt: %d -> %d\n", input_vcnt, output_vcnt);
	mesh.updateAllHandles();
	printf("running time: %dms\n", (int)(end - start));
}