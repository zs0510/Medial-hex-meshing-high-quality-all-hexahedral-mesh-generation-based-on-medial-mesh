#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

#ifndef EPSILON
#define EPSILON 1E-6F
#endif // !EPSILON

#ifndef D_PI
#define D_PI 6.2831852f
#endif

// Author: Zhang Sheng
// Time: 2021.8.27

class Parameterization {

public:
	Parameterization(MeshKernel::SurfaceMesh& _mesh, size_t boundary_mode = 2, size_t vertex_index_offset = 64) :
		mesh(_mesh), BoundaryMode(boundary_mode % 3), IndexOffset(vertex_index_offset) {
	};
	void Execute();

private:
	void initCotWeight();
	void findAllBoundary();
	void setIndexOffset();
	void setOriginBoundary();
	void setSquareBoundary();
	void setCircleBoundary();
	size_t BoundaryMode;// 0: origin sharp boundary; 1: squard; 2: circle
	size_t IndexOffset;// the first vertex_index
	MeshKernel::SurfaceMesh& mesh;
	std::vector<MeshKernel::iGameVertexHandle> vex_boundary;
	std::unordered_map<MeshKernel::iGameEdgeHandle, float> weights;
	std::vector<Eigen::Vector3f> positions;// 保存所有顶点的位置
	std::unordered_map<MeshKernel::iGameVertexHandle, size_t> indices;// 保存所有顶点的索引
	Eigen::VectorXf B_U, B_V;

};

void Parameterization::Execute() {
	size_t vcnt = mesh.allvertices().size(), idx = 0;
	B_U.resize(vcnt);
	B_V.resize(vcnt);
	positions.resize(vcnt, Eigen::Vector3f::Zero());
	for (auto v : mesh.allvertices()) {
		assert(idx < vcnt);
		indices[v.first] = idx;
		positions[idx++] = Eigen::Vector3f(v.second.x(), v.second.y(), v.second.z());
	}
	initCotWeight();
	findAllBoundary();
	switch (BoundaryMode) {
	case 0:
		printf("the func is needed to complete\n");
		return;
		break;
	case 1:
		setSquareBoundary();
		break;
	case 2:
		setCircleBoundary();
		break;
	}
	std::vector<Eigen::Triplet<float>> triplets;
	for (auto v : mesh.allvertices()) {// 构建稀疏矩阵
		if (mesh.isOnBoundary(v.first)) {
			triplets.emplace_back(indices[v.first], indices[v.first], 1);
		} else {
			float sum_weight = 0.f;
			for (auto eh : mesh.NeighborEh(v.first)) {
				auto u = (mesh.edges(eh).vh1() == v.first) ? mesh.edges(eh).vh2() : mesh.edges(eh).vh1();
				float weight = weights[eh];
				triplets.emplace_back(indices[v.first], indices[u], -weight);
				sum_weight += weight;
			}
			triplets.emplace_back(indices[v.first], indices[v.first], sum_weight);
			B_U[indices[v.first]] = 0.f;
			B_V[indices[v.first]] = 0.f;
		}
	}

	// 初始化 A 矩阵
	Eigen::SparseMatrix<float> A(vcnt, vcnt);
	A.setFromTriplets(triplets.begin(), triplets.end());

	Eigen::SparseMatrix<float> AT = A.transpose();
	Eigen::SparseMatrix<float> ATA = AT * A;

	// 方程组左边
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<float>> MatricesCholesky(ATA);

	// 方程组右边
	Eigen::VectorXf AT_BU = AT * B_U;
	Eigen::VectorXf AT_BV = AT * B_V;

	// 解
	Eigen::VectorXf U = MatricesCholesky.solve(AT_BU);
	Eigen::VectorXf V = MatricesCholesky.solve(AT_BV);

	for (auto v : mesh.allvertices()) {
		MeshKernel::iGameVertex& vex = mesh.vertices(v.first);
		vex.setPosition(U[indices[v.first]], V[indices[v.first]], 0.f);
	}
	printf("parameterization success\n");
}

void Parameterization::initCotWeight() {
	for (auto edge : mesh.alledges()) {
		auto vh = edge.second.vh1();
		auto uh = edge.second.vh2();
		auto v = mesh.vertices(vh);
		auto u = mesh.vertices(uh);
		float cot_weight = 0.f;
		for (auto face : mesh.NeighborFh(edge.first)) {
			int idx = 0;
			auto wh = mesh.faces(face).vh(idx++);
			while (wh == vh || wh == uh) {
				assert(idx < 3);
				wh = mesh.faces(face).vh(idx++);
			}
			auto w = mesh.vertices(wh);
			Eigen::Vector3f v_pos(v.x(), v.y(), v.z());
			Eigen::Vector3f u_pos(u.x(), u.y(), u.z());
			Eigen::Vector3f w_pos(w.x(), w.y(), w.z());
			Eigen::Vector3f vec1 = (v_pos - w_pos).normalized();
			Eigen::Vector3f vec2 = (u_pos - w_pos).normalized();
			float cos_theta = vec1.dot(vec2);
			float cot_theta = cos_theta / std::sqrt(1 - cos_theta * cos_theta);
			cot_weight += std::max(EPSILON, cot_theta);
		}
		weights[edge.first] = cot_weight * 0.5f;
	}
	//printf("init cot weight success\n");
}

void Parameterization::findAllBoundary() {
	MeshKernel::iGameVertexHandle first_vex, pre_vex, cur_vex;
	for (auto v : mesh.allvertices()) {
		if (mesh.isOnBoundary(v.first)) {
			pre_vex = first_vex = v.first;
			MeshKernel::iGameVertexHandle a, b;
			bool tag = true;
			for (auto adjvh : mesh.NeighborVh(v.first)) {// 找到其两个相邻边界顶点
				if (mesh.isOnBoundary(adjvh)) {
					if (tag) {
						a = adjvh;
						tag = false;
					} else {
						b = adjvh;
						break;
					}
				}
			}
			auto va = mesh.vertices(a);
			auto vb = mesh.vertices(b);
			Eigen::Vector3f v1(v.second.x(), v.second.y(), v.second.z());
			Eigen::Vector3f v2(va.x(), va.y(), va.z());
			Eigen::Vector3f v3(vb.x(), vb.y(), vb.z());
			Eigen::Vector3f tmp = (v2 - v1).cross(v3 - v1);
			if (tmp[2] > 0) {
				cur_vex = a;
			} else {
				cur_vex = b;
			}
			break;
		}
	}
	vex_boundary.push_back(pre_vex);
	//printf("boundary: %d\n", int(pre_vex));
	while (cur_vex != first_vex) {
		vex_boundary.push_back(cur_vex);
		//printf("boundary: %d\n", int(cur_vex));
		for (auto vh : mesh.NeighborVh(cur_vex)) {
			if (vh != pre_vex && mesh.isOnBoundary(vh)) {
				pre_vex = cur_vex;
				cur_vex = vh;
				break;
			}
		}
	}
	setIndexOffset();
	//printf("find boundary success: %d\n", int(vex_boundary.size()));
}

void Parameterization::setIndexOffset() {
	auto sz = vex_boundary.size();
	IndexOffset %= sz;
	if (IndexOffset == 0) return;
	std::vector<MeshKernel::iGameVertexHandle> tmp(sz);
	for (int i = 0; i < sz; ++i) {
		tmp[i] = vex_boundary[(IndexOffset + i) % sz];
	}
	vex_boundary = tmp;
}

void Parameterization::setOriginBoundary() {
	
}

void Parameterization::setSquareBoundary() {
	size_t boundary_sz = vex_boundary.size();
	std::vector<float> chord_length(boundary_sz + 1, 0);
	for (int i = 1; i <= boundary_sz; ++i) {
		Eigen::Vector3f pre = positions[indices[vex_boundary[i - 1]]], cur;
		if (i == boundary_sz) cur = positions[indices[vex_boundary[0]]];
		else cur = positions[indices[vex_boundary[i]]];
		chord_length[i] = chord_length[i - 1] + (cur - pre).norm();
	}
	//printf("cul chord_length success\n");
	float qtr1 = chord_length.back() / 4.0f;
	float qtr2 = qtr1 * 2.0f;
	float qtr3 = qtr1 * 3.0f;
	std::vector<int> idxs(5, -1);// index of end point
	for (int i = 1; i < chord_length.size(); ++i) {
		if (idxs[1] == -1 && chord_length[i] > qtr1) {
			float diff1 = abs(chord_length[i] - qtr1);
			float diff2 = (i - 1 >= 0) ? abs(chord_length[i - 1] - qtr1) : std::numeric_limits<float>::max();
			idxs[1] = (diff1 < diff2) ? i : (i - 1);
		} else if (idxs[2] == -1 && chord_length[i] > qtr2) {
			float diff1 = abs(chord_length[i] - qtr2);
			float diff2 = (i - 1 >= 0) ? abs(chord_length[i - 1] - qtr2) : std::numeric_limits<float>::max();
			idxs[2] = (diff1 < diff2) ? i : (i - 1);
		} else if (idxs[3] == -1 && chord_length[i] > qtr3) {
			float diff1 = abs(chord_length[i] - qtr3);
			float diff2 = (i - 1 >= 0) ? abs(chord_length[i - 1] - qtr3) : std::numeric_limits<float>::max();
			idxs[3] = (diff1 < diff2) ? i : (i - 1);
		}
	}
	idxs[0] = 0, idxs[4] = boundary_sz;
	//printf("init idxs success\n");
	for (int i = 0; i < 4; ++i) {
		// 映射 index 在 [beg, end) 上的点
		int beg = idxs[i], end = idxs[i + 1];
		//assert(end < chord_length.size());
		float length = chord_length[end];
		int index = indices[vex_boundary[beg]];
		//assert(index < mesh.allvertices().size());
		switch (i) {
		case 0:
			B_U[index] = 0.0f;
			B_V[index] = 0.0f;
			for (int j = beg + 1; j < end; ++j) {
				index = indices[vex_boundary[j]];
				//assert(index < mesh.allvertices().size());
				B_U[index] = (chord_length[j] - chord_length[beg]) / (chord_length[end] - chord_length[beg]);
				B_V[index] = 0.0f;
			}
			break;
		case 1:
			B_U[index] = 1.0f;
			B_V[index] = 0.0f;
			for (int j = beg + 1; j < end; ++j) {
				index = indices[vex_boundary[j]];
				//assert(index < mesh.allvertices().size());
				B_U[index] = 1.0f;
				B_V[index] = (chord_length[j] - chord_length[beg]) / (chord_length[end] - chord_length[beg]);
			}
			break;
		case 2:
			B_U[index] = 1.0f;
			B_V[index] = 1.0f;
			for (int j = beg + 1; j < end; ++j) {
				index = indices[vex_boundary[j]];
				//assert(index < mesh.allvertices().size());
				B_U[index] = 1.0f - (chord_length[j] - chord_length[beg]) / (chord_length[end] - chord_length[beg]);
				B_V[index] = 1.0f;
			}
			break;
		case 3:
			B_U[index] = 0.0f;
			B_V[index] = 1.0f;
			for (int j = beg + 1; j < end; ++j) {
				index = indices[vex_boundary[j]];
				//assert(index < mesh.allvertices().size());
				B_U[index] = 0.0f;
				B_V[index] = 1.0f - (chord_length[j] - chord_length[beg]) / (chord_length[end] - chord_length[beg]);
			}
			break;
		}
	}
	//printf("set boundary success\n");
}

void Parameterization::setCircleBoundary() {
	int sz = vex_boundary.size();
	std::vector<float> chord_length(sz + 1, 0);
	for (int i = 1; i < sz; ++i) {
		Eigen::Vector3f pre = positions[indices[vex_boundary[i - 1]]];
		Eigen::Vector3f cur = positions[indices[vex_boundary[i]]];
		chord_length[i] = chord_length[i - 1] + (cur - pre).norm();
	}
	chord_length[sz] = chord_length[sz - 1] + (positions[indices[vex_boundary[0]]] - positions[indices[vex_boundary[sz - 1]]]).norm();
	float R = 0.5f;
	for (int i = 0; i < sz; ++i) {
		int index = indices[vex_boundary[i]];
		float theta = chord_length[i] / chord_length.back() * D_PI;
		B_U[index] = R + R * cos(theta);
		B_V[index] = R + R * sin(theta);
	}
}