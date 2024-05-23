#pragma once
#include "Kernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

// Author: Zhang Sheng
// Time: 2021.8.27
// Reference: 
//	Paper: As-Rigid-As-Possible Shape Interpolation
//	Code: https://github.com/USTC-GCL-F/AMMesh/tree/main/src/hw8

class MeshInterpolation_ARAP {// support 2D Mesh only(z == 0)
public:
	MeshInterpolation_ARAP(MeshKernel::SurfaceMesh& source_mesh, MeshKernel::SurfaceMesh& target_mesh, size_t num_of_interpolation_mesh) :
		source(source_mesh), target(target_mesh), delta_t(1.f / (num_of_interpolation_mesh + 1)) {
	};
	void Execute();
private:
	MeshKernel::SurfaceMesh& source, target;
	float delta_t;
};

void MeshInterpolation_ARAP::Execute() {
	assert(source.VertexSize() == target.VertexSize());
	assert(source.FaceSize() == target.FaceSize());
	auto fcnt = source.FaceSize();
	auto vcnt = source.VertexSize();
	std::vector<float> area(fcnt), angle(fcnt);
	std::vector<Eigen::Vector3f> src_pos(vcnt), tag_pos(vcnt);
	std::vector<Eigen::Vector3f> y_sub(fcnt), x_sub(fcnt);
	std::vector<Eigen::Matrix2f> S(fcnt);
	
	for (size_t i = 0; i < vcnt; ++i) {
		auto v1 = source.vertices(MeshKernel::iGameVertexHandle(i));
		src_pos[i] = Eigen::Vector3f(v1.x(), v1.y(), v1.z());
		auto v2 = target.vertices(MeshKernel::iGameVertexHandle(i));
		tag_pos[i] = Eigen::Vector3f(v2.x(), v2.y(), v2.z());
	}

	for (size_t i = 0; i < fcnt; ++i) {
		auto face = source.faces(MeshKernel::iGameFaceHandle(i));
		auto vex = face.getVertexHandle();
		assert(size_t(vex[0]) < vcnt && size_t(vex[1]) < vcnt && size_t(vex[2]) < vcnt);
		area[i] = (src_pos[vex[1]] - src_pos[vex[0]]).cross(src_pos[vex[2]] - src_pos[vex[0]]).norm() * 0.5f;
		// (y_j - y_k, y_k - y_i, y_i - y_j)
		y_sub[i] = Eigen::Vector3f(src_pos[vex[1]].y() - src_pos[vex[2]].y(), src_pos[vex[2]].y() - src_pos[vex[0]].y(), src_pos[vex[0]].y() - src_pos[vex[1]].y());
		// (x_k - x_j, x_i - x_k, x_j - x_i)
		x_sub[i] = Eigen::Vector3f(src_pos[vex[2]].x() - src_pos[vex[1]].x(), src_pos[vex[0]].x() - src_pos[vex[2]].x(), src_pos[vex[1]].x() - src_pos[vex[0]].x());

		Eigen::Vector3f u(tag_pos[vex[0]].x(), tag_pos[vex[1]].x(), tag_pos[vex[2]].x());
		Eigen::Vector3f v(tag_pos[vex[0]].y(), tag_pos[vex[1]].y(), tag_pos[vex[2]].y());

		Eigen::Matrix2f J;// sur_pos * J = tar_pos
		J << y_sub[i].dot(u), x_sub[i].dot(u), y_sub[i].dot(v), x_sub[i].dot(v);
		J /= (area[i] * 2);

		Eigen::JacobiSVD<Eigen::Matrix2f> svd(J, Eigen::ComputeFullV | Eigen::ComputeFullU);

		Eigen::Matrix2f V = svd.matrixV();
		Eigen::Matrix2f U = svd.matrixU();
		Eigen::Matrix2f R = U * V.transpose();
		angle[i] = atan2(R(1, 0), R(1, 1));// atan2(y,x) = 以坐标原点为起点，指向(y,x)的射线在坐标平面上与x轴正方向之间的角度

		Eigen::Matrix2f sigma;
		sigma(0, 1) = sigma(1, 0) = 0;
		sigma(0, 0) = svd.singularValues()[0];
		sigma(1, 1) = svd.singularValues()[1];
		S[i] = V * sigma * V.transpose();
	}

	// 固定一个点
	float fix_x = src_pos[vcnt - 1].x();
	float fix_y = src_pos[vcnt - 1].y();

	std::vector<Eigen::Triplet<float>> triplets;

	for (size_t i = 0; i < fcnt; ++i) {
		auto face = source.faces(MeshKernel::iGameFaceHandle(i));
		auto vex = face.getVertexHandle();
		for (size_t j = 0; j < 3; ++j) {
			for (size_t k = 0; k < 3; ++k) {
				if (vex[j] == (vcnt - 1) || vex[k] == (vcnt - 1)) continue;
				triplets.emplace_back(2 * vex[j], 2 * vex[k], y_sub[i][j] * y_sub[i][k] / (2.f * area[i] * area[i]));
				triplets.emplace_back(2 * vex[j], 2 * vex[k], x_sub[i][j] * x_sub[i][k] / (2.f * area[i] * area[i]));
				triplets.emplace_back(2 * vex[j] + 1, 2 * vex[k] + 1, y_sub[i][j] * y_sub[i][k] / (2.f * area[i] * area[i]));
				triplets.emplace_back(2 * vex[j] + 1, 2 * vex[k] + 1, x_sub[i][j] * x_sub[i][k] / (2.f * area[i] * area[i]));
			}
		}
	}

	Eigen::SparseMatrix<float> A(2 * (vcnt - 1), 2 * (vcnt - 1));
	// setFromTriplets 会将 triplets 中相同位置的值累加
	A.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
	solver.analyzePattern(A);
	solver.factorize(A);

	
	float t = delta_t;
	Eigen::Matrix2f I = Eigen::Matrix2f::Identity();
	while (t < 1.f) {
		Eigen::VectorXf b(2 * (vcnt - 1));
		b.setZero();
		for (size_t i = 0; i < fcnt; ++i) {
			auto face = source.faces(MeshKernel::iGameFaceHandle(i));
			auto vex = face.getVertexHandle();
			float theta = angle[i] * t;
			Eigen::Matrix2f R;
			R << cos(theta), -sin(theta), sin(theta), cos(theta);
			Eigen::Matrix2f A = R * ((1.f - t) * I + t * S[i]);// theoretical mapping
			for (size_t j = 0; j < 3; ++j) {
				if (vex[j] == vcnt - 1) continue;
				for (size_t k = 0; k < 3; ++k) {
					if (vex[k] == vcnt - 1) {
						b[2 * vex[j]] += -y_sub[i][j] * y_sub[i][k] * fix_x / (2.f * area[i] * area[i]);
						b[2 * vex[j]] += -x_sub[i][j] * x_sub[i][k] * fix_x / (2.f * area[i] * area[i]);
						b[2 * vex[j] + 1] += -y_sub[i][j] * y_sub[i][k] * fix_y / (2.f * area[i] * area[i]);
						b[2 * vex[j] + 1] += -x_sub[i][j] * x_sub[i][k] * fix_y / (2.f * area[i] * area[i]);
					}
				}
				b[2 * vex[j]] += y_sub[i][j] * A(0, 0) / area[i];
				b[2 * vex[j]] += x_sub[i][j] * A(0, 1) / area[i];
				b[2 * vex[j] + 1] += y_sub[i][j] * A(1, 0) / area[i];
				b[2 * vex[j] + 1] += x_sub[i][j] * A(1, 1) / area[i];
			}
		}
		Eigen::VectorXf xy = solver.solve(b);
		for (int i = 0; i < vcnt - 1; ++i) {
			auto& v = source.vertices(MeshKernel::iGameVertexHandle(i));
			v.setPosition(xy[2 * i], xy[2 * i + 1], 0.f);
		}
		std::string filename = "interpolation " + (std::to_string(t)).substr(0, 4) + ".obj";
		std::cout << "interpolation " + (std::to_string(t)).substr(0, 4) + " success" << std::endl;
		MeshKernel::IO io;
		io.WriteFile(source, filename);
		t += delta_t;
	}
}
