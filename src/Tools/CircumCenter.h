#pragma once
#include <vector>
#include <Eigen/Core>

using namespace std;

namespace CircumCenter {

	double Helen(Eigen::Vector3d& A, Eigen::Vector3d& B, Eigen::Vector3d& C);
	bool solve(Eigen::Vector3d& A, Eigen::Vector3d& B, Eigen::Vector3d& C, Eigen::Vector3d& D, Eigen::Vector3d& O);
	double get_radius_and_center(Eigen::Vector3d& A, Eigen::Vector3d& B, Eigen::Vector3d& C, Eigen::Vector3d& D, Eigen::Vector3d& O);

	bool get_circum_center(vector<Eigen::Vector3d>& points, Eigen::Vector3d& center, double& radius) {

		Eigen::Vector3d O(0, 0, 0);
		if (points.size() != 4) return false;
		auto& A = points[0];
		auto& B = points[1];
		auto& C = points[2];
		auto& D = points[3];
		if (!solve(A, B, C, D, O)) {
			return false;// 四点共面
		}
		solve(A, B, D, C, O);
		solve(A, C, D, B, O);
		solve(B, C, D, A, O);
		center = O;
		radius = get_radius_and_center(A, B, C, D, O);
		return true;
	}

	

	double Helen(Eigen::Vector3d& A, Eigen::Vector3d& B, Eigen::Vector3d& C) {

		auto AB = B - A;
		double c = AB.norm();
		auto AC = C - A;
		double b = AC.norm();
		auto BC = C - B;
		double a = BC.norm();
		double p = (a + b + c) / 2;
		return sqrt(p * (p - a) * (p - b) * (p - c));
	}

	bool solve(Eigen::Vector3d& A, Eigen::Vector3d& B, Eigen::Vector3d& C, Eigen::Vector3d& D, Eigen::Vector3d& O) {

		auto AB = B - A;
		auto BC = C - B;
		auto n = AB.cross(BC);
		double ret = (D - A).dot(n);
		if (fabs(ret) < 1e-8f) return false;// 四点共面
		double S = Helen(A, B, C);
		auto tmp_d = D * S;
		O += tmp_d;
		return true;
	}

	double get_radius_and_center(Eigen::Vector3d& A, Eigen::Vector3d& B, Eigen::Vector3d& C, Eigen::Vector3d& D, Eigen::Vector3d& O) {

		auto AB = B - A;
		auto BC = C - B;
		auto n = AB.cross(BC);
		auto AD = D - A;
		if (AD.dot(n) < 0) n *= -1;
		double h = AD.dot(n) / n.norm();
		double V4 = h * Helen(A, B, C);
		double s = Helen(A, B, C) + Helen(A, B, D) + Helen(A, C, D) + Helen(B, C, D);
		double radius = V4 / s;
		O /= s;
		return radius;
	}

};


