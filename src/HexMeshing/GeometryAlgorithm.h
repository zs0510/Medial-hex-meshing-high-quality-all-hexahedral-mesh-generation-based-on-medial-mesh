#pragma once
#include <vector>
#include "Kernel/Mesh.h"
#include "Tools/MeshMath.h"

using namespace std;

namespace GeometryAlgorithm {

	void normal_matching(vector<Vex>& vertices, const vector<vector<int>>& v_neigbor_faces, const vector<Vex>& faces_center, const vector<Vec>& faces_normal) {

		int vcnt = vertices.size();
		int fcnt = faces_center.size();

		if (v_neigbor_faces.size() != vcnt || faces_normal.size() != fcnt) {
			std::cerr << "Size not compatible!!!" << std::endl;
			return;
		}

		for (int it = 0; it < 10; ++it) {
			for (int vh = 0; vh < vcnt; ++vh) {
				Vec move(0, 0, 0);
				auto& adj_faces = v_neigbor_faces[vh];
				for (auto& fh : adj_faces) {
					double wj = faces_normal[fh].dot(faces_center[fh] - vertices[vh]);
					move += (faces_normal[fh] * wj);
				}
				move /= adj_faces.size();
				vertices[vh] += move;
			}
		}

	}

	bool get_intersection_triangle(const Vector3d& src, const Vector3d& dir, const vector<Vector3d>& triangle, Vector3d& res) {

		if (triangle.size() != 3) {
			std::cerr << "Your triangle is non-triangle!!!" << std::endl;
			return false;
		}

		// Refer to https://www.cnblogs.com/graphics/archive/2010/08/09/1795348.html
		double t, u, v;
		Vector3d E1 = triangle[1] - triangle[0];
		Vector3d E2 = triangle[2] - triangle[0];
		// ≤Ê≥À
		Vector3d P = dir.Cross(E2);
		// µ„≥À
		double det = P.Dot(E1);
		Vector3d T;
		if (det > 0) {
			T = src - triangle[0];
		}
		else {
			T = triangle[0] - src;
			det = -det;
		}
		// Œ¥œ‡Ωª
		if (det < 0)
			return false;

		// Calculate u and make sure u <= 1
		u = T.Dot(P);
		if (u < 0.0f || u > det)
			return false;

		// Q
		Vector3d Q = T.Cross(E1);

		// Calculate v and make sure u + v <= 1
		v = dir.Dot(Q);
		if (v < 0.0f || u + v > det)
			return false;

		//// Calculate t, scale parameters, ray intersects triangle
		//t = E2.x() * Q.x() + E2.y() * Q.y() + E2.z() * Q.z();
		float fInvDet = 1.0f / det;
		//t *= fInvDet;
		u *= fInvDet;
		v *= fInvDet;
		res = u * E1 + v * E2 + triangle[0];
		Vector3d dir1 = res - src;
		dir1.Normalize();
		if (dir1.Dot(dir) >= 0.f) {
			//isec.happened = true;
			//isec.distance = (isec.pos - ray.pos).Length();
			return true;
		}

		return false;

	}

	void print_quality_msj(MeshKernel::VolumeMesh& hexmesh) {

		int ccnt = hexmesh.csize();
		vector<int> counts(11);
		double mj_min = 2.0;
		double mj_max = -2.0;
		double mj_sum = 0.0;
		std::cout << "***** SJ begin *****\n";
		for (auto& cp : hexmesh.allcells()) {
			double mj = MeshMath::get_min_scaled_Jacobian(hexmesh, cp.first);
			mj_sum += mj;
			mj_min = std::min(mj_min, mj);
			mj_max = std::max(mj_max, mj);
			for (int i = 0; i <= 10; ++i) {
				if (mj <= 0.1 * i) {
					counts[i]++;
					break;
				}
			}
			printf("%.3f, ", mj);
		}
		std::cout << "\n\n***** SJ end *****\n\n";
		std::cout << "\n\n[HexMesh Quality]: mj_min = " << mj_min << ", mj_max = " << mj_max << ", mj_avg = " << (mj_sum / ccnt) << std::endl;
		for (int i = 0; i <= 10; ++i) {
			if (i == 0) std::cout << "[-1.0, 0): ";
			else std::cout << "[" << 0.1 * (i - 1) << ", " << 0.1 * i << "): ";

			//std::cout << "Count = " << counts[i] << ", ratio = " << counts[i] * 1.0 / ccnt << std::endl;
			printf("\tcount = %d, \tratio = %.2f%%\n", counts[i], counts[i] * 100.0 / ccnt);
		}

	}

};