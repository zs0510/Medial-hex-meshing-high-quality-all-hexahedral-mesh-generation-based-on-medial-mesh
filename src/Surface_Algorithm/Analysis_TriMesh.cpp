#include "Analysis_TriMesh.h"

Vex Analysis_TriMesh::getCircleCenter(FH fh) {
	auto vhs = mesh.faces(fh).getVertexHandle();
	auto& v1 = mesh.vertices(vhs[0]);
	auto& v2 = mesh.vertices(vhs[1]);
	auto& v3 = mesh.vertices(vhs[2]);
	auto vec12 = v2 - v1, vec13 = v3 - v1;
	auto vec12X13 = vec12 % vec13;
	auto vec12X13X12 = vec12X13 % vec12;
	auto vec13X12X13 = vec13 % vec12X13;
	auto ret = v1 + (vec13X12X13 * vec12.norm2() + vec12X13X12 * vec13.norm2()) / (2.f * vec12X13.norm2());
	return ret;
}

void Analysis_TriMesh::calculateGaussCurvature(std::vector<double>& curvatures) {
	
	int vcnt = mesh.vsize(), fcnt = mesh.fsize();
	curvatures.resize(vcnt);
	std::vector<Vex> circum_center(fcnt);

	// isBotuse[i] = pair: face_i is botuse triangle, pair.first is obtuse vertex, pair.second is area
	std::unordered_map<int, std::pair<int, double>> isBotuse;

//#pragma omp parallel for
	for (int i = 0; i < fcnt; ++i) {
		FH fh(i);
		auto& face = mesh.faces(fh);
		auto vhs = face.getVertexHandle();
		auto& pos1 = mesh.vertices(vhs[0]);
		auto& pos2 = mesh.vertices(vhs[1]);
		auto& pos3 = mesh.vertices(vhs[2]);
		auto vec12 = pos2 - pos1;
		auto vec13 = pos3 - pos1;
		double area = 0.5f * std::fabs((vec12 % vec13).norm());
		vec12 = vec12.normalized();
		vec13 = vec13.normalized();
		double theta1 = std::acos(vec12 * vec13);
		if (theta1 > M_PI / 2.f) {
			isBotuse[fh] = std::pair<int, double>(vhs[0], area);
		} else {
			auto vec21 = vec12 * (-1);
			auto vec23 = (pos3 - pos2).normalized();
			double theta2 = std::acos(vec21 * vec23);
			if (theta2 > M_PI / 2.f) {
				isBotuse[fh] = std::pair<int, double>(vhs[1], area);
			} else {
				auto vec31 = vec13 * (-1);
				auto vec32 = vec23 * (-1);
				double theta3 = std::acos(vec31 * vec32);
				if (theta3 > M_PI / 2.f) {
					isBotuse[fh] = std::pair<int, double>(vhs[2], area);
				} else {
					// is non-botuse, cal circum center
					circum_center[fh] = getCircleCenter(fh);
				}
			}
		}

	}

	// cal curvature
//#pragma omp parallel for
	for (int i = 0; i < vcnt; ++i) {
		VH vh(i);
		auto& v = mesh.vertices(vh);
		double area_mixed = 0.f;
		double angle_sum = 0.f;
		for (auto fh : mesh.NeighborFh(vh)) {
			auto vhs = (mesh.faces(fh)).getVertexHandle();
			VH vh1(-1), vh2(-1);
			for (auto i : vhs) {
				if (i == vh) continue;
				if (vh1 == -1) vh1 = i;
				else vh2 = i;
			}
			assert(vh1 != -1 && vh2 != -1);
			auto& v1 = mesh.vertices(vh1);
			auto& v2 = mesh.vertices(vh2);

			// cal area
			if (isBotuse.count(fh)) {
				auto p = isBotuse[fh];
				if (p.first == vh) area_mixed += (p.second * 0.5f);
				else area_mixed += (p.second * 0.25f);
			} else {
				auto p1 = (v + v1) * 0.5f;
				auto p2 = (v + v2) * 0.5f;
				area_mixed += 0.5f * (p1 - circum_center[fh]).norm() * (p1 - v).norm();
				area_mixed += 0.5f * (p2 - circum_center[fh]).norm() * (p2 - v).norm();
			}

			// cal angle
			auto vec1 = (v1 - v).normalized();
			auto vec2 = (v2 - v).normalized();
			angle_sum += std::acos(vec1 * vec2);
		}
		assert(area_mixed > 0.f);
		curvatures[vh] = (D_PI - angle_sum) / area_mixed;
		//std::cout << "vh " << vh << ": curvature = " << curvatures[vh] << std::endl;
	}

}

void Analysis_TriMesh::calculateTriangleQuality(std::vector<double>& face_quality) {
	face_quality.reserve(mesh.fsize());
	for (auto& fp : mesh.allfaces()) {
		auto vhs = fp.second.getVertexHandle();
		double min_angle = 180.f;
		for (int i = 0; i < 3; ++i) {
			auto& v0 = mesh.vertices(vhs[i]);
			auto& v1 = mesh.vertices(vhs[(i + 1) % 3]);
			auto& v2 = mesh.vertices(vhs[(i + 2) % 3]);
			auto vec01 = (v1 - v0).normalized();
			auto vec02 = (v2 - v0).normalized();
			double angle = std::acos(vec01 * vec02) / M_PI * 180;
			min_angle = std::min(min_angle, angle);
		}
		double val = std::fmin(1, (60 - min_angle) / 60);
		//std::cout << "val: " << val << std::endl;
		face_quality.push_back(val);
	}

}