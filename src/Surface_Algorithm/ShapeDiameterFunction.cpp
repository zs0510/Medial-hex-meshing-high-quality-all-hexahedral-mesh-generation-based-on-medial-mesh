#include "ShapeDiameterFunction.h"

void ShapeDiameterFunction::get_faces_sdf(std::vector<double>& FaceDiameter, bool smooth_flag, bool normalize_flag) {

	calc_faces_sdf();

	if (smooth_flag) {
		calc_sigma();
		smooth_sdf();
	}

	if (normalize_flag) {
		normalize_sdf();
		FaceDiameter = normalized_faces_sdf;
	} else {
		FaceDiameter = faces_sdf;
	}
	
}

void ShapeDiameterFunction::calc_faces_sdf() {

	time_t beg = clock();

	mesh.genAllFacesNormal();
	int fcnt = mesh.fsize();
	faces_sdf.resize(fcnt);

#pragma omp parallel for
	for (int fid = 0; fid < fcnt; ++fid) {
		FH fh(fid);
		auto& face = mesh.faces(fh);
		Vex center = mesh.getFaceCenter(fh);
		Vec normal = face.getNormal();
		std::vector<Vex> sample_points;
		sample_on_circle(center, normal, sample_points);
		if (sample_points.empty()) continue;
		double length_sum = 0, weight_sum = 0;
		Vex origin = center - normal * 1e-8;// 向内移动，避免出现精度问题
		for (auto& sp : sample_points) {
			Vec dir = (sp - origin).normalized();
			Ray ray(Vector3d(origin.x(), origin.y(), origin.z()), Vector3d(dir.x(), dir.y(), dir.z()));
			Intersection inter = bvh_tree.getIntersection(ray);
			if (inter.happened) {
				FH inter_fh(inter.index);
				Vec inter_normal = mesh.faces(inter_fh).getNormal();
				if (inter_normal.dot(normal) < 0.f) {
					// For each such ray we check the normal at the point of intersection,
					//	and ignore intersections where the normal at the intersection points in the same direction as the point - of - origin
					//	(the same direction is defined as an angle difference less than 90◦).
					double angle = std::acos(dir.dot(normal * -1));
					double weight = 1.0 / angle;
					length_sum += weight * inter.distance;
					weight_sum += weight;
				}
			}
		}

		faces_sdf[fid] = length_sum / weight_sum;

	}

	time_t end = clock();

	printf("Calculate SDF success, cost %dms\n", int(end - beg));

}

void ShapeDiameterFunction::calc_sigma() {

	faces_center.resize(mesh.fsize());

#pragma omp parallel for
	for (int fid = 0; fid < mesh.fsize(); ++fid) {
		FH fh(fid);
		faces_center[fid] = mesh.getFaceCenter(fh);
	}

	int cnt = 0;
	for (auto& fp : mesh.allfaces()) {
		auto fhs = mesh.NeighborFh(fp.first);// share common edge
		for (auto fh : fhs) {
			sigma_center += (faces_center[fp.first] - faces_center[fh]).norm();
			sigma_sdf += std::fabs(faces_sdf[fp.first] - faces_sdf[fh]);
			cnt++;
		}
	}
	sigma_center /= cnt;
	sigma_sdf /= cnt;

	printf("[Parameter]: Sigma Center = %.4f, Sigma SDF = %.4f\n", sigma_center, sigma_sdf);

}

void ShapeDiameterFunction::smooth_sdf(int iter_times) {

	time_t beg = clock();

	int fcnt = mesh.fsize();
	double double_sigmacenter2 = 2.0 * sigma_center * sigma_center;
	double double_sigmasdf2 = 2.0 * sigma_sdf * sigma_sdf;

	for (int it = 0; it < iter_times; ++it) {

		std::vector<double> new_faces_sdf(fcnt);

#pragma omp parallel for
		for (int fid = 0; fid < fcnt; ++fid) {
			FH fh(fid);
			double sdf_sum = 0, weight_sum = 0;
			const auto& adjfhs = mesh.NeighborFh(fh);
			for (auto& adjfh : adjfhs) {
				double delta_center = (faces_center[fh] - faces_center[adjfh]).norm();
				double delta_sdf = std::fabs(faces_sdf[fh] - faces_sdf[adjfh]);
				double weight_center = std::exp(-(delta_center * delta_center) / double_sigmacenter2);
				double weight_sdf = std::exp(-(delta_sdf * delta_sdf) / double_sigmasdf2);
				double weight = weight_center * weight_sdf;
				sdf_sum += faces_sdf[adjfh] * weight;
				weight_sum += weight;
			}
			new_faces_sdf[fh] = sdf_sum / weight_sum;
		}

		faces_sdf = new_faces_sdf;

	}

	time_t end = clock();

	printf("Smooth SDF success, cost %dms\n", int(end - beg));

}

void ShapeDiameterFunction::sample_on_circle(const Vex& center, const Vex& normal, std::vector<Vex>& sample_points) {

	// cone parameters: opening_angle = 120, rays_count = 40; 
	sample_points.clear();

	double d = -(normal.dot(center));
	Vex tmp(center.x() + 0.01, center.y() + 0.01, center.z() + 0.01);
	if (normal.z() != static_cast<double>(0.0)) {
		tmp.setZ( -(tmp.x() * normal.x() + tmp.y() * normal.y() + d) / normal.z() );
	} else if (normal.y() != static_cast<double>(0.0)) {
		tmp.setY( -(tmp.x() * normal.x() + tmp.z() * normal.z() + d) / normal.y() );
	} else if (normal.z() != static_cast<double>(0.0)) {
		tmp.setX( -(tmp.y() * normal.y() + tmp.z() * normal.z() + d) / normal.x() );
	} else {
		std::cerr << "[Error]: Normal is ZERO!!!\n";
		return;
	}
	Vec dir_u = (tmp - center).normalized();
	Vec dir_v = (normal.cross(dir_u)).normalized();
	Vec dir_aux0 = (dir_u + dir_v).normalized();
	Vec dir_aux1 = (dir_u - dir_v).normalized();

	std::vector<Vec> dirs = { dir_u, dir_v, dir_aux0, dir_aux1 };

	double radius = 1.73205081;// sqrt(3)

	//sample_points.emplace_back(center);// 圆心

	int sample_dense = 5;
	for (int i = 1; i <= sample_dense; ++i) {
		for (auto& dir : dirs) {
			sample_points.emplace_back(center + dir * radius * i / sample_dense);
			sample_points.emplace_back(center - dir * radius * i / sample_dense);
		}
	}

	for (auto& sp : sample_points) {// 形成120°的圆锥
		sp -= normal;
	}

}

void ShapeDiameterFunction::update_sdf_max_min() {

	sdf_min = std::numeric_limits<double>::max(), sdf_max = 0;
	for (auto& sdf : faces_sdf) {
		sdf_min = std::fmin(sdf, sdf_min);
		sdf_max = std::fmax(sdf, sdf_max);
	}

}

void ShapeDiameterFunction::normalize_sdf() {

	int fcnt = mesh.fsize();
	normalized_faces_sdf.resize(fcnt);
	update_sdf_max_min();
	double alpha = 4.0;// alpha is a normalizing parameter, which is set to 4
	double inverse_log_alpha = 1.0 / std::log(alpha + 1.0);
	double sdf_range = sdf_max - sdf_min;

#pragma omp parallel for
	for (int fid = 0; fid < fcnt; ++fid) {
		normalized_faces_sdf[fid] = std::log( (faces_sdf[fid] - sdf_min) / sdf_range * alpha + 1.0) * inverse_log_alpha;
		/*if (normalized_faces_sdf[fid] > 1.0) {
			std::cerr << "Error: sdf is invalid! ";
			std::cerr << "SDF_MAX = " << sdf_max << ", SDF_I = " << faces_sdf[fid] << std::endl;
		}*/
	}

}