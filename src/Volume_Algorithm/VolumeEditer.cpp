#include "VolumeEditer.h"


double VolumeEditer::get_min_sj(MeshKernel::VolumeMesh& hexmesh, std::vector<Eigen::Vector3d> vertices) {

	if (vertices.size() != 8) {
		std::cerr << "Error: support hex element only!!! ___ get_min_sj ___ \n";
		return -1;
	}

	// HexaLab Quality
	std::vector<Eigen::Vector3f> points;
	for (auto& v : vertices) {
		Eigen::Vector3f p(v.x(), v.y(), v.z());
		points.emplace_back(p);
	}
	double msj = HexaLab::QualityMeasureFun::scaled_jacobian(points[0], points[1], points[2], points[3], 
		points[4], points[5], points[6], points[7], "");

	return msj;

	//// Whale-Sea Gu's method
	//double jacobian = 2.0;
	//double current_jacobian = 0.0;

	//Eigen::Vector3d L_0, L_2, L_3;

	//L_0 = (vertices[1] - vertices[0]).normalized();
	//L_2 = (vertices[3] - vertices[0]).normalized();// P_3 - P_0
	//L_3 = (vertices[4] - vertices[0]).normalized();// P_4 - P_0

	//current_jacobian = (L_2.cross(L_3)).dot(L_0);
	//if (current_jacobian < jacobian) jacobian = current_jacobian;

	//L_0 = (vertices[2] - vertices[1]).normalized();// P_2 - P_1
	//L_2 = (vertices[0] - vertices[1]).normalized();// P_0 - P_1
	//L_3 = (vertices[5] - vertices[1]).normalized();// P_5 - P_1

	//current_jacobian = (L_2.cross(L_3)).dot(L_0);
	//if (current_jacobian < jacobian) jacobian = current_jacobian;

	//L_0 = (vertices[3] - vertices[2]).normalized();// P_3 - P_2
	//L_2 = (vertices[1] - vertices[2]).normalized();// P_1 - P_2
	//L_3 = (vertices[6] - vertices[2]).normalized();// P_6 - P_2

	//current_jacobian = (L_2.cross(L_3)).dot(L_0);
	//if (current_jacobian < jacobian) jacobian = current_jacobian;

	//L_0 = (vertices[0] - vertices[3]).normalized();// P_0 - P_3
	//L_2 = (vertices[2] - vertices[3]).normalized();// P_2 - P_3
	//L_3 = (vertices[7] - vertices[3]).normalized();// P_7 - P_3

	//current_jacobian = (L_2.cross(L_3)).dot(L_0);
	//if (current_jacobian < jacobian) jacobian = current_jacobian;

	//L_0 = (vertices[7] - vertices[4]).normalized();// P_7 - P_4
	//L_2 = (vertices[5] - vertices[4]).normalized();// P_5 - P_4
	//L_3 = (vertices[0] - vertices[4]).normalized();// P_0 - P_4

	//current_jacobian = (L_2.cross(L_3)).dot(L_0);
	//if (current_jacobian < jacobian) jacobian = current_jacobian;

	//L_0 = (vertices[4] - vertices[5]).normalized();// P_4 - P_5
	//L_2 = (vertices[6] - vertices[5]).normalized();// P_6 - P_5
	//L_3 = (vertices[1] - vertices[5]).normalized();// P_1 - P_5

	//current_jacobian = (L_2.cross(L_3)).dot(L_0);
	//if (current_jacobian < jacobian) jacobian = current_jacobian;

	//L_0 = (vertices[5] - vertices[6]).normalized();// P_5 - P_6
	//L_2 = (vertices[7] - vertices[6]).normalized();// P_7 - P_6
	//L_3 = (vertices[2] - vertices[6]).normalized();// P_2 - P_6

	//current_jacobian = (L_2.cross(L_3)).dot(L_0);
	//if (current_jacobian < jacobian) jacobian = current_jacobian;

	//L_0 = (vertices[6] - vertices[7]).normalized();// P_6 - P_7
	//L_2 = (vertices[4] - vertices[7]).normalized();// P_4 - P_7
	//L_3 = (vertices[3] - vertices[7]).normalized();// P_3 - P_7

	//current_jacobian = (L_2.cross(L_3)).dot(L_0);
	//if (current_jacobian < jacobian) jacobian = current_jacobian;

	//return jacobian;


	//// Zhang Sheng's method
	//Eigen::Vector3d center(0, 0, 0);
	//for (auto& v : vertices) {
	//	center = center + v;
	//}
	//center /= 8.0;

	//double msj = 2;
	//int swap_cnt = 0;

	//for (int i = 0; i < 8; ++i) {

	//	int pre, next, oppo;

	//	if (i < 4) {
	//		pre = (i + 3) % 4;
	//		next = (i + 1) % 4;
	//		oppo = i + 4;

	//	} else {
	//		pre = i - 1;
	//		if (pre == 3) pre = 7;
	//		next = i + 1;
	//		if (next == 8) next = 4;
	//		oppo = i - 4;
	//	}
	//	// v1 v2 v3 形成的三角法向量应与 (center - v0) 同向
	//	auto v1 = vertices[pre];
	//	auto v2 = vertices[next];
	//	auto v3 = vertices[oppo];
	//	//auto center = (v1 + v2 + v3) / 3.0;
	//	auto vec = (center - vertices[i]).normalized();
	//	auto vec12 = v2 - v1, vec13 = v3 - v1;
	//	auto normal = (vec12.cross(vec13)).normalized();
	//	double cosine = vec.dot(normal);
	//	if (cosine < 0) {
	//		std::swap(v2, v3);
	//		swap_cnt++;
	//	}
	//	std::vector<Eigen::Vector3d> ev(3);
	//	ev[0] = (v1 - vertices[i]).normalized();
	//	ev[1] = (v2 - vertices[i]).normalized();
	//	ev[2] = (v3 - vertices[i]).normalized();
	//	Eigen::Matrix3d J;
	//	J << ev[0][0], ev[1][0], ev[2][0],
	//		ev[0][1], ev[1][1], ev[2][1],
	//		ev[0][2], ev[1][2], ev[2][2];
	//	//std::cout << J << endl;
	//	/*
	//	For an element, the full range of the Scaled Jacobian value is from −1 to + 1.
	//	And fscaled_jacobian = 1 if the element is an ideal element,
	//	fscaled_jacobian = −1 if the element is a worst distorted element
	//	*/
	//	msj = std::min(msj, J.determinant());

	//}
	//std::cout << "msj = " << msj << ", swap_cnt = " << swap_cnt << "\n";
	//if (swap_cnt == 8) return -1.0;
	//return msj;

}

double VolumeEditer::get_grafting_energy(std::vector<Eigen::Vector3d> vertices) {

	double energy = std::numeric_limits<double>::max();

	if (vertices.size() != 8) return energy;

	double alpha = 0.6, beta = 0.4;

	Eigen::Vector3d c1 = (vertices[0] + vertices[1] + vertices[2] + vertices[3]) / 4.0;
	Eigen::Vector3d c2 = (vertices[4] + vertices[5] + vertices[6] + vertices[7]) / 4.0;
	Eigen::Vector3d dir = (c2 - c1).normalized();

	// calculate the distance energy
	double energy_dis = 0;
	for (int i = 0; i < 4; ++i) {
		energy_dis += (vertices[i + 4] - vertices[i]).squaredNorm();
	}
	energy_dis *= alpha;

	// calculate the direction energy
	double energy_dir = 0;
	for (int i = 0; i < 4; ++i) {
		energy_dir += ((vertices[i + 4] - vertices[i]).normalized()).dot(dir);
	}
	energy_dir = beta / energy_dir;

	std::cout << "Grafting energy: energy_distance = " << energy_dis << ", energy_dir = " << energy_dir << ", the dis/dir = " << (energy_dis / energy_dir) << "\n";

	return energy_dis + energy_dir;

}

void VolumeEditer::add_cell(MeshKernel::VolumeMesh& hexmesh, FH fh1, FH fh2) {

	if (!hexmesh.isValid(fh1) || !hexmesh.isValid(fh2)) {
		std::cerr << "Error: input invalid face handle!!!\n";
		return;
	}

	if (hexmesh.isConnected(fh1, fh2)) {// 通过判断是否有相邻的边来判断是否连接
		add_cell_connected(hexmesh, fh1, fh2);
	} else {
		add_cell_disconnected(hexmesh, fh1, fh2);
	}

}

void VolumeEditer::add_cell_connected(MeshKernel::VolumeMesh& hexmesh, FH fh1, FH fh2) {

	// 找到公共边
	EH common_eh(-1);
	auto face1 = hexmesh.faces(fh1), face2 = hexmesh.faces(fh2);
	auto f_ehs1 = face1.getEdgeHandle(), f_ehs2 = face2.getEdgeHandle();
	for (auto& f_eh1 : f_ehs1) {
		for (auto& f_eh2 : f_ehs2) {
			if (f_eh2 == f_eh1) {
				common_eh = f_eh1;
			}
		}
		if (common_eh != -1) break;
	}

	if (common_eh == -1) {
		std::cerr << "Error: the facet no common edge!!!\n";
		return;
	}

	/*std::vector<VH> common_vhs;
	auto common_edge = hexmesh.edges(common_eh);
	common_vhs.emplace_back(common_edge.vh1());
	common_vhs.emplace_back(common_edge.vh2());*/

	std::vector<std::vector<VH>> vhs{ face1.getVertexHandle() , face2.getVertexHandle() };
	auto& common_edge = hexmesh.edges(common_eh);
	std::vector<VH> common_vhs{ common_edge.vh1() , common_edge.vh2() };
	std::vector<std::vector<VH>> adjvhs(2);
	for (int i = 0; i < 2; ++i) {
		for (auto& vh : vhs[i]) {
			if (vh == common_vhs[0] || vh == common_vhs[1]) continue;
			if (hexmesh.isConnected(vh, common_vhs[0])) adjvhs[0].push_back(vh);
			if (hexmesh.isConnected(vh, common_vhs[1])) adjvhs[1].push_back(vh);
		}
	}
	assert(adjvhs[0].size() == 2 && adjvhs[1].size() == 2 && " adjvhs size is wrong! ");
	std::vector<VH> new_vhs(2, VH(-1));
	for (int i = 0; i < 2; ++i) {
		auto& v0 = hexmesh.vertices(common_vhs[i]);
		auto& v1 = hexmesh.vertices(adjvhs[i][0]);
		auto& v2 = hexmesh.vertices(adjvhs[i][1]);
		auto vec = (v1 - v0) + (v2 - v0);
		auto new_v = v0 + vec;
		new_vhs[i] = hexmesh.AddVertex(new_v);
	}
	std::vector<VH> new_cell = { common_vhs[0], adjvhs[0][0], new_vhs[0], adjvhs[0][1], common_vhs[1], adjvhs[1][0], new_vhs[1], adjvhs[1][1] };
	auto new_ch = hexmesh.AddCell(new_cell);



}

void VolumeEditer::add_cell_disconnected(MeshKernel::VolumeMesh& hexmesh, FH fh1, FH fh2) {

	

	auto& face1 = hexmesh.faces(fh1);
	auto& face2 = hexmesh.faces(fh2);
	const auto& f_vhs2 = face2.getVertexHandle();
	std::vector<VH> hex_vhs = face1.getVertexHandle();
	hex_vhs.insert(hex_vhs.end(), f_vhs2.begin(), f_vhs2.end());
	if (hex_vhs.size() != 8) {
		std::cerr << "Error: support hex only!!!\n";
		return;
	}
	// 调整 hex_vhs 的顺序
	Vex face_center = hexmesh.getFaceCenter(fh1), cell_center(0, 0, 0);
	for (auto& vh : hex_vhs) {
		cell_center = cell_center + hexmesh.vertices(vh);
	}
	cell_center /= 8.0;
	Vec dir = (cell_center - face_center).normalized();
	Vec normal = ((hexmesh.vertices(hex_vhs[1]) - hexmesh.vertices(hex_vhs[0])).cross(hexmesh.vertices(hex_vhs[2]) - hexmesh.vertices(hex_vhs[0]))).normalized();
	if (dir.dot(normal) < 0) {// 保证第一个面的法向量是朝向体质心的
		std::reverse(hex_vhs.begin(), hex_vhs.begin() + 4);
		std::reverse(hex_vhs.begin() + 4, hex_vhs.end());
	}

	std::vector<std::vector<int>> total_mode = {
		{ 0, 1, 2, 3, 4, 5, 6, 7 }, { 0, 1, 2, 3, 4, 7, 6, 5 },
		{ 0, 1, 2, 3, 5, 4, 7, 6 }, { 0, 1, 2, 3, 5, 6, 7, 4 },
		{ 0, 1, 2, 3, 6, 5, 4, 7 }, { 0, 1, 2, 3, 6, 7, 4, 5 },
		{ 0, 1, 2, 3, 7, 4, 5, 6 }, { 0, 1, 2, 3, 7, 6, 5, 4 },
	};

	/*double max_msj = -1;
	int best_mode = -1;

	for (int i = 0; i < 8; ++i) {
		std::vector<Eigen::Vector3d> positions;
		for (auto& idx : total_mode[i]) {
			auto& v = hexmesh.vertices(hex_vhs[idx]);
			Eigen::Vector3d pos(v.x(), v.y(), v.z());
			positions.emplace_back(pos);
		}
		double msj = get_min_sj(hexmesh, positions);
		if (msj > max_msj) {
			max_msj = msj;
			best_mode = i;
		}
	}

	if (max_msj > 0) {
		std::vector<VH> hexvhs;
		for (auto& idx : total_mode[best_mode]) {
			hexvhs.emplace_back(hex_vhs[idx]);
		}
		hexmesh.AddCell(hexvhs);
		std::cout << "Add cell success. The scaled Jocabian of new cell is " << max_msj << "\n";
	} else {
		std::cerr << "Error: opertor is cancel because it will generate a invert cell!!!\t msj = " << max_msj << ".\n";
	}*/


	//// From Paper: Filling Triangular Mesh Model with All-Hex Mesh by Volume Subdivision Fitting
	//double min_energy = std::numeric_limits<double>::max();
	//int best_mode = -1;
	//for (int i = 0; i < 8; ++i) {
	//	std::vector<Eigen::Vector3d> positions;
	//	for (auto& idx : total_mode[i]) {
	//		auto& v = hexmesh.vertices(hex_vhs[idx]);
	//		Eigen::Vector3d pos(v.x(), v.y(), v.z());
	//		positions.emplace_back(pos);
	//	}
	//	double energy = get_grafting_energy(positions);
	//	if (energy < min_energy) {
	//		min_energy = energy;
	//		best_mode = i;
	//	}
	//}


	// From Sheng
	
	// 求得第二个面的法向量
	Vec normal_f2 = ((hexmesh.vertices(hex_vhs[5]) - hexmesh.vertices(hex_vhs[4])).cross((hexmesh.vertices(hex_vhs[6]) - hexmesh.vertices(hex_vhs[4])))).normalized();
	double a = normal_f2.x(), b = normal_f2.y(), c = normal_f2.z();
	Vex center_f2 = hexmesh.getFaceCenter(fh2);
	double d = -(a * center_f2.x() + b * center_f2.y() + c * center_f2.z());

	// 求出四个交点
	std::vector<Vex> intersections(4);
	for (int i = 0; i < 4; ++i) {
		auto srcv = hexmesh.vertices(hex_vhs[i]);
		double t = -(a * srcv.x() + b * srcv.y() + c * srcv.z() + d)
			/ (a * dir.x() + b * dir.y() + c * dir.z());
		intersections[i] = srcv + dir * t;
		/*std::cout << "Intersection[" << i << "]: (" << intersections[i].x() << ", " << intersections[i].y() << ", "
			<< intersections[i].z() << ").\n";*/
	}

	double min_energy = std::numeric_limits<double>::max();
	int best_mode = -1;
	for (int i = 0; i < 8; ++i) {
		
		double energy = 0;
		for (int j = 0; j < 4; ++j) {
			Vec diff = intersections[j] - hexmesh.vertices(hex_vhs[total_mode[i][j + 4]]);
			
			energy += diff.norm2();
		}

		if (energy < min_energy) {
			min_energy = energy;
			best_mode = i;
		}

		//std::cout << "\ti = " << i << ", energy = " << energy << std::endl;

	}



	if (best_mode != -1) {
		std::vector<VH> hexvhs;
		for (auto& idx : total_mode[best_mode]) {
			hexvhs.emplace_back(hex_vhs[idx]);
		}
		hexmesh.AddCell(hexvhs);
		std::cout << "Add cell success. The grafting energy of new cell is " << min_energy << ".\n";
	} else {
		std::cerr << "Error: opertor is cancel because it will generate a invert cell!!!\t msj = " << min_energy << ".\n";
	}
	

}