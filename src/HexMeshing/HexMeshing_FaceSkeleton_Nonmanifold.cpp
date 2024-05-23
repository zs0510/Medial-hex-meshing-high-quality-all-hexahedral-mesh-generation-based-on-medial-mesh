#include "HexMeshing_FaceSkeleton_Nonmanifold.h"

void HexMeshing_FaceSkeleton_Nonmanifold::hexmeshing() {

	bool input_ok = false;// 获取面片域的个数
	int dcnt = QInputDialog::getInt(nullptr, "", "Please input domain size",
		1, 0, 30, 1, &input_ok);
	if (!input_ok) return;

	//Project_KD_Tree_SkeletalMesh* kdtree_sklmesh = new Project_KD_Tree_SkeletalMesh(sklmesh);// 建立KD树, 从而可以获取任一位置的半径

	MeshKernel::iGameKdTreeMedialMesh* kdtree_sklmesh = new MeshKernel::iGameKdTreeMedialMesh(sklmesh);

	for (int di = 0; di < dcnt; ++di) {

		//std::vector<std::pair<Vex, Vex>> pos2_pos3;// 保存二维到三维点的映射关系

		//// 1. 读入映射关系
		//{
		//	std::string line;
		//	std::string filename = "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(di) + ".mapping.txt";
		//	std::ifstream mapping_inputfile(filename, std::ios::in);

		//	int pair_count = 0;
		//	{
		//		line.clear();
		//		getline(mapping_inputfile, line);
		//		std::stringstream linestream;
		//		linestream.str(line);
		//		linestream >> pair_count;
		//		std::cout << "Reading " << pair_count << " pairs position.\n";
		//	}
		//	for (int i = 0; i < pair_count; ++i) {
		//		line.clear();
		//		getline(mapping_inputfile, line);
		//		std::stringstream linestream;
		//		linestream.str(line);
		//		double x3, y3, z3, x2, y2, z2;
		//		linestream >> x3 >> y3 >> z3 >> x2 >> y2 >> z2;
		//		Vex pos3(x3, y3, z3);
		//		Vex pos2(x2, y2, z2);
		//		pos2_pos3.emplace_back(pos2, pos3);
		//	}
		//}

		//// 2. 读入四边形网格

		//std::string filename_2d = "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(di) + ".2d.obj";
		//Project_KD_Tree_2D_To_3D kd_tree(filename_2d, pos2_pos3);//.负责将二维顶点映射回三维

		//MeshKernel::IO io;
		//MeshKernel::SurfaceMesh quad_mesh // HexMesh_td 输出的四边形网格
		//	= io.ReadOffFile("C:\\My Files\\Graphics\\model_data\\!test_data\\quad_mesh" + std::to_string(di) + ".off");

		//for (auto& vp : quad_mesh.allvertices()) {
		//	auto& v = quad_mesh.vertices(vp.first);
		//	v = kd_tree.get_project_pos3(v);// 将 2D 坐标转换到 3D
		//}

		//MeshKernel::IO io_3d_writer;
		//io_3d_writer.WriteObjFile(quad_mesh, "C:\\My Files\\Graphics\\model_data\\!test_data\\quad_mesh" + std::to_string(di) + "3d.obj");

		MeshKernel::IO io;
		MeshKernel::SurfaceMesh quad_mesh
			= io.ReadOffFile("C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(di) + "_quad.off");

		std::vector<Vec> quad_normals;
		normal_filtering(quad_mesh, quad_normals);// 已是三维空间的网格

		// 3. 初始化每个点的半径
		std::vector<double> quad_radius(quad_mesh.vsize());
		for (auto& vp : quad_mesh.allvertices()) {
			quad_radius[vp.first] = kdtree_sklmesh->get_radius(vp.second) * radius_scale;
		}

		// 4. 拉伸出六面体网格
		std::vector<VH> up_vhs(quad_mesh.vsize());
		std::vector<VH> down_vhs(quad_mesh.vsize());

		for (auto& vp : quad_mesh.allvertices()) {
			Vec vn = quad_normals[vp.first];
			up_vhs[vp.first] = hexmesh.AddVertex(vp.second + vn * quad_radius[vp.first] * 0.5);
			down_vhs[vp.first] = hexmesh.AddVertex(vp.second - vn * quad_radius[vp.first] * 0.5);
		}
		for (auto& quad_fp : quad_mesh.allfaces()) {
			auto vhs = quad_fp.second.getVertexHandle();
			std::vector<VH> up_face, down_face;
			for (auto& vh : vhs) {
				up_face.push_back(up_vhs[vh]);
				down_face.push_back(down_vhs[vh]);
			}
			up_face.insert(up_face.end(), down_face.begin(), down_face.end());
			hexmesh.AddCell(up_face);
		}


	}

}

void HexMeshing_FaceSkeleton_Nonmanifold::normal_filtering(MeshKernel::SurfaceMesh& quad_mesh, vector<Vec>& quad_normals) {

	int vcnt = quad_mesh.vsize();

	std::vector<Vex> positions(vcnt);
	for (auto& vp : quad_mesh.allvertices()) {
		positions[vp.first] = vp.second;
	}

	vector<Vector3f> abtree_vs;
	for (auto& fp : sklmesh.allfaces()) {
		const auto& f_vhs = fp.second.getVertexHandle();
		for (int i = 2; i < f_vhs.size(); ++i) {
			vector<VH> trivhs = { f_vhs[0], f_vhs[i - 1], f_vhs[i] };
			for (auto& vh : trivhs) {
				auto& v = sklmesh.vertices(vh);
				abtree_vs.emplace_back(Vector3f(v.x(), v.y(), v.z()));
			}
		}
	}
	AABB_Tree* abtree = new AABB_Tree(abtree_vs);

	for (int it = 0; it < 15; ++it) {
		std::vector<Vex> pre_positions = positions;
#pragma omp parallel for
		for (int vid = 0; vid < vcnt; ++vid) {
			VH vh(vid);
			if (quad_mesh.isOnBoundary(vh)) continue;// 边界顶点不更新
			Vec pos(0, 0, 0);
			const auto& adjvhs = quad_mesh.NeighborVh(vh);
			if (adjvhs.empty()) continue;
			for (auto& adjvh : adjvhs) {
				pos += pre_positions[adjvh];
			}
			pos /= adjvhs.size();
			Vector3f _pos(pos.x(), pos.y(), pos.z());
			Vector3f nearest;
			abtree->findNearstPoint(_pos, nearest);
			positions[vh] = Vex(nearest.x(), nearest.y(), nearest.z());
			//positions[vh] = pos;
		}
	}

	for (auto& vp : quad_mesh.allvertices()) {
		auto& v = quad_mesh.vertices(vp.first);
		v = positions[vp.first];
	}


	quad_mesh.genAllVerticesNormal();
	quad_normals.resize(quad_mesh.vsize());
	for (auto& vp : quad_mesh.allvertices()) {
		quad_normals[vp.first] = Vec(vp.second.getNormalX(), vp.second.getNormalY(), vp.second.getNormalZ());
	}

	//std::vector<Vec> pre_quad_normals;
	//for (int it = 0; it < 3; ++it) {
	//	pre_quad_normals = quad_normals;
	//	for (int i = 0; i < quad_normals.size(); ++i) {
	//		VH vh(i);
	//		//if (quad_mesh.isOnBoundary(vh)) continue;// 边界顶点法向量不更新
	//		Vec normal(0, 0, 0);
	//		const auto& adjvhs = quad_mesh.NeighborVh(vh);
	//		if (adjvhs.empty()) continue;
	//		for (auto& adjvh : adjvhs) {
	//			normal += pre_quad_normals[adjvh];
	//		}
	//		normal /= adjvhs.size();
	//		normal.normalize();
	//		quad_normals[i] = normal;
	//	}
	//}

}