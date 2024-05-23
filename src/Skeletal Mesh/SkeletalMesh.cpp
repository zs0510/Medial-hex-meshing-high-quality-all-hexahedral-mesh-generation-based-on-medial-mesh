#include "SkeletalMesh.h"

bool SkeletalMesh::readMaFile(std::string filename) {

	std::ifstream inputfile(filename, std::ios::in);
	std::string line;
	int vsz, esz, fsz;
	SkeletalMesh* new_mesh = new SkeletalMesh();

	std::cout << "Reading " << filename << " File" << std::endl;

	{
		std::stringstream strs;
		getline(inputfile, line);
		strs.str(line);
		strs >> vsz >> esz >> fsz;

	}
	
	char c;
	double x, y, z, r;
	int v0, v1, v2;

	for (int i = 0; i < vsz; ++i) {
		line.clear();
		getline(inputfile, line);
		std::stringstream strs;
		strs.str(line);
		strs >> c >> x >> y >> z >> r;
		MeshKernel::iGameVertex vex(x, y, z);
		int vsz_before = new_mesh->vsize();
		auto vh = new_mesh->AddVertex(vex);
		if (vh != vsz_before) std::cout << "Skeletal Mesh: exist vertices have same position.\n";
		new_mesh->radius[vh] = r;
	}

	for (int i = 0; i < esz; ++i) {
		line.clear();
		getline(inputfile, line);
		std::stringstream strs;
		strs.str(line);
		strs >> c >> v0 >> v1;
		new_mesh->AddEdge(VH(v0), VH(v1));
	}

	for (int i = 0; i < fsz; ++i) {
		line.clear();
		getline(inputfile, line);
		std::stringstream strs;
		strs.str(line);
		strs >> c >> v0 >> v1 >> v2;
		std::vector<VH> vhs = {VH(v0), VH(v1), VH(v2)};
		new_mesh->AddFace(vhs);
	}

	inputfile.close();
	*this = *new_mesh;

	std::cout << "Read Skeletal Mesh success. vertices = " << this->vsize() <<
		", edges = " << this->esize() << ", faces = " << this->fsize() << std::endl;

	return true;
}

bool SkeletalMesh::writeMaFile(std::string filename) {

	std::ofstream off(filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
		off.close();
		return false;
	}
	updateAllHandles();
	off << vsize() << " " << esize() << " " << fsize() << std::endl;
	for (int i = 0; i < vsize(); ++i) {
		VH vh(i);
		auto& v = vertices_[vh];
		off << "v " << v.x() << " " << v.y() << " " << v.z() << " " << radius[vh] << std::endl;
	}
	for (auto e : edges_) {
		off << "e " << e.second.vh1() << " " << e.second.vh2() << std::endl;
	}
	for (auto f : faces_) {
		const auto& vhs = f.second.getVertexHandle();
		off << "f " << vhs[0] << " " << vhs[1] << " " << vhs[2] << std::endl;
	}
	off.close();
	std::cout << "Writing file success. #V = " << vsize() << ", #E = " << esize() << ", #F = " << fsize() << std::endl;
	return true;
}


void SkeletalMesh::updateAllHandles() {// 使用多线程前最好调用此函数

	int vcnt = VertexSize(), fcnt = FaceSize();
	std::vector<MeshKernel::iGameVertex> newVertices;
	std::vector<std::vector<MeshKernel::iGameVertexHandle>> newEdges;
	std::vector<std::vector<MeshKernel::iGameVertexHandle>> newFaces;
	std::unordered_map<int, int> mp;// old id to new id
	std::vector<double> newRadius;
	int idx = 0;

	for (auto& ep : alledges()) {
		std::vector<MeshKernel::iGameVertexHandle> vhs = { ep.second.vh1(), ep.second.vh2() };
		for (auto& vh : vhs) {
			if (!mp.count(vh)) {
				mp[vh] = idx++;
				newVertices.push_back(vertices_[vh]);
				newRadius.push_back(radius[vh]);
			}
			vh = VH(mp[vh]);
		}
		newEdges.push_back(vhs);
	}

	for (auto& fp : allfaces()) {
		auto vhs = fp.second.getVertexHandle();
		for (auto& vh : vhs) {
			if (!mp.count(vh)) {
				mp[vh] = idx++;
				newVertices.push_back(vertices_[vh]);
				newRadius.push_back(radius[vh]);
			}
			vh = VH(mp[vh]);
		}
		newFaces.push_back(vhs);
	}

	SkeletalMesh new_mesh;
	for (int i = 0; i < newVertices.size(); ++i) {
		auto new_vh = new_mesh.AddVertex(newVertices[i]);
		new_mesh.radius[new_vh] = newRadius[i];
	}
	for (int i = 0; i < newEdges.size(); ++i) {
		new_mesh.AddEdge(newEdges[i][0], newEdges[i][1]);
	}
	for (int i = 0; i < newFaces.size(); ++i) {
		new_mesh.AddFace(newFaces[i]);
	}
	new_mesh.scale_factor = this->scale_factor;

	*this = new_mesh;

}


void SkeletalMesh::scale2Uint() {

	initBBox();

	scale_factor = 1 / (this->BBoxMax - this->BBoxMin).norm();

	for (auto vp : vertices_) {
		auto& v = this->vertices(vp.first);
		v *= scale_factor;
		this->radius[vp.first] *= scale_factor;
	}
	
	initBBox();
	std::cout << "Scale success. Diagonal length is " << (this->BBoxMax - this->BBoxMin).norm() << std::endl;
}

void SkeletalMesh::scale2Origin() {

	scale_factor = 1 / scale_factor;
	for (auto vp : vertices_) {
		auto& v = this->vertices(vp.first);
		v *= scale_factor;
		this->radius[vp.first] *= scale_factor;
	}
	
}

bool SkeletalMesh::readSkelFile(std::string filename) {

	std::ifstream inputfile(filename, std::ios::in);
	std::string line;
	int vsz, esz;
	SkeletalMesh* new_mesh = new SkeletalMesh();

	std::cout << "Reading " << filename << " File" << std::endl;

	{
		std::stringstream strs;
		getline(inputfile, line);// 首行, 列名
		line.clear();
		getline(inputfile, line);// 第二行, 顶点个数
		strs.str(line);
		strs >> vsz;
	}
	char c;
	double x, y, z, r;
	int vid, v0, v1, v2;
	vector<vector<int>> neighbors;
	for (int i = 0; i < vsz; ++i) {
		line.clear();
		getline(inputfile, line);
		std::stringstream strs;
		strs.str(line);
		strs >> vid >> x >> y >> z >> r >> esz;
		MeshKernel::iGameVertex vex(x, y, z);
		auto vh = new_mesh->AddVertex(vex);
		new_mesh->radius[vh] = r;
		vector<int> neighbor;
		for (int j = 0; j < esz; ++j) {
			strs >> vid;
			neighbor.emplace_back(vid);
		}
		neighbors.emplace_back(neighbor);
	}

	for (int i = 0; i < vsz; ++i) {
		VH vh1(i);
		for (auto& j : neighbors[i]) {
			if (j > i) continue;
			VH vh2(j);
			new_mesh->AddEdge(vh1, vh2);
		}
	}

	inputfile.close();
	*this = *new_mesh;

	std::cout << "Read Skeletal Mesh success. vertices = " << this->vsize() <<
		", edges = " << this->esize() << ", faces = " << this->fsize() << std::endl;

	return true;

}

bool SkeletalMesh::writeSkelFile(std::string filename) {

	std::ofstream off(filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
		off.close();
		return false;
	}
	updateAllHandles();
	off << "ID Cx Cy Cz RADIUS #NEIGHBORS NEIGHBORS_LIST" << std::endl;
	off << vsize() << std::endl;
	for (int i = 0; i < vsize(); ++i) {
		VH vh(i);
		auto& v = vertices_[vh];
		auto adjvhs = NeighborVh(vh);
		off << i << " " << v.x() << " " << v.y() << " " << v.z() << " " << radius[vh] << " " << adjvhs.size();
		for (auto& adjvh : adjvhs) {
			off << " " << adjvh;
		}
		off << std::endl;
	}
	off.close();
	std::cout << "Writing file success. #V = " << vsize() << ", #E = " << esize() << ", #F = " << fsize() << std::endl;
	return true;
}
