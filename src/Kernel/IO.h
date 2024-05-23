#pragma once
#include "Mesh.h"
#include <fstream>
#include <string>
#include <sstream>
#include <array>

using namespace std;

namespace MeshKernel {
	class IO {
	public:
		// 旧的，暂时留存
		VolumeMesh iGameReadVolumeFile(const std::string&);// 读六面体网格
		void iGameWriteVolumeFile(const VolumeMesh& _mesh, const std::string&);// 写六面体网格

		TetMesh iGameReadTetMeshFile(const std::string&);// 读四面体网格
		TetMesh ReadTetMeshFromVTK(const std::string&);

		SurfaceMesh ReadObjFile(const std::string& _InputFile);
		SurfaceMesh ReadOffFile(const std::string& _InputFile);
		//SurfaceMesh ReadStlFile(const std::string& _InputFile, int& sides_num);
		//TetMesh ReadTetMeshFile(const std::string& _InputFile);
		bool WriteObjFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
		bool WriteOffFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
		//bool WriteStlFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
		//bool WriteMeshFile(const VolumeMesh& _mesh, const std::string& _OutputFile);
		//bool WriteVtkFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
		std::string WriteOffString(const SurfaceMesh& _mesh);
		TetMesh ReadMeshFileFromStr(const std::string& data);
		std::vector<std::string> SplitFileName(const std::string& fileName);

		VolumeMesh ReadHexMesh_HexaLab(const std::string);
		VolumeMesh ReadVtkFile_Volume(const std::string& _InputFile);

	private:
		void ReOrderiGameVertexHandle(const SurfaceMesh& _mesh);
		std::vector<iGameVertexHandle> reorderedvh_;                        // 重排顶点
		std::unordered_map<iGameVertexHandle, std::size_t> newvh_;          // 新的顶点handle
	};
}

MeshKernel::TetMesh MeshKernel::IO::ReadTetMeshFromVTK(const std::string& filename) {

	std::ifstream inputfile(filename, std::ios::in);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> cells;
	std::string line;
	std::cout << "Is Reading : " << filename << " File." << std::endl;
	int points_count, cells_count;
	double x, y, z;
	int vcnt, vid;

	do {
		line.clear();
		getline(inputfile, line);
		if (line[0] == 'P') {
			std::stringstream linestream;
			linestream.str(line);
			linestream >> line >> points_count;
			std::cout << "Points size = " << points_count << std::endl;
			break;
		}
	} while (inputfile);

	vertices.reserve(points_count);

	for (int i = 0; i < points_count; ++i) {
		line.clear();
		getline(inputfile, line);
		std::stringstream linestream;
		linestream.str(line);
		linestream >> x >> y >> z;
		Vex v(x, y, z);
		vertices.emplace_back(v);
	}

	do {
		line.clear();
		getline(inputfile, line);
		if (line[0] == 'C') {
			std::stringstream linestream;
			linestream.str(line);
			linestream >> line >> cells_count;
			std::cout << "Cells size = " << cells_count << std::endl;
			break;
		}
	} while (inputfile);

	cells.reserve(cells_count);

	for (int i = 0; i < cells_count; ++i) {
		line.clear();
		getline(inputfile, line);
		std::stringstream linestream;
		linestream.str(line);
		linestream >> vcnt;
		vector<VH> vhs;
		for (int j = 0; j < vcnt; ++j) {
			linestream >> vid;
			vhs.push_back(VH(vid));
		}
		cells.emplace_back(vhs);
	}

	inputfile.close();

	return TetMesh(vertices, cells);

}

MeshKernel::VolumeMesh MeshKernel::IO::ReadHexMesh_HexaLab(const std::string path) {

	string header;

	std::vector<Vex> vertices;
	std::vector<std::vector<VH>> indices;// hex

	ifstream stream(path, ifstream::in | ifstream::binary);

	int precision;
	int dimension;

	MeshKernel::VolumeMesh mesh;

	while (stream.good()) {
		// Read a line
		stream >> header;
		if (header.compare("MeshVersionFormatted") == 0) {
			stream >> precision;
			//HL_ASSERT_LOG(stream >> precision, "ERROR: malformed mesh file. Unexpected value after %s tag.\n", header.c_str());
		} else if (header.compare("Dimension") == 0) {
			stream >> dimension;
			//HL_ASSERT_LOG(stream >> dimension, "ERROR: malformed mesh file. Unexpected value after %s tag.\n", header.c_str());
		} else if (header.compare("Vertices") == 0) {
			int vertices_count;
			/*HL_ASSERT_LOG(stream >> vertices_count, "ERROR: malformed mesh file. Unexpected value after %s tag.\n", header.c_str());
			HL_LOG("[Loader] Reading %d vertices...\n", vertices_count);*/
			stream >> vertices_count;
			printf("[Read HexMesh] Reading %d vertices...\n", vertices_count);
			vertices.reserve(vertices_count);
			for (int i = 0; i < vertices_count; ++i) {
				Vex v;
				float x;
				//HL_ASSERT_LOG(stream >> v.x() >> v.y() >> v.z() >> x, "ERROR: malformed mesh file. Unexpected vertex data format at vert %i.\n", i);
				stream >> v.x() >> v.y() >> v.z() >> x;
				vertices.push_back(v);
			}
		} else if (header.compare("Quadrilaterals") == 0 || header.compare("Quads") == 0) {
			int quads_count;
			stream >> quads_count;
			for (int i = 0; i < quads_count; ++i) {
				int idx[4];
				int x;
				stream >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> x;
			}
		} else if (header.compare("Hexahedra") == 0) {
			int hexas_count;
			/*HL_ASSERT_LOG(stream >> hexas_count, "ERROR: malformed mesh file. Unexpected tag after hexahedras tag.\n");
			HL_LOG("[Loader] Reading %d hexas...\n", hexas_count);*/
			stream >> hexas_count;
			printf("[Read HexMesh] Reading %d hexas...\n", hexas_count);
			indices.reserve(hexas_count);
			for (int h = 0; h < hexas_count; ++h) {
				int idx[8];
				std::vector<VH> vhs;
				int x;
				// irrational to rational vertex ordering!
				//stream >> idx[0] >> idx[1] >> idx[3] >> idx[2] >> idx[4] >> idx[5] >> idx[7] >> idx[6] >> x;// HexaLab
				stream >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> idx[4] >> idx[5] >> idx[6] >> idx[7] >> x;// iGame
				for (int i = 0; i < 8; ++i) {
					vhs.push_back(VH(idx[i] - 1));
				}
				indices.push_back(vhs);

			}
		} else if (header.compare("Triangles") == 0) {
			int tri_count;
			stream >> tri_count;
			for (int i = 0; i < tri_count; ++i) {
				int idx[3];
				int x;
				stream >> idx[0] >> idx[1] >> idx[2] >> x;
			}
		} else if (header.compare("Edges") == 0) {
			int edge_count;
			stream >> edge_count;
			for (int i = 0; i < edge_count; ++i) {
				int idx[2];
				int x;
				stream >> idx[0] >> idx[1] >> x;
			}
		} else if (header.compare("Corners") == 0) {
			int corner_count;
			stream >> corner_count;
			for (int i = 0; i < corner_count; ++i) {
				int c;
				stream >> c;
			}
		} else if (header.compare("End") == 0) {
			break;
		} else if (header[0] == '#') {
			stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		} else {
			printf("[Loader] Unexpected header \"%s\"\n ", header.c_str());
			printf("[Loader] ERROR: malformed mesh file.Unexpected header tag.\n");
			//return mesh;
		}
	}

	// Make sure at least vertex and hexa index data was read
	if (vertices.size() == 0) {
		printf("ERROR: mesh does not contain any vertex!\n");
		return mesh;
	} else if (indices.size() == 0) {
		printf("ERROR: mesh does not contain any thetra or hexa!\n");
		return mesh;
	}
	
	mesh.InitMesh(vertices, indices);
	
	return mesh;

}

MeshKernel::SurfaceMesh MeshKernel::IO::ReadObjFile(const std::string& _InputFile) {
	std::ifstream inputfile(_InputFile, std::ios::in);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> faces;
	std::vector<std::vector<double>> normals;
	std::vector<std::vector<double>> uvs;
	std::unordered_map<int, int> V2N;// vertex to normal
	std::unordered_map<int, int> V2T;// vertex to uv
	std::string line;

	std::cout << "Reading " << _InputFile << " File" << std::endl;
	// std::cout << inputfile.good() << std::endl;
	while (inputfile) {
		line.clear();
		getline(inputfile, line);
		if (line[0] == '#') {
			continue;// 注释
		}
		std::stringstream linestream;
		linestream.str(line);

		std::string flag;
		linestream >> flag;
		if (flag == "v") {
			double x, y, z;
			linestream >> x >> y >> z;
			vertices.push_back(iGameVertex(x, y, z));
		} else if (flag == "f") {
			// f 1575/1514/1569 1581/1520/1575 1576/1515/1570
			// f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
			std::vector<std::string> vex;
			std::string tmp;
			while (linestream >> tmp) vex.push_back(tmp);
			auto n = vex.size();
			//std::cout << n << std::endl;
			// 创建面
			std::vector<iGameVertexHandle> face(n);
			for (size_t i = 0; i < n; i++) {
				size_t idx = 0;
				while (idx < vex[i].length() && std::isdigit(vex[i][idx])) idx++;
				int vh = std::stoi(vex[i].substr(0, idx)) - 1;// 注意 obj 文件 v 从1开始
				face[i] = (iGameVertexHandle)(vh);
				//if (idx != vex[i].length()) {// 指定了纹理坐标
				//	// 注意：纹理坐标存在复用！法向量存在复用！顶点、纹理坐标与法向量个数均可互相不等
				//	size_t beg = idx;
				//	if (beg < vex[i].length() && !std::isdigit(vex[i][beg])) beg++;
				//	size_t end = beg;
				//	if (end < vex[i].length() && std::isdigit(vex[i][end])) end++;
				//	if (!V2T.count(vh)) {
				//		int uv_idx = std::stoi(vex[i].substr(beg, end - beg));// 注意 obj 文件 vt 从0开始
				//		V2T[vh] = uv_idx;
				//	}
				//	if (end != vex[i].length() && !V2N.count(vh)) {// 指定了法向量
				//		beg = end;
				//		if (beg < vex[i].length() && !std::isdigit(vex[i][beg])) beg++;
				//		end = beg;
				//		if (end < vex[i].length() && std::isdigit(vex[i][end])) end++;
				//		int n_idx = std::stoi(vex[i].substr(beg, end - beg));// 注意 obj 文件 vn 从0开始
				//		V2N[vh] = n_idx;
				//	}
				//}
			}
			faces.push_back(face);
		} else if (flag == "vt") {
			double u, v;
			linestream >> u >> v;
			uvs.push_back({ u, v });
		} else if (flag == "vn") {
			double x, y, z;
			linestream >> x >> y >> z;
			normals.push_back({ x, y, z });
		}
	}
	//printf("read file success, fcnt: %d, vcnt: %d, vtcnt: %d, vncnt: %d\n", faces.size(), vertices.size(), uvs.size(), normals.size());
	if (!normals.empty()) {
		int ncnt = normals.size();
		for (int i = 0; i < vertices.size(); ++i) {
			int nidx = V2N[i];
			//if (nidx < 0 || nidx >= ncnt) printf("error: nidx = %d\n", nidx);// debug 用
			assert(nidx >= 0 && nidx < ncnt);
			vertices[i].setNormal(normals[nidx]);
		}
	}
	
	auto mesh = SurfaceMesh(vertices, faces);
	inputfile.close();
	return mesh;
}

MeshKernel::VolumeMesh MeshKernel::IO::ReadVtkFile_Volume(const std::string& _InputFile) {
	/*
	ZQF`s Order
		   v
	3----------2
	|\     ^   |\
	| \    |   | \
	|  \   |   |  \
	|   7------+---6
	|   |  +-- |-- | -> u : (+) is the barycenter of the cube
	0---+---\--1   |
	 \  |    \  \  |
	  \ |     \  \ |
	   \|      w  \|
		4----------5
	face = { {0,3,2,1},{0,4,7,3},{0,1,5,4},{4,5,6,7},{1,2,6,5},{2,3,7,6} }
*/
/*
	VTK File Format Order
		   v
	7----------6
	|\     ^   |\
	| \    |   | \
	|  \   |   |  \
	|   4------+---5
	|   |  +-- |-- | -> u : (+) is the barycenter of the cube
	3---+---\--2   |
	 \  |    \  \  |
	  \ |     \  \ |
	   \|      w  \|
		0----------1
*/

	std::ifstream hexfile(_InputFile.c_str());

	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> cells;

	if (hexfile.is_open())     {
		std::cout << "Reading " << _InputFile << " File" << std::endl;

		std::string str;
		try         {
			do             {
				hexfile >> str;
			} while (str != "DATASET");
			hexfile >> str;
			if (str != "UNSTRUCTURED_GRID")             {
				//should stop reading, but not implement here
				1 + 1 == 2;
			}
			do             {
				hexfile >> str;
			} while (str != "POINTS");

			size_t nv, ne, num_ne_type;

			hexfile >> nv >> str;

			vertices.resize(nv);
			for (size_t i = 0; i < nv; i++)             {
				hexfile >> vertices[i].x() >> vertices[i].y() >> vertices[i].z();
			}

			hexfile >> str >> ne >> num_ne_type;

			cells.reserve(ne);
			std::vector<iGameVertexHandle> vhs_per_cell;
			vhs_per_cell.resize(8);
			int index[8] = { 4, 5, 1, 0, 7, 6, 2, 3 };
			int tmp = 0;
			for (size_t i = 0; i < ne; i++)             {
				//处理掉第一个numPoints Flag
				hexfile >> tmp; // 8

				for (int j = 0; j < 8; j++)                 {
					hexfile >> tmp;
					vhs_per_cell[index[j]] = iGameVertexHandle(tmp);
				}
				cells.push_back(vhs_per_cell);
			}
		}         catch (...)         {
			hexfile.close();
		}
		hexfile.close();
		std::cout << "Reading " << " Finished" << std::endl;
	}
	auto mesh = MeshKernel::VolumeMesh(vertices, cells);

	return mesh;
}


std::string MeshKernel::IO::WriteOffString(const SurfaceMesh& _mesh) {
	std::stringstream ret;
	ReOrderiGameVertexHandle(_mesh);
	ret << "OFF" << std::endl;
	ret << _mesh.allvertices().size() << " " << _mesh.allfaces().size() << " " << _mesh.alledges().size() << std::endl;
	for (iGameVertexHandle vh : reorderedvh_) {
		iGameVertex v(_mesh.vertices(vh));
		ret << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}
	auto allf = _mesh.allfaces();
	for (auto f : allf) {
		ret << f.second.size();
		for (int i = 0; i < f.second.size(); ++i) {
			ret << " " << newvh_[f.second.vh(i)];
		}
		ret << std::endl;
	}
	return std::string(ret.str());
}
MeshKernel::SurfaceMesh MeshKernel::IO::ReadOffFile(const std::string& _InputFile) {
	std::ifstream inputfile(_InputFile, std::ios::in);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> faces;
	std::vector<std::vector<double>> normals;
	std::vector<std::vector<double>> uvs;
	std::unordered_map<int, int> V2N;// vertex to normal
	std::unordered_map<int, int> V2T;// vertex to uv
	std::string line;
	int v_size, f_size, e_size;

	std::cout << "Reading " << _InputFile << " File" << std::endl;
	line.clear();
	getline(inputfile, line);
	std::stringstream linestream;

	if (line == "OFF") {
		line.clear();
		getline(inputfile, line);
		linestream.str(line);
		linestream >> v_size >> f_size >> e_size;
	}
	for (int i = 0; i < v_size; i++) {
		line.clear();
		getline(inputfile, line);
		std::stringstream linestream;
		linestream.str(line);
		double x, y, z;
		linestream >> x >> y >> z;
		vertices.push_back(iGameVertex(x, y, z));
		//std::cout << x << " " << y << " " << z << std::endl;
	}
	for (int i = 0; i < f_size; i++) {
		line.clear();
		getline(inputfile, line);
		std::stringstream linestream;
		linestream.str(line);
		int v;
		linestream >> v;
		//std::cout << v << " " << x << " " << y << " " << z << std::endl;
		std::vector<iGameVertexHandle> face(v);
		for (int i = 0; i < v; i++) {
			int temp;
			linestream >> temp;
			face[i] = (iGameVertexHandle)(temp);
		}
		faces.push_back(face);
	}
	printf("read file success, fcnt: %d, vcnt: %d, vtcnt: %d, vncnt: %d\n", faces.size(), vertices.size(), uvs.size(), normals.size());
	if (!normals.empty()) {
		int ncnt = normals.size();
		for (int i = 0; i < vertices.size(); ++i) {
			int nidx = V2N[i];
			//if (nidx < 0 || nidx >= ncnt) printf("error: nidx = %d\n", nidx);// debug 用
			assert(nidx >= 0 && nidx < ncnt);
			vertices[i].setNormal(normals[nidx]);
		}
	}
	
	auto mesh = SurfaceMesh(vertices, faces);
	inputfile.close();
	return mesh;
}


//MeshKernel::SurfaceMesh MeshKernel::IO::ReadStlFile(const std::string& _InputFile, int& sides_num) {
//
//	char const* file_path = _InputFile.c_str();
//
//	size_t file_size;
//
//	FILE* file = fopen(file_path, "rb");
//	if (!file) {
//		fprintf(stderr, "file %s does not exist\n", file_path);
//		exit(EXIT_FAILURE);
//	}
//
//	fseek(file, 0, SEEK_END);
//	file_size = ftell(file);
//	rewind(file);
//
//	char* file_contents = (char*)malloc(file_size);
//	if (fread(file_contents, 1, file_size, file) != file_size) {
//		fprintf(stderr, "I/O error while reading %s\n", file_path);
//		exit(EXIT_FAILURE);
//	}
//	fclose(file);
//
//	size_t buf_size;
//	if (one_stl_buf_size(&buf_size, file_contents, file_size, ONE_STL_NVVV)) {
//		fprintf(stderr, "file \"%s\" is not valid stl data\n", file_path);
//		exit(EXIT_FAILURE);
//	}
//
//	float* buf = (float*)malloc(buf_size);
//	size_t trig_count = one_stl_parse(buf, file_contents, ONE_STL_NVVV);
//
//	printf("triangle count: %ld\n", trig_count);
//
//	typedef struct {
//		float normal[3];
//		float v0[3];
//		float v1[3];
//		float v2[3];
//	} trig_nvvv_s;
//	trig_nvvv_s* triangles = (trig_nvvv_s*)buf;
//	float EPS = 0;
//	//for (size_t i = 0; i < trig_count; ++i){
//	//    EPS=std::max(EPS,(triangles + i)->v0[0]);
//	//    EPS=std::max(EPS,(triangles + i)->v0[1]);
//	//    EPS=std::max(EPS,(triangles + i)->v0[2]);
//	//}
//	//EPS/=100000;
//	//EPS = 0;
//
//	struct PointEps {
//		PointEps() {}
//		PointEps(float a, float b, float c, float EPS) {
//			v0[0] = a; v0[1] = b;
//			v0[2] = c;
//			this->EPS = EPS;
//		}
//		float v0[3];
//		float EPS;
//		bool operator < (const PointEps& other) const {
//			for (int i = 0; i < 3; i++) {
//				if (abs(v0[i] - other.v0[i]) > EPS) {
//					return  v0[i] + EPS < other.v0[i];
//				}
//			}
//			return false;
//		}
//	};
//
//	std::map<PointEps, int>mp;
//	std::vector<iGameVertex> vertices;
//	std::vector<std::vector<iGameVertexHandle>> faces;
//	std::function<int(float, float, float)> getHandleId = [&](float x, float y, float z) {
//		if (mp.count(PointEps(x, y, z, EPS))) {
//			return mp[PointEps(x, y, z, EPS)];
//		}
//		int cnt = mp.size();
//		mp[PointEps(x, y, z, EPS)] = cnt;
//		vertices.push_back(MeshKernel::iGameVertex{ x,y,z });
//		return cnt;
//	};
//	std::vector<std::vector<float> >normal;
//	for (size_t i = 0; i < trig_count; ++i) {
//		faces.push_back({ (iGameVertexHandle)getHandleId((triangles + i)->v0[0],(triangles + i)->v0[1],(triangles + i)->v0[2]),
//						(iGameVertexHandle)getHandleId((triangles + i)->v1[0],(triangles + i)->v1[1],(triangles + i)->v1[2]),
//						(iGameVertexHandle)getHandleId((triangles + i)->v2[0],(triangles + i)->v2[1],(triangles + i)->v2[2])
//			});
//		normal.push_back(std::vector<float>{ (triangles + i)->normal[0], (triangles + i)->normal[1], (triangles + i)->normal[2] });
//	}
//	std::cout << vertices.size() << " " << faces.size() << std::endl;
//	auto mesh = SurfaceMesh(vertices, faces);
//	std::cout << trig_count << " " << faces.size() << std::endl;
//	//for(int i=0;i<faces.size();i++){
//	//    mesh.faces(iGameFcaeHandle(i)).setNormal(normal[i][0], normal[i][1], normal[i][2]);
//	//}   
//	free(file_contents);
//	free(buf);
//	return mesh;
//}


bool MeshKernel::IO::WriteObjFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
	std::ofstream outputfile(_OutputFile, std::ios::out);
	ReOrderiGameVertexHandle(_mesh);
	for (iGameVertexHandle vh : reorderedvh_) {
		iGameVertex v(_mesh.vertices(vh));
		outputfile << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}
	auto allf = _mesh.allfaces();
	for (auto f : allf) {
		outputfile << "f";
		for (int i = 0; i < f.second.size(); ++i) {
			outputfile << " " << newvh_[f.second.vh(i)] + 1;
		}
		outputfile << std::endl;
	}
	outputfile.close();
	return true;
}

bool MeshKernel::IO::WriteOffFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
	std::ofstream outputfile(_OutputFile, std::ios::out);
	ReOrderiGameVertexHandle(_mesh);
	outputfile << "OFF" << std::endl;
	outputfile << _mesh.allvertices().size() << " " << _mesh.allfaces().size() << " " << 0 << std::endl;
	for (iGameVertexHandle vh : reorderedvh_) {
		iGameVertex v(_mesh.vertices(vh));
		outputfile << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}
	auto allf = _mesh.allfaces();
	for (auto f : allf) {
		outputfile << f.second.size();
		for (int i = 0; i < f.second.size(); ++i) {
			outputfile << " " << newvh_[f.second.vh(i)];
		}
		outputfile << std::endl;
	}
	outputfile.close();
	return true;
}

//bool MeshKernel::IO::WriteStlFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
//	STLWrite c;
//	for (auto i : _mesh.allfaces()) {
//		printf("size : %d\n", i.second.getSortediGameVertexHandle().size());
//
//		if (i.second.getSortediGameVertexHandle().size() == 3) {
//			for (int j = 0; j < 3; j++) {
//				auto point = _mesh.vertices(i.second.vh(j));
//				c.add(point.x(), point.y(), point.z());
//			}
//		} 		else if (i.second.getSortediGameVertexHandle().size() == 4) {
//			for (int j = 0; j < 3; j++) {
//				auto point = _mesh.vertices(i.second.vh(j));
//				c.add(point.x(), point.y(), point.z());
//			}
//
//			auto point = _mesh.vertices(i.second.vh(0));
//			c.add(point.x(), point.y(), point.z());
//			point = _mesh.vertices(i.second.vh(2));
//			c.add(point.x(), point.y(), point.z());
//			point = _mesh.vertices(i.second.vh(3));
//			c.add(point.x(), point.y(), point.z());
//
//		}
//	}
//	c.write(_OutputFile);
//	return true;
//}

std::vector<std::string> MeshKernel::IO::SplitFileName(const std::string& fileName)
{
	// JFR DO NOT CHANGE TO std::vector<std::string> s(3), it segfaults while
	// destructor si called
	std::vector<std::string> s;
	s.resize(3);
	if (fileName.size()) {
		// returns [path, baseName, extension]
		int idot = (int)fileName.find_last_of('.');
		int islash = (int)fileName.find_last_of("/\\");
		if (idot == (int)std::string::npos) idot = -1;
		if (islash == (int)std::string::npos) islash = -1;
		if (idot > 0) s[2] = fileName.substr(idot);
		if (islash > 0) s[0] = fileName.substr(0, islash + 1);
		s[1] =
			fileName.substr(s[0].size(), fileName.size() - s[0].size() - s[2].size());
	}
	return s;
}

//bool MeshKernel::IO::WriteVtkFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
//	// 还不支持约束线的输出
//
//	std::vector<std::array<double, 3>>verts;
//	std::vector<std::vector<uint32_t>> faces;
//	std::vector<std::array<int, 2>> constrains;
//
//	for (int i = 0; i < _mesh.allvertices().size(); i++) {
//		std::array<double, 3> vert;
//		vert[0] = _mesh.vertices(MeshKernel::iGameVertexHandle(i)).x();
//		vert[1] = _mesh.vertices(MeshKernel::iGameVertexHandle(i)).y();
//		vert[2] = _mesh.vertices(MeshKernel::iGameVertexHandle(i)).z();
//		verts.push_back(vert);
//	}
//
//	for (auto fp : _mesh.allfaces()) {
//		auto fh = fp.first;
//		auto face = fp.second;
//		std::vector<uint32_t> oneface;
//		for (int i = 0; i < face.getVertexHandle().size(); i++) {
//			oneface.push_back(face.getVertexHandle()[i]);
//		}
//		faces.push_back(oneface);
//	}
//
//	int numVertices = verts.size();
//
//	FILE* fp = fopen(_OutputFile.c_str(), "w");
//
//	std::vector<std::string> s = SplitFileName(_OutputFile);
//	s.resize(3);
//
//	if (!fp) {
//		fprintf(stdout, "Unable to open file '%s'\n", _OutputFile.c_str());
//		return false;
//	}
//
//	fprintf(fp, "# vtk DataFile Version 2.0\n");
//	fprintf(fp, "%s, Created by labb \n", s[1].c_str());
//	fprintf(fp, "ASCII\n");
//	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
//
//	// write mesh vertices
//	fprintf(fp, "POINTS %d double\n", numVertices);
//	for (std::size_t i = 0; i < verts.size(); i++)
//		fprintf(fp, "%.16g %.16g %.16g\n", verts[i][0], verts[i][1], verts[i][2]);
//	fprintf(fp, "\n");
//
//	int numElements = faces.size();
//	//int totalNumInt = numElements * 5;
//	int totalNumInt = 0;
//	for (std::size_t i = 0; i < faces.size(); i++) {
//		if (faces[i].size() == 3) totalNumInt += 4;
//		else if (faces[i].size() == 4) totalNumInt += 5;
//	}
//
//	for (std::size_t i = 0; i < constrains.size(); ++i) {
//		++numElements;
//		totalNumInt += 3;
//	}
//
//	// print vertex indices in ascii or binary
//	fprintf(fp, "CELLS %d %d\n", numElements, totalNumInt);
//	for (std::size_t i = 0; i < constrains.size(); ++i) {
//		fprintf(fp, "%d", 2);
//		fprintf(fp, " %ld", constrains[i][0]);
//		fprintf(fp, " %ld", constrains[i][1]);
//		fprintf(fp, "\n");
//	}
//
//	for (std::size_t i = 0; i < faces.size(); i++) {
//		fprintf(fp, "%d", faces[i].size());
//		for (int j = 0; j < faces[i].size(); j++)
//			fprintf(fp, " %ld", faces[i][j]);
//		fprintf(fp, "\n");
//	}
//	fprintf(fp, "\n");
//
//	// print element types in ascii or binary
//	fprintf(fp, "CELL_TYPES %d\n", numElements);
//	for (std::size_t i = 0; i < numElements - faces.size(); i++) {
//		fprintf(fp, "%d\n", 3);
//	}
//	for (std::size_t i = 0; i < faces.size(); i++) {
//		if (faces[i].size() == 3) fprintf(fp, "%d\n", 5);
//		else if (faces[i].size() == 4) fprintf(fp, "%d\n", 9);
//	}
//	fclose(fp);
//	return true;
//}

void MeshKernel::IO::ReOrderiGameVertexHandle(const SurfaceMesh& _mesh) {
	auto allv = _mesh.allvertices();
	int idx = 0;
	for (auto v : allv) {
		reorderedvh_.push_back(v.first);
		newvh_[v.first] = idx++;
	}
}


//MeshKernel::TetMesh MeshKernel::IO::ReadMeshFileFromStr(const std::string& data) {
//	std::stringstream ss(data);
//	std::vector<iGameVertex> vertices;
//	std::vector<std::vector<iGameVertexHandle>> surface_faces;
//	std::vector<std::vector<std::vector<iGameVertexHandle>>> cells;
//	std::vector<std::vector<iGameVertexHandle>> faces;
//	std::string str;
//	enum State
//	{
//		USELESS = 1, VERTEX, FACE, TET, EDGE
//	}state = USELESS;
//	while (std::getline(ss, str)) {
//		int len = str.length();
//		std::vector<std::string>info;
//		std::string s;
//		for (int i = 0; i < len; i++) {
//			if (str[i] == '#')break;
//			if (str[i] == ' ') {
//				if (s.length() > 0)
//					info.push_back(s);
//				s = "";
//			}             else {
//				s.push_back(str[i]);
//			}
//		}
//		if (s.length() > 0)
//			info.push_back(s);
//		if (info.size() == 0)continue;
//		if (info[0] == "MeshVersionFormatted") {
//			state = USELESS;
//		}         else if (info[0] == "Dimension") {
//			std::getline(ss, str);
//			state = USELESS;
//		}         else if (info[0] == "Vertices") {
//			std::getline(ss, str);
//			state = VERTEX;
//		}         else if (info[0] == "Tetrahedra") {
//			std::getline(ss, str);
//			state = TET;
//		}         else if (info[0] == "Triangles") {
//			std::getline(ss, str);
//			state = FACE;
//		}         else if (info[0] == "iGameEdges") {
//			std::getline(ss, str);
//			state = EDGE;
//		}         else if (info[0] == "End") {
//			state = USELESS;
//		}         else {
//			if (state == USELESS) {
//				continue;
//			}             else if (state == VERTEX) {
//				vertices.push_back(MeshKernel::iGameVertex(std::stod(info[0])
//					, std::stod(info[1]), std::stod(info[2])));
//			}             else if (state == TET) {
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//				cells.push_back(faces);
//				faces.clear();
//			}             else if (state == FACE) {
//				surface_faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//			}
//		}
//	}
//	return TetMesh(vertices, cells, surface_faces);
//}
//MeshKernel::TetMesh MeshKernel::IO::ReadTetMeshFile(const std::string& _InputFile) {
//	FILE* fp = fopen(_InputFile.c_str(), "r");
//	char str[100];
//	enum State
//	{
//		USELESS = 1, VERTEX, FACE, TET, EDGE
//	}state = USELESS;
//	std::vector<iGameVertex> vertices;
//	std::vector<std::vector<iGameVertexHandle>> surface_faces;
//	std::vector<std::vector<std::vector<iGameVertexHandle>>> cells;
//	std::vector<std::vector<iGameVertexHandle>> faces;
//	while (fscanf(fp, "%[^\n]\n", str) != EOF) {
//		int len = strlen(str);
//		std::vector<std::string>info;
//		std::string s;
//		for (int i = 0; i < len; i++) {
//			if (str[i] == '#')break;
//			if (str[i] == ' ') {
//				if (s.length() > 0)
//					info.push_back(s);
//				s = "";
//			} 			else {
//				s.push_back(str[i]);
//			}
//		}
//		if (s.length() > 0)
//			info.push_back(s);
//		if (info.size() == 0)continue;
//		if (info[0] == "MeshVersionFormatted") {
//			state = USELESS;
//		} 		else if (info[0] == "Dimension") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = USELESS;
//		} 		else if (info[0] == "Vertices") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = VERTEX;
//		} 		else if (info[0] == "Tetrahedra") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = TET;
//		} 		else if (info[0] == "Triangles") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = FACE;
//		} 		else if (info[0] == "iGameEdges") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = EDGE;
//		} 		else if (info[0] == "End") {
//			state = USELESS;
//		} 		else {
//			if (state == USELESS) {
//				continue;
//			} 			else if (state == VERTEX) {
//				vertices.push_back(MeshKernel::iGameVertex(std::stod(info[0])
//					, std::stod(info[1]), std::stod(info[2])));
//			} 			else if (state == TET) {
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//				cells.push_back(faces);
//				faces.clear();
//			} 			else if (state == FACE) {
//				surface_faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//			}
//		}
//	}
//	std::cout << "read file success" << std::endl;
//	return TetMesh(vertices, cells, surface_faces);
//}

void MeshKernel::IO::iGameWriteVolumeFile(const VolumeMesh& _mesh, const std::string& _filename) {


	std::ofstream off(_filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
		off.close();
		return;
	}

	// write header
	off << "MeshVersionFormatted 1" << std::endl;
	off << "Dimension 3" << std::endl;
	

	// write vertices
	uint64_t n_vertices(_mesh.vsize());
	off << "Vertices" << std::endl;
	off << n_vertices << std::endl;

	for (uint64_t v_it = 0; v_it < n_vertices; ++v_it) {
		auto& v = _mesh.vertices((MeshKernel::iGameVertexHandle)v_it);
		off <<  v.x() << " "
			<<  v.y() << " "
			<<  v.z() << " "
			<< "-1" << std::endl;// reference
	}
	
	// write hexahedra
	uint64_t n_cells(_mesh.csize());
	off << "Hexahedra" << std::endl;
	off << n_cells << std::endl;

	for (auto& cp : _mesh.allcells()) {
		auto vhs = cp.second.getVertexHandle();
		for (auto& vh : vhs) {
			off << (int)vh + 1 << " ";
		}
		off << "1" << std::endl;
	}

	off << "End" << std::endl;
	off.close();

}


MeshKernel::VolumeMesh MeshKernel::IO::iGameReadVolumeFile(const std::string& _InputFile) {

	std::ifstream inputfile(_InputFile, std::ios::in);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> cells;
	std::string line;
	std::cout << "Is Reading : " << _InputFile << " File." << std::endl;

	std::vector<int> invalid_vhs;// 记录无效的顶点 id

	while (inputfile) {
		line.clear();
		getline(inputfile, line);// 读入每一行
		//std::cout << line << "\n";
		if (line[0] == '#') continue;// 注释 的情况

		else if (line[0] == 'V') {
			// 开始读点
			line.clear();
			getline(inputfile, line);// 读入下一行
			int num = stoi(line);// 点的个数
			std::cout << "该模型文件中点的个数为 : " << num << std::endl;
			for (int i = 0; i < num; ++i) {
				// 依次读入每一个点
				line.clear();
				getline(inputfile, line);
				std::stringstream linestream;
				linestream.str(line);
				double x, y, z;
				int tag = -1;
				linestream >> x >> y >> z >> tag;
				iGameVertex vv(x, y, z);
				//printf("v: %.4f, %.4f, %.4f\n", x, y, z);
				vertices.push_back(vv);
				if (tag == 0)
					invalid_vhs.push_back(vertices.size() - 1);
			}
		}

		else if (line[0] == 'H') {
			// 开始读体
			line.clear();
			getline(inputfile, line);
			std::stringstream linestream;
			int num = stoi(line);
			std::cout << "该模型文件中体的个数为 : " << num << std::endl;
			for (int i = 0; i < num; i++) {
				// 依次读入每一个体
				line.clear();
				getline(inputfile, line);
				linestream.str(line);
				std::vector<iGameVertexHandle> cellH;
				int n = 0;
				for (auto& c : line) if (c == ' ') n++;
				while (n--) {
					int t;
					linestream >> t;
					cellH.push_back(iGameVertexHandle(t - 1));// 加入点的时候需要 - 1 
				}
				cells.push_back(cellH);
			}
		}
	}
	inputfile.close();
	auto volume_mesh = VolumeMesh(vertices, cells);
	return volume_mesh;
}

MeshKernel::TetMesh MeshKernel::IO::iGameReadTetMeshFile(const std::string& _InputFile) {

	std::ifstream inputfile(_InputFile, std::ios::in);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> cells;
	std::string line;
	std::cout << "Is Reading : " << _InputFile << " File." << std::endl;

	std::vector<int> invalid_vhs;// 记录无效的顶点 id

	while (inputfile) {
		line.clear();
		getline(inputfile, line);// 读入每一行
		//std::cout << line << "\n";
		if (line[0] == '#') continue;// 注释 的情况

		else if (line[0] == 'V') {
			// 开始读点
			line.clear();
			getline(inputfile, line);// 读入下一行
			int num = stoi(line);// 点的个数
			std::cout << "该模型文件中点的个数为 : " << num << std::endl;
			for (int i = 0; i < num; ++i) {
				// 依次读入每一个点
				line.clear();
				getline(inputfile, line);
				std::stringstream linestream;
				linestream.str(line);
				double x, y, z;
				int tag = -1;
				linestream >> x >> y >> z >> tag;
				iGameVertex vv(x, y, z);
				//printf("v: %.4f, %.4f, %.4f\n", x, y, z);
				vertices.push_back(vv);
			}
		}

		else if (line[0] == 'T' && line[1] == 'e') {
			// 开始读体
			line.clear();
			getline(inputfile, line);
			std::stringstream linestream;
			int num = stoi(line);
			std::cout << "该模型文件中体的个数为 : " << num << std::endl;
			for (int i = 0; i < num; i++) {
				// 依次读入每一个体
				line.clear();
				getline(inputfile, line);
				linestream.str(line);
				std::vector<iGameVertexHandle> cellH;
				int n = 0;
				//for (auto& c : line) if (c == ' ') n++;
				n = 4;
				while (n--) {
					int t;
					linestream >> t;
					cellH.push_back(iGameVertexHandle(t - 1));// 加入点的时候需要 - 1 
				}
				cells.push_back(cellH);
			}
		}
	}
	inputfile.close();
	auto volume_mesh = TetMesh(vertices, cells);
	std::cout << "the number of vertices is " << volume_mesh.VertexSize() << std::endl;
	std::cout << "the number of cells is " << volume_mesh.CellSize() << std::endl;
	//std::cout << "该模型文件中无效点的个数为 : " << invalid_cnt << std::endl;
	return volume_mesh;

}

MeshKernel::TetMesh MeshKernel::IO::ReadMeshFileFromStr(const std::string& data) {
	std::stringstream ss(data);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> surface_faces;
	std::vector<std::vector<std::vector<iGameVertexHandle>>> cells;
	std::vector<std::vector<iGameVertexHandle>> faces;
	std::string str;
	enum State
	{
		USELESS = 1, VERTEX, FACE, TET, EDGE
	}state = USELESS;
	while (std::getline(ss, str)) {
		int len = str.length();
		std::vector<std::string>info;
		std::string s;
		for (int i = 0; i < len; i++) {
			if (str[i] == '#')break;
			if (str[i] == ' ') {
				if (s.length() > 0)
					info.push_back(s);
				s = "";
			}             else {
				s.push_back(str[i]);
			}
		}
		if (s.length() > 0)
			info.push_back(s);
		if (info.size() == 0)continue;
		if (info[0] == "MeshVersionFormatted") {
			state = USELESS;
		}         else if (info[0] == "Dimension") {
			std::getline(ss, str);
			state = USELESS;
		}         else if (info[0] == "Vertices") {
			std::getline(ss, str);
			state = VERTEX;
		}         else if (info[0] == "Tetrahedra") {
			std::getline(ss, str);
			state = TET;
		}         else if (info[0] == "Triangles") {
			std::getline(ss, str);
			state = FACE;
		}         else if (info[0] == "Edges") {
			std::getline(ss, str);
			state = EDGE;
		}         else if (info[0] == "End") {
			state = USELESS;
		}         else {
			if (state == USELESS) {
				continue;
			}             else if (state == VERTEX) {
				vertices.push_back(MeshKernel::iGameVertex(std::stod(info[0])
					, std::stod(info[1]), std::stod(info[2])));
			}             else if (state == TET) {
				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
					(iGameVertexHandle)(std::stoi(info[0]) - 1)
						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
				});
				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
					(iGameVertexHandle)(std::stoi(info[0]) - 1)
						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
				});
				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
					(iGameVertexHandle)(std::stoi(info[0]) - 1)
						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
				});
				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
					(iGameVertexHandle)(std::stoi(info[1]) - 1)
						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
				});
				cells.push_back(faces);
				faces.clear();
			}             else if (state == FACE) {
				surface_faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
					(iGameVertexHandle)(std::stoi(info[0]) - 1)
						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
				});
			}
		}
	}
	return TetMesh(vertices, cells, surface_faces);
}