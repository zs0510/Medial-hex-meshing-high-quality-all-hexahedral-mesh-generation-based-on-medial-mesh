#include "InpIO.h"

MeshKernel::VolumeMesh InpIO::get_volume_mesh_from_inp_file(std::string filename) {

	MeshKernel::VolumeMesh mesh;

	std::ifstream inputfile(filename, std::ios::in);
	std::string line;
	std::vector<Vex> vertices;
	std::vector<std::vector<VH>> cells;
	std::unordered_map<int, VH> nid2vh;

	enum ReadMode {// BOUNDLINE ��δʵ�� 
		NODE, ELEMENT, INFO
	};

	ReadMode read_mode = ReadMode::INFO;

	while (inputfile) {

		line.clear();
		getline(inputfile, line);// ����ÿһ��

		if (line.empty()) continue;
		else if (line[0] == '*') {// ע��
			if (line[1] == 'N') read_mode = ReadMode::NODE;
			else if (line[1] == 'E') read_mode = ReadMode::ELEMENT;

		} else {

			std::stringstream linestream;
			linestream.str(line);

			if (read_mode == ReadMode::NODE) {
				// ����ڵ�
				int nid;
				double x, y, z;
				char c;
				linestream >> nid >> c >> x >> c >> y >> c >> z;
				nid2vh[nid] = VH(vertices.size());
				Vex v(x, y, z);
				vertices.emplace_back(v);

			} else if (read_mode == ReadMode::ELEMENT){
				// ���뵥Ԫ
				std::vector<int> nodes;
				std::vector<VH> vhs;
				int nid;
				char c;
				linestream >> nid >> c;// ��Ԫ ID
				while (linestream >> nid) {
					nodes.emplace_back(nid);
					if (nodes.size() != 8) {
						linestream >> c;
					}
				}
				if (nodes.size() == 7) {// ���һ���ڵ� ID ����������һ����
					line.clear();
					getline(inputfile, line);// ������һ��
					linestream.clear();
					linestream.str(line);
					linestream >> nid;
					nodes.emplace_back(nid);
				}

				bool is_all_valid = true;
				for (auto& i : nodes) {
					if (nid2vh.count(i)) {
						vhs.emplace_back(nid2vh[i]);

					} else {
						is_all_valid = false;
						break;

					}
				}

				if (is_all_valid && vhs.size() == 8) {
					cells.emplace_back(vhs);

				} else if (vhs.size() != 8) {
					std::cout << "[Inp File IO]: Error vhs size!!!\n";

				} else if (!is_all_valid) {
					std::cout << "[Inp File IO]: Not all nid has been initialized!!!\n";

				}

			}

		}

	}

	std::cout << "[Inp File IO]: Read " << vertices.size() << " vertices, " << cells.size() << " cells.\n";

	for (auto& v : vertices) {
		mesh.AddVertex(v);
	}
	for (auto& cell : cells) {
		mesh.AddCell(cell);
	}

	return mesh;

}