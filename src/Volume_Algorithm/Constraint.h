#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include "Kernel/Mesh.h"

using namespace std;

void read_constraint_data(string filename, unordered_map<int, bool>& res) {

	std::ifstream inputfile(filename, std::ios::in);
	std::string line;
	res.clear();
	int cnt = 0, id = -1;

	while (inputfile) {
		line = "";
		id = -1;
		getline(inputfile, line);
		if (line.length() > 0) {
			id = stoi(line);
			id--;
			if (!res.count(id)) {
				res[id] = true;
				cnt++;
				//std::cout << "id: " << id << std::endl;
			}
		}
	
	}
	std::cout << "Read constraint vertices size = " << cnt << std::endl;

}

void write_constraint_data(string filename, unordered_map<int, bool>& res) {

	std::ofstream off(filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
		off.close();
		return;
	}

	for (auto& vp : res) {
		off << vp.first + 1 << std::endl;
	}
	off.close();

}

void updateAllHandles_kepp_constraint(MeshKernel::VolumeMesh& mesh, unordered_map<int, bool>& force_mp, unordered_map<int, bool>& spc_mp) {
	
	std::vector<MeshKernel::iGameVertex> newVertices;
	std::vector<std::vector<MeshKernel::iGameVertexHandle>> newCells;
	std::unordered_map<int, int> mp;// old id to new id
	int idx = 0;
	for (auto& cp : mesh.allcells()) {
		auto vhs = cp.second.getVertexHandle();
		for (auto& vh : vhs) {
			if (!mp.count(vh)) {
				mp[vh] = idx++;// vh = old_id, idx = new_id
				newVertices.push_back(mesh.vertices(vh));
			}
			vh = iGameVertexHandle(mp[vh]);// new_id
		}
		newCells.push_back(vhs);
	}

	std::unordered_map<int, bool> _force_mp = force_mp, _spc_mp = spc_mp;
	force_mp.clear();
	spc_mp.clear();
	for (auto& hp : _force_mp) {
		force_mp[mp[hp.first]] = true;
	}
	for (auto& hp : _spc_mp) {
		spc_mp[mp[hp.first]] = true;
	}

	mesh = MeshKernel::VolumeMesh(newVertices, newCells);

}
