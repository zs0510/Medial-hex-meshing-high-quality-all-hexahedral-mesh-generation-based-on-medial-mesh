#include "QuadDomain.h"

void read_quad_info(std::vector<Streamline>& streamlines, std::vector<SingularPoint>& points, std::vector<QuadDomain>& quad_domains) {

	streamlines.clear();
	points.clear();
	quad_domains.clear();
	int global_slid = 0, global_spid = 0;// 全局的流线id与奇异点id
	std::vector<std::pair<Streamline, int>> tmp_streamlines;
	std::vector<std::pair<SingularPoint, int>> tmp_singularpoints;
	bool input_ok = false;
	size_t dcnt = QInputDialog::getInt(nullptr, "", "Please input domain size",
		2, 1, 30, 1, &input_ok);
	if (!input_ok) return;

	for (int i = 0; i < dcnt; ++i) {

		std::vector<std::pair<Vex, Vex>> pos2_pos3;// 保存二维到三维点的映射关系

		// 1. 读入映射关系
		{
			std::string line;
			/*QString title = "Load mapping file";
			QString dir = "C:\\My Files\\Graphics\\model_data\\!test_data\\";
			QString filter = "(*.mapping.txt)";
			QString filename = QFileDialog::getOpenFileName(nullptr, title, dir, filter);
			
			if (filename.isEmpty()) break;
			std::ifstream mapping_inputfile(filename.toStdString(), std::ios::in);*/

			std::string filename = "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(i) + ".mapping.txt";
			std::ifstream mapping_inputfile(filename, std::ios::in);

			int pair_count = 0;
			{
				line.clear();
				getline(mapping_inputfile, line);
				std::stringstream linestream;
				linestream.str(line);
				linestream >> pair_count;
				std::cout << "Reading " << pair_count << " pairs position.\n";
			}
			for (int i = 0; i < pair_count; ++i) {
				line.clear();
				getline(mapping_inputfile, line);
				std::stringstream linestream;
				linestream.str(line);
				double x3, y3, z3, x2, y2, z2;
				linestream >> x3 >> y3 >> z3 >> x2 >> y2 >> z2;
				Vex pos3(x3, y3, z3);
				Vex pos2(x2, y2, z2);
				pos2_pos3.emplace_back(pos2, pos3);
			}
		}



		// 2. 读四边形区域信息
		{
			/*QString title_cdt = "Load cdt mesh file";
			QString dir_cdt = "C:\\My Files\\Graphics\\model_data\\!test_data\\";
			QString filter_cdt = "(*.2d.obj)";
			QString filename_cdt = QFileDialog::getOpenFileName(nullptr, title_cdt, dir_cdt, filter_cdt);
			if (filename_cdt.isEmpty()) return;*/
			std::string filename_2d = "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(i) + ".2d.obj";

			Project_KD_Tree_2D_To_3D kd_tree(filename_2d, pos2_pos3);


			std::string line;
			/*QString title = "Load quad info file";
			QString dir = "C:\\My Files\\Graphics\\model_data\\!test_data\\";
			QString filter = "(*.mlinequaInfo.txt)";
			QString filename = QFileDialog::getOpenFileName(nullptr, title, dir, filter);
			if (filename.isEmpty()) break;*/

			//std::string quad_file = "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(i) + ".mlinequaInfo.txt";
			std::string quad_file = "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(i) + ".mlinequaInfo.txt";
			std::ifstream quad_inputfile(quad_file, std::ios::in);

			std::unordered_map<int, int> streamline_loacl_to_global;
			std::unordered_map<int, int> singularpoint_local_to_global;
			int quad_cnt, sline_cnt, spoint_cnt;
			{
				line.clear();
				getline(quad_inputfile, line);
				std::stringstream linestream;
				linestream.str(line);
				linestream >> quad_cnt >> sline_cnt >> spoint_cnt;
				std::cout << "\nReading " << quad_file << ": quad_cnt = " << quad_cnt <<
					", sline_cnt = " << sline_cnt << ", spoint_cnt = " << spoint_cnt << std::endl;
			}

			for (int j = 0; j < quad_cnt; ++j) {// 读四边形区域划分信息
				line.clear();
				getline(quad_inputfile, line);
				std::stringstream linestream;
				linestream.str(line);
				int side, slid, dir;
				linestream >> side;
				if (side != 4) {
					std::cerr << "\nDoamine side != 4 !!!\n";
				}
				QuadDomain qd;
				for (int k = 0; k < side; ++k) {
					linestream >> slid >> dir;
					if (!streamline_loacl_to_global.count(slid)) {
						streamline_loacl_to_global[slid] = global_slid++;
					}
					qd.streamlines.push_back(streamline_loacl_to_global[slid]);
					qd.directions.push_back(dir);// 区域的奇异点未记录
				}
				quad_domains.push_back(qd);

			}

			for (int j = 0; j < sline_cnt; ++j) {// 读流线信息
				line.clear();
				getline(quad_inputfile, line);
				std::stringstream _linestream;
				_linestream.str(line);
				int slid, tmp, pcnt, beg_spid, end_spid;
				_linestream >> slid >> beg_spid >> end_spid >> tmp
					>> tmp >> tmp >> tmp >> pcnt;
				//std::cout << "Streamline#" << slid << ", points size = " << pcnt << ".\n";
				if (!singularpoint_local_to_global.count(beg_spid)) {
					singularpoint_local_to_global[beg_spid] = global_spid++;
				}
				if (!singularpoint_local_to_global.count(end_spid)) {
					singularpoint_local_to_global[end_spid] = global_spid++;
				}
				Streamline sl;
				sl.beg_sid = singularpoint_local_to_global[beg_spid];
				sl.end_sid = singularpoint_local_to_global[end_spid];
				double x, y, z;
				for (int k = 0; k < pcnt; ++k) {
					line.clear();
					getline(quad_inputfile, line);
					std::stringstream linestream;
					linestream.str(line);
					linestream >> x >> y >> z;
					Vex project_res = kd_tree.get_project_pos3(Vex(x, y, z));// 把二维点映射回三维

					sl.points.push_back(project_res);

				}
				if (!streamline_loacl_to_global.count(slid)) continue;
				tmp_streamlines.emplace_back(sl, streamline_loacl_to_global[slid]);
			}

			for (int j = 0; j < spoint_cnt; ++j) {// 读奇异点信息
				line.clear();
				getline(quad_inputfile, line);
				std::stringstream linestream;
				linestream.str(line);
				int spid, degree;
				double x, y, z;
				linestream >> spid >> degree >> x >> y >> z;
				Vex project_res = kd_tree.get_project_pos3(Vex(x, y, z));

				SingularPoint sp;
				sp.pos = project_res;
				//std::cout << "Singular point: (" << sp.pos.x() << ", " << sp.pos.y() << ", " << sp.pos.z() << ")" << std::endl;
				//points.push_back(sp);
				if (singularpoint_local_to_global.count(spid) == 0) std::cerr << "Error: exist local spid no global spid!\n";
				tmp_singularpoints.emplace_back(sp, singularpoint_local_to_global[spid]);
				for (int k = 0; k < degree; ++k) {// 读入该点流线id
					// 暂时不使用
					line.clear();
					getline(quad_inputfile, line);
					std::stringstream linestream;
					linestream.str(line);

				}
			}

		}

	}

	sort(tmp_streamlines.begin(), tmp_streamlines.end(), [&](std::pair<Streamline, int>& p0, std::pair<Streamline, int>& p1) {
		return p0.second < p1.second;
		});

	sort(tmp_singularpoints.begin(), tmp_singularpoints.end(), [&](std::pair<SingularPoint, int>& p0, std::pair<SingularPoint, int>& p1) {
		return p0.second < p1.second;
		});

	for (auto& sl : tmp_streamlines) {
		streamlines.push_back(sl.first);
	}

	for (auto& sl : streamlines) {
		sl.init_length_ratio();
	}

	for (auto& sp : tmp_singularpoints) {
		points.push_back(sp.first);
	}

	for (auto& qd : quad_domains) {

		// 设置流线对立关系
		std::unordered_map<int, bool> visited;
		for (int i = 0; i < qd.streamlines.size(); ++i) {
			if (visited.count(i)) continue;
			auto& sl0 = streamlines[qd.streamlines[i]];
			int sp00 = sl0.beg_sid, sp01 = sl0.end_sid;
			int j = 0;
			for (; j < 4; ++j) {
				if (i == j || visited.count(j)) continue;
				auto& sl1 = streamlines[qd.streamlines[j]];
				int sp10 = sl1.beg_sid, sp11 = sl1.end_sid;
				if (sp10 == sp00 || sp10 == sp01 || sp11 == sp00 || sp11 == sp01) continue;
				sl0.opposite_slid.push_back(qd.streamlines[j]);
				sl1.opposite_slid.push_back(qd.streamlines[i]);
				visited[i] = visited[j] = true;
				break;
			}
			if (j == 4) {
				std::cerr << "We cann't find the opposite streamline!\n";
			}
		}

		// 对流线重新排序
		for (int i = 1; i < 4; ++i) {
			std::vector<int> new_qdsl = qd.streamlines;
			std::vector<int> new_qddir = qd.directions;
			auto& op_slid = streamlines[new_qdsl[0]].opposite_slid;
			if (std::find(op_slid.begin(), op_slid.end(), qd.streamlines[i]) == op_slid.end()) {// 设置相邻的为第二条流线
				new_qdsl[1] = qd.streamlines[i];
				new_qddir[1] = qd.directions[i];
				for (int j = 1; j < 4; ++j) {
					if (i == j) continue;
					if (std::find(op_slid.begin(), op_slid.end(), qd.streamlines[j]) == op_slid.end()) {
						new_qdsl[3] = qd.streamlines[j];
						new_qddir[3] = qd.streamlines[j];
					} else {
						new_qdsl[2] = qd.streamlines[j];
						new_qdsl[2] = qd.streamlines[j];
					}
				}
				qd.streamlines = new_qdsl;
				qd.directions = new_qddir;
				break;
			}
		}

		for (int i = 0; i < 4; ++i) {
			std::cout << "[" << streamlines[qd.streamlines[i]].beg_sid << ", "
				<< streamlines[qd.streamlines[i]].end_sid << "]";
			if (i != 3) std::cout << ", ";
			else std::cout << "\n";
		}
		
		// 更新奇异点
		for (int i = 0; i < qd.streamlines.size(); ++i) {
			int global_slid_ = qd.streamlines[i];
			if (qd.directions[i]) {
				qd.singularpoints.push_back(streamlines[global_slid_].beg_sid);
			} else {
				qd.singularpoints.push_back(streamlines[global_slid_].end_sid);
			}
		}


	}

	std::cout << "Read quad info success.\n Streamlines size = " << streamlines.size() << ", singular points size = " << points.size()
		<< ", quad domains size = " << quad_domains.size() << std::endl;


}
