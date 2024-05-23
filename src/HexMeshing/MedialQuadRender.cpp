#include "MedialQuadRender.h"

void MedialQuadRender::get_streamlines(SkeletalMesh& sklmesh, vector<vector<Vex>>& res) {

	bool input_ok = false;// 获取面片域的个数
	int dcnt = QInputDialog::getInt(nullptr, "", "Please input domain size",
		1, 0, 30, 1, &input_ok);
	if (!input_ok) return;

	for (int di = 0; di < dcnt; ++di) {

		std::vector<std::pair<Vex, Vex>> pos2_pos3;// 保存二维到三维点的映射关系

		// 1. 读入映射关系
		{
			std::string line;
			std::string filename = "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(di) + ".mapping.txt";
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

		// 2. 读入二维流线信息
		std::vector<std::vector<Vex>> streamlines;

		{
			std::string line;
			std::string filename = "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(di) + ".quaInfo.txt";
			std::ifstream streamlines_inputfile(filename, std::ios::in);
			int quad_cnt, sline_cnt, nodes_cnt;
			std::unordered_set<int> fine_slid;
			{
				line.clear();
				getline(streamlines_inputfile, line);
				std::stringstream linestream;
				linestream.str(line);
				linestream >> quad_cnt >> sline_cnt >> nodes_cnt;
				std::cout << "Quad cnt = " << quad_cnt 
					<< ", sline cnt = " << sline_cnt 
					<< ", nodes cnt = " << nodes_cnt << std::endl;
				for (int qua_id = 0; qua_id < quad_cnt; ++qua_id) {// 四边区域与其流线关系
					line.clear();
					getline(streamlines_inputfile, line);
					std::stringstream linestream;
					linestream.str(line);
					int sl_cnt, sl_id;
					linestream >> sl_cnt;
					while (sl_cnt--) {
						linestream >> sl_id;
						fine_slid.insert(sl_id);
					}
				}

			}

			for (int sl_id = 0; sl_id < sline_cnt; ++sl_id) {
				int sline_idx, points_cnt;
				{
					line.clear();
					getline(streamlines_inputfile, line);
					std::stringstream linestream;
					linestream.str(line);
					linestream >> sline_idx >> points_cnt;
					std::cout << "sline_idx = " << sline_idx << ", points_cnt = " << points_cnt << std::endl;
				}
				std::vector<Vex> streamline;
				for (int pid = 0; pid < points_cnt; ++pid) {
					line.clear();
					getline(streamlines_inputfile, line);
					std::stringstream linestream;
					linestream.str(line);
					double x2, y2;
					linestream >> x2 >> y2;
					Vex point2d(x2, y2, 0);
					streamline.emplace_back(point2d);
				}
				if (fine_slid.count(sline_idx)) {
					streamlines.emplace_back(streamline);
				}
				
			}

		}

		// 3. 将二维流线转换至三维
		std::string filename_2d = "C:\\My Files\\Graphics\\model_data\\!test_data\\domain" + std::to_string(di) + ".2d.obj";
		Project_KD_Tree_2D_To_3D kd_tree(filename_2d, pos2_pos3);//.负责将二维顶点映射回三维

		for (auto& line : streamlines) {
			for (auto& vex : line) {
				//std::cout << "(" << vex.x() << "," << vex.y() << "," << vex.z() << ") --> ";
				vex = kd_tree.get_project_pos3(vex);// 将 2D 坐标转换到 3D
				//std::cout << "(" << vex.x() << "," << vex.y() << "," << vex.z() << ")\n";
			}
		}
		
		res.insert(res.end(), streamlines.begin(), streamlines.end());
		

	}

}