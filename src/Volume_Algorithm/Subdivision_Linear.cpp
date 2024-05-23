#include "Subdivision_Linear.h"

void Subdivision_Linear::subdivision(int iter_times, bool update_handle_flag) {

	std::time_t time_beg = clock();

	//if (iter_times > 1) update_handle_flag = true;

	for (uint16_t it = 0; it < iter_times; ++it) {

		auto old_vcnt = mesh.vsize(), old_ccnt = mesh.csize();

		std::time_t it_beg = clock();
		
		if (update_handle_flag) {
			add_new_cell();
			mesh.updateAllHandles();// �޸������˽ṹ��Ӧ�õ��ô˷���
		} else {
			add_new_cell_non_update_handle();
		}

		std::time_t it_end = clock();

		std::cout << "Subdivision: It " << it << " cost " << int(it_end - it_beg) << "ms. " <<
			"before vsize = " << old_vcnt << ", csize = " << old_ccnt << "; after vsize = " << mesh.vsize() << ", csize = " << mesh.csize() << std::endl;

	}
	std::time_t time_end = clock();
	std::cout << "Subdivision(Ordinary) success, cost " << int(time_end - time_beg) << "ms." << std::endl;

}

void Subdivision_Linear::add_new_cell() {

	int ccnt = mesh.csize(), fcnt = mesh.fsize(), ecnt = mesh.esize(), vcnt = mesh.vsize();

	std::vector<Vex> cell_v(ccnt);// ���
	std::vector<Vex> face_v(fcnt);// ���
	std::vector<Vex> edge_v(ecnt);// �ߵ�
	std::vector<Vex> face_c(fcnt);// ��¼ÿ���������

	// calculate cell points
#pragma omp parallel for
	for (int i = 0; i < ccnt; ++i) {
		cell_v[i] = mesh.getCellCenter((CH)i);
	}

	// calculate face points
#pragma omp parallel for
	for (int i = 0; i < fcnt; ++i) {

		FH fh(i);
		face_v[i] = face_c[i] = mesh.getFaceCenter(fh);
	}

	// calculate edge points
#pragma omp parallel for
	for (int i = 0; i < ecnt; ++i) {
		EH eh(i);
		edge_v[i] = mesh.getEdgeMidpoint(eh);// Ϊ�˷�����ϸ��
	}

	std::vector<VH> cell2vh, face2vh, edge2vh;
	cell2vh.reserve(ccnt);
	face2vh.reserve(fcnt);
	edge2vh.reserve(ecnt);
	// ���������������
	for (size_t i = 0; i < ccnt; ++i) {
		cell2vh.push_back(mesh.AddVertex(cell_v[i]));
	}
	for (size_t i = 0; i < fcnt; ++i) {
		face2vh.push_back(mesh.AddVertex(face_v[i]));
	}
	for (size_t i = 0; i < ecnt; ++i) {
		edge2vh.push_back(mesh.AddVertex(edge_v[i]));
	}

	// ����� Cell, ɾ���� Cell
	auto ori_cells = mesh.allcells();
	auto ori_edges = mesh.alledges();

	for (auto& cp : ori_cells) {
		auto cell_fhs = cp.second.getFaceHandle();
		auto cell_ehs = cp.second.getEdgeHandle();
		auto cell_vhs = cp.second.getVertexHandle();
		for (auto& cell_vh : cell_vhs) {// ���� Cell �е�ÿ���������� 1 ����

			std::vector<EH> adjehs;
			for (auto& cell_eh : cell_ehs) {
				auto& cell_e = mesh.edges(cell_eh);
				if (cell_e.vh1() == cell_vh || cell_e.vh2() == cell_vh) {
					adjehs.push_back(cell_eh);
					if (adjehs.size() == 3) break;
				}
			}
			if (adjehs.size() != 3) {
				std::cerr << "Subdivision Error: adjehs size is wrong." << std::endl;
				return;
			}
			std::vector<FH> adjfhs;
			for (int i = 0; i < 3; ++i) {
				auto tmp_eh1 = adjehs[i], tmp_eh2 = adjehs[(i + 1) % 3];
				for (auto& cell_fh : cell_fhs) {
					auto& face = mesh.faces(cell_fh);
					auto face_ehs = face.getEdgeHandle();
					int count = 0;
					for (auto& face_eh : face_ehs) {
						if (face_eh == tmp_eh1 || face_eh == tmp_eh2) count++;
					}
					if (count == 2) {
						adjfhs.push_back(cell_fh);
						break;
					}
				}
			}
			if (adjfhs.size() != 3) {
				std::cerr << "Subdivision Error: adjfhs size is wrong." << std::endl;
				return;
			}

			std::vector<VH> new_cell = {
				cell_vh, edge2vh[adjehs[1]], face2vh[adjfhs[1]], edge2vh[adjehs[2]],
				edge2vh[adjehs[0]], face2vh[adjfhs[0]], cell2vh[cp.first], face2vh[adjfhs[2]]
			};
			mesh.AddCell(new_cell);

		}
	}

	for (auto& ep : ori_edges) {
		mesh.DeleteEdge(ep.first);
	}

}

void Subdivision_Linear::add_new_cell_non_update_handle() {

	std::unordered_map<EH, VH> edge2vh;
	std::unordered_map<FH, VH> face2vh;
	std::unordered_map<CH, VH> cell2vh;
	
	for (auto& ep : mesh.alledges()) {
		edge2vh[ep.first] = mesh.AddVertex(mesh.getEdgeMidpoint(ep.first));
	}

	for (auto& fp : mesh.allfaces()) {
		face2vh[fp.first] = mesh.AddVertex(mesh.getFaceCenter(fp.first));
	}

	for (auto& cp : mesh.allcells()) {
		cell2vh[cp.first] = mesh.AddVertex(mesh.getCellCenter(cp.first));
	}

	// ����� Cell, ɾ���� Cell
	auto ori_cells = mesh.allcells();
	auto ori_edges = mesh.alledges();

	for (auto& cp : ori_cells) {
		auto cell_fhs = cp.second.getFaceHandle();
		auto cell_ehs = cp.second.getEdgeHandle();
		auto cell_vhs = cp.second.getVertexHandle();
		for (auto& cell_vh : cell_vhs) {// ���� Cell �е�ÿ���������� 1 ����

			std::vector<EH> adjehs;
			for (auto& cell_eh : cell_ehs) {
				auto& cell_e = mesh.edges(cell_eh);
				if (cell_e.vh1() == cell_vh || cell_e.vh2() == cell_vh) {
					adjehs.push_back(cell_eh);
					if (adjehs.size() == 3) break;
				}
			}
			if (adjehs.size() != 3) {
				std::cerr << "Subdivision Error: adjehs size is wrong." << std::endl;
				return;
			}
			std::vector<FH> adjfhs;
			for (int i = 0; i < 3; ++i) {
				auto tmp_eh1 = adjehs[i], tmp_eh2 = adjehs[(i + 1) % 3];
				for (auto& cell_fh : cell_fhs) {
					auto& face = mesh.faces(cell_fh);
					auto face_ehs = face.getEdgeHandle();
					int count = 0;
					for (auto& face_eh : face_ehs) {
						if (face_eh == tmp_eh1 || face_eh == tmp_eh2) count++;
					}
					if (count == 2) {
						adjfhs.push_back(cell_fh);
						break;
					}
				}
			}
			if (adjfhs.size() != 3) {
				std::cerr << "Subdivision Error: adjfhs size is wrong." << std::endl;
				return;
			}

			std::vector<VH> new_cell = {
				cell_vh, edge2vh[adjehs[1]], face2vh[adjfhs[1]], edge2vh[adjehs[2]],
				edge2vh[adjehs[0]], face2vh[adjfhs[0]], cell2vh[cp.first], face2vh[adjfhs[2]]
			};
			mesh.AddCell(new_cell);

		}
	}

	for (auto& ep : ori_edges) {
		mesh.DeleteEdge(ep.first);
	}


}