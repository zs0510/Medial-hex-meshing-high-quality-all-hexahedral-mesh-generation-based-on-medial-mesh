#include "Subdivision_CC.h"

void Subdivision_CC::subdivision_ordinary(uint16_t iter_times) {

	std::time_t time_beg = clock();
	

	for (uint16_t it = 0; it < iter_times; ++it) {

		auto old_vcnt = mesh.vsize(), old_ccnt = mesh.csize();

		std::time_t it_beg = clock();
		add_new_cell();
		mesh.updateAllHandles();// 修改了拓扑结构就应该调用此方法
		std::time_t it_end = clock();

		std::cout << "Subdivision: It " << it << " cost " << int(it_end - it_beg) << "ms. " << 
			"before vsize = " << old_vcnt << ", csize = " << old_ccnt << "; after vsize = " << mesh.vsize() << ", csize = " << mesh.csize() << std::endl;
		
	}
	std::time_t time_end = clock();
	std::cout << "Subdivision(Ordinary) success, cost " << int(time_end - time_beg) << "ms." << std::endl;
	
}

void Subdivision_CC::add_new_cell() {

	int ccnt = mesh.csize(), fcnt = mesh.fsize(), ecnt = mesh.esize(), vcnt = mesh.vsize();

	std::vector<Vex> cell_v(ccnt);// 体点
	std::vector<Vex> face_v(fcnt);// 面点
	std::vector<Vex> edge_v(ecnt);// 边点
	std::vector<Vex> face_c(fcnt);// 记录每个面的中心

	// calculate cell points
#pragma omp parallel for
	for (int i = 0; i < ccnt; ++i) {
		cell_v[i] = mesh.getCellCenter((CH)i);
	}

	// calculate face points
#pragma omp parallel for
	for (int i = 0; i < fcnt; ++i) {

		FH fh(i);
		if (mesh.isOnBoundary(fh)) {
			// CC Surface Subdivision
			face_v[i] = face_c[i] = mesh.getFaceCenter(fh);
		} else {
			// CC Volume Subdivision
			Vex vex(0, 0, 0);
			for (auto ch : mesh.NeighborCh(fh)) {
				vex += cell_v[ch];
			}
			face_c[i] = mesh.getFaceCenter(fh);
			vex += face_c[i] * 2;
			face_v[i] = vex / 4;
		}
		//face_v[i] = face_c[i] = mesh.getFaceCenter(fh);
	}

	// calculate edge points
#pragma omp parallel for
	for (int i = 0; i < ecnt; ++i) {

		EH eh(i);
		if (mesh.isOnBoundary(eh)) {
			Vex vex(0, 0, 0);
			auto fhs = mesh.NeighborFh(eh);
			int bdy_cnt = 0;
			for (auto& fh : fhs) {
				if (mesh.isOnBoundary(fh)) {
					vex += face_c[fh];
					bdy_cnt++;
				}
			}
			assert(bdy_cnt != 0);
			vex /= bdy_cnt;
			edge_v[eh] = (vex + mesh.getEdgeMidpoint(eh)) / 2;

		} else {

			auto chs = mesh.NeighborCh(eh);
			assert(!chs.empty());
			Vex c_avg(0, 0, 0);
			for (auto& ch : chs) {
				c_avg += cell_v[ch];
			}
			c_avg /= chs.size();

			auto fhs = mesh.NeighborFh(eh);
			assert(!fhs.empty());
			Vex f_avg(0, 0, 0);
			for (auto& fh : fhs) {
				// centroid of face
				f_avg += face_c[fh];
			}
			f_avg /= fhs.size();

			Vex e_mid = mesh.getEdgeMidpoint(eh);

			edge_v[i] = (c_avg + f_avg * 2 + e_mid) / 4;
		}
		//edge_v[i] = mesh.getEdgeMidpoint(eh);// 为了方正地细分

	}

	// update vertex positions
	 {
#pragma omp parallel for
		for (int i = 0; i < vcnt; ++i) {

			VH vh(i);
			auto& v = mesh.vertices(vh);

			if (mesh.isOnBoundary(vh)) {

				Vex f_avg(0, 0, 0);
				auto adjfhs = mesh.NeighborFh(vh);
				int cnt = 0;
				for (auto& fh : adjfhs) {
					if (mesh.isOnBoundary(fh)) {
						f_avg += face_c[fh];
						cnt++;
					}
				}
				f_avg /= cnt;

				Vex e_avg(0, 0, 0);
				auto adjehs = mesh.NeighborEh(vh);
				cnt = 0;
				for (auto& eh : adjehs) {
					if (mesh.isOnBoundary(eh)) {
						e_avg += mesh.getEdgeMidpoint(eh);
						cnt++;
					}
				}
				e_avg /= cnt;

				v = (f_avg + e_avg * 2 + v) / 4;

			} else {

				Vex c_avg(0, 0, 0);
				auto adjchs = mesh.NeighborCh(vh);
				for (auto& ch : adjchs) {
					c_avg += cell_v[ch];
				}
				c_avg /= adjchs.size();

				Vex f_avg(0, 0, 0);
				auto adjfhs = mesh.NeighborFh(vh);
				for (auto& fh : adjfhs) {
					f_avg += face_v[fh];
				}
				f_avg /= adjfhs.size();

				Vex e_avg(0, 0, 0);
				auto adjehs = mesh.NeighborEh(vh);
				for (auto& eh : adjehs) {
					e_avg += edge_v[eh];
				}
				e_avg /= adjehs.size();

				v = (c_avg + f_avg * 3 + e_avg * 3 + v) / 8;

			}

		}
	}

	std::vector<VH> cell2vh, face2vh, edge2vh;
	cell2vh.reserve(ccnt);
	face2vh.reserve(fcnt);
	edge2vh.reserve(ecnt);
	// 添加所有新增顶点
	for (size_t i = 0; i < ccnt; ++i) {
		cell2vh.push_back(mesh.AddVertex(cell_v[i]));
	}
	for (size_t i = 0; i < fcnt; ++i) {
		face2vh.push_back(mesh.AddVertex(face_v[i]));
	}
	for (size_t i = 0; i < ecnt; ++i) {
		edge2vh.push_back(mesh.AddVertex(edge_v[i]));
	}

	// 添加新 Cell, 删除旧 Cell
	auto ori_cells = mesh.allcells();
	auto ori_edges = mesh.alledges();

	for (auto& cp : ori_cells) {
		auto cell_fhs = cp.second.getFaceHandle();
		auto cell_ehs = cp.second.getEdgeHandle();
		auto cell_vhs = cp.second.getVertexHandle();
		for (auto& cell_vh : cell_vhs) {// 对于 Cell 中的每个顶点增加 1 个体

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