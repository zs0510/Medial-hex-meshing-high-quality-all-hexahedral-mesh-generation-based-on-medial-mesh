#include"Mesh.h"
#include"string"

// Mesh 定义
namespace MeshKernel { 
	
	/*=========================公共操作===============================*/
	void Mesh::initBBox() {
		double bbox_min_x = 99999999, bbox_min_y = 99999999, bbox_min_z = 99999999;
		double bbox_max_x = -99999999, bbox_max_y = -99999999, bbox_max_z = -99999999;
		for (auto& vp : vertices_) {
			bbox_min_x = std::min(bbox_min_x, vp.second.x());
			bbox_min_y = std::min(bbox_min_y, vp.second.y());
			bbox_min_z = std::min(bbox_min_z, vp.second.z());
			bbox_max_x = std::max(bbox_max_x, vp.second.x());
			bbox_max_y = std::max(bbox_max_y, vp.second.y());
			bbox_max_z = std::max(bbox_max_z, vp.second.z());
		}
		BBoxMin = iGameVertex(bbox_min_x, bbox_min_y, bbox_min_z);
		BBoxMax = iGameVertex(bbox_max_x, bbox_max_y, bbox_max_z);
		/*printf("iGame mesh: bounding-box: (%.3f, %.3f, %.3f) --> (%.3f, %.3f, %.3f)\n", 
			BBoxMin.x(), BBoxMin.y(), BBoxMin.z(), BBoxMax.x(), BBoxMax.y(), BBoxMax.z());*/
	}

	bool Mesh::isConnected(iGameVertexHandle vh1, iGameVertexHandle vh2) {
		for (auto vh : NeighborVh(vh1)) {
			if (vh == vh2) return true;
		}
		return false;
	}

	bool Mesh::isConnected(iGameEdgeHandle eh1, iGameEdgeHandle eh2) {
		if (!isValid(eh1) || !isValid(eh2)) return false;
		auto e1 = edges_[eh1];
		auto e2 = edges_[eh2];
		auto vh1 = e1.vh1(), vh2 = e1.vh2();
		auto vh3 = e2.vh1(), vh4 = e2.vh2();
		return (vh1 == vh3 || vh1 == vh4 || vh2 == vh3 || vh2 == vh4);
	}

	bool Mesh::isConnected(iGameFaceHandle fh1, iGameFaceHandle fh2) {
		if (!isValid(fh1) || !isValid(fh2)) return false;
		auto& face1 = faces_[fh1];
		auto& face2 = faces_[fh2];
		// 注意别用顶点，因为四边形是否邻接使用是否存在公共边来判定
		const auto& ehs1 = face1.getEdgeHandle();
		for (const auto& eh : face2.getEdgeHandle()) {
			if (std::find(ehs1.begin(), ehs1.end(), eh) != ehs1.end()) {
				return true;
			}
		}
		return false;
	}

	double Mesh::getLength(iGameEdgeHandle eh) {
		assert(edges_.count(eh));
		auto& edge = edges_[eh];
		auto& v1 = vertices_[edge.vh1()];
		auto& v2 = vertices_[edge.vh2()];
		return (v1 - v2).norm();
	}

	Vec Mesh::getFaceNormal(iGameFaceHandle fh) {
		assert(faces_.count(fh));
		auto& face = faces_[fh];
		int tricnt = 0;
		Vec res(0, 0 ,0);
		const auto& vhs = face.getVertexHandle();
		for (int i = 2; i < vhs.size(); ++i) {// 支持多边形的法向量均匀求解
			auto& v0 = vertices_[vhs[0]];
			auto& v1 = vertices_[vhs[i - 1]];
			auto& v2 = vertices_[vhs[i]];
			Vec vec01 = (v1 - v0).normalized();
			Vec vec02 = (v2 - v0).normalized();
			res += (vec01 % vec02).normalized();
			tricnt++;
		}
		res /= tricnt;
		return res.normalized();
	}

	void Mesh::genLength(iGameEdgeHandle eh) {
		auto& edge = edges_[eh];
		auto& v1 = vertices_[edge.vh1()];
		auto& v2 = vertices_[edge.vh2()];
		edge.setLength((v2 - v1).norm());
	}

	void Mesh::genAllEdgesLength() {
		for (auto& ep : this->edges_) {
			genLength(ep.first);
		}
	}

	iGameEdgeHandle Mesh::getEdgeHandle(iGameVertexHandle vh1, iGameVertexHandle vh2) {
		iGameEdgeHandle ret(-1);
		int vh_sum = vh1 + vh2;
		if (!isValid(vh1) || !isValid(vh2)) {
			std::cerr << "Error: There exist invalid vh.\n";
			return ret;
		}
		for (auto& eh : NeighborEh(vh1)) {
			auto& e = edges_[eh];
			if (vh_sum == e.vh1() + e.vh2()) {
				ret = eh;
				//std::cout << "Input: vh1 = " << vh1 << ", vh2 = " << vh2 << ", Output: vh1 = " << e.vh1() << ", vh2 = " << e.vh2() << std::endl;
				break;
			}
		}
		if (ret == -1) std::cerr << "Error: Return a invalid eh.\n";
		return ret;
	}

	iGameVertex Mesh::getEdgeMidpoint(iGameEdgeHandle eh) {
		if (!isValid(eh)) {
			std::cerr << "Edge handle is invalid!" << std::endl;
			return iGameVertex(0, 0, 0);
		}
		auto& edge = edges_[eh];
		auto& v1 = vertices_[edge.vh1()];
		auto& v2 = vertices_[edge.vh2()];
		return (v1 + v2) / 2;
	}

	iGameVertex Mesh::getFaceCenter(iGameFaceHandle fh) {
		if (!isValid(fh)) {
			std::cerr << "Face handle is invalid!" << std::endl;
			return iGameVertex(0, 0, 0);
		}
		auto& face = faces_[fh];
		auto vhs = face.getVertexHandle();
		iGameVertex center(0, 0, 0);
		for (auto& vh : vhs) {
			center += vertices_[vh];
		}
		return center / vhs.size();
	}


	/*=========================添加元素===============================*/
	iGameFaceHandle Mesh::AddFace(const std::vector<iGameVertexHandle>& _vhs) {
		//std::cout << " 现在加的面中的点为 : " << std::endl;
		//for (auto& vh : _vhs) {
		//	std::cout << vh << " : ";
		//}
		//std::cout << std::endl;
		std::vector<iGameEdgeHandle> ehs(_vhs.size());
		// 得到该面边的handle
		for (int i = 0; i < _vhs.size(); ++i) {
			if (i == 0) {
				ehs[i] = AddEdge(_vhs[_vhs.size() - 1], _vhs[i]);
			}
			else {
				ehs[i] = AddEdge(_vhs[i], _vhs[i - 1]);
			}
		}
		iGameFace f(_vhs, ehs);
		// 如果该面已经存在，则返回面的handle
		if (Face2Fh_.count(f)) {
			//std::cout << "该面已经存在 序号为 : " << Face2Fh_[f] << "以及面中的各个边的序号为 : " << std::endl;
			iGameFaceHandle fh = Face2Fh_[f];
			//for (int i = 0; i < ehs.size(); i++) {
			//	std::cout << ehs[i] << " ";
			//}
			//for (int i = 0; i < 4; i++) {
			//	std::cout << faces_[fh].eh(i) << " ";
			//}
			//std::cout << std::endl;
			// 仍然更新 handle 和 面直接的映射
			faces_[fh] = f;											// 建立handle与该面之间的映射
			return Face2Fh_[f];
		}
		// 否则
		else {
			iGameFaceHandle fh = GenFaceHandle();						// 生成一个新的handle
			//std::cout << " 新的面的序号为 : " << fh << " 和面中的各个边的序号为 : " << std::endl;
			//for (int i = 0; i < ehs.size(); i++) {
			//	std::cout << ehs[i] << " ";
			//}
			faces_[fh] = f;											// 建立handle与该面之间的映射
			Face2Fh_[f] = fh;										// 建立该面与handle之间的映射
			AddFace2Neighbor(fh);									// 将该面添加至面所包含的点和边的相邻面中
			//for (int i = 0; i < 4; i++) {
			//	std::cout << faces_[fh].eh(i) << " ";
			//}
			//std::cout << std::endl;
			return Face2Fh_[f];										// 返回该面的handle
		}
	}
	iGameEdgeHandle Mesh::AddEdge(const iGameVertexHandle& vh1, const iGameVertexHandle& vh2) {
		iGameEdge e(vh1, vh2);
		//std::cout << "边的两个端点的序号为 : " << vh1 << " : " << vh2 << std::endl;
		// 如果该边已经存在，则返回该边的handle
		if (Edge2Eh_.count(e)) {
			//std::cout << "此时该边已经存在 且边的序号为 : " << Edge2Eh_[e] << std::endl;
			return Edge2Eh_[e];
		}
		// 否则
		else {
			iGameEdgeHandle eh = GenEdgeHandle();						// 生成一个新的handle
			//std::cout << "生成边的序号为 : " << eh << std::endl;
			edges_[eh] = e;											// 建立handle与该边之间的映射
			Edge2Eh_[e] = eh;										// 建立该边与handle之间的映射
			AddEdge2Neighbor(eh);									// 将该边记录为其两个顶点的邻接边
			return Edge2Eh_[e];										// 返回该边的handle
		}
	}
	iGameVertexHandle Mesh::AddVertex(const iGameVertex& _v) {
		//if (Vertex2Vh_.count(_v)) return Vertex2Vh_[_v];
		iGameVertexHandle vh = GenVertexHandle();
		vertices_[vh] = _v;
		Vertex2Vh_[_v] = vh;
		return Vertex2Vh_[_v];
		
	}

	/*=========================删除元素===============================*/
	iGameVertexHandle Mesh::DeleteVertex(iGameVertexHandle _vh) {
		if (!vertices_.count(_vh)) return iGameVertexHandle(-1);        // 如果该点不存在，则返回-1
		else {
			// 删除相邻边元素（删除边时自然删除了相邻面）
			auto ve = NeighborEhOfVertex_[_vh];
			for (iGameEdgeHandle eh : ve) {
				DeleteEdge(eh);
			}
			// 删除顶点元素和以其为根据的所有邻接关系
			Vertex2Vh_.erase(vertices_[_vh]);
			vertices_.erase(_vh);
			NeighborEhOfVertex_.erase(_vh);
			NeighborFhOfVertex_.erase(_vh);
			return _vh;
		}

	}
	iGameEdgeHandle Mesh::DeleteEdge(iGameEdgeHandle _eh) {
		if (!edges_.count(_eh)) return iGameEdgeHandle(-1);             // 如果该边不存在，则返回-1
		else {
			// 删除邻接关系
			iGameEdge e(edges_[_eh]);
			for (int i = 0; i < 2; ++i) {
				iGameVertexHandle ev = e.vh(i);
				NeighborEhOfVertex_[ev].erase(_eh);                // 删除点邻接边
			}
			// 删除相邻面元素
			auto ef = NeighborFhOfEdge_[_eh];
			for (iGameFaceHandle fh : ef ) {
				DeleteFace(fh);
			}
			// 删除边元素和以其为根据的所有邻接关系
			Edge2Eh_.erase(edges_[_eh]);
			edges_.erase(_eh);
			NeighborFhOfEdge_.erase(_eh);
			return _eh;
		}
	}
	iGameFaceHandle Mesh::DeleteFace(iGameFaceHandle _fh) {
		if (!faces_.count(_fh)) return iGameFaceHandle(-1);             // 如果该面不存在，则返回-1
		else {                                                     // 如果该面存在，则返回删除的这个面的handle
			// 删除邻接关系
			iGameFace f(faces_[_fh]);
			for (int i = 0; i < f.size(); ++i) {
				iGameVertexHandle fv = f.vh(i);
				iGameEdgeHandle fe = f.eh(i);
				NeighborFhOfVertex_[fv].erase(_fh);               // 删除点邻接面
				NeighborFhOfEdge_[fe].erase(_fh);                 // 删除边邻接面
			}
			// 删除面元素
			Face2Fh_.erase(faces_[_fh]);
			faces_.erase(_fh);
			return _fh;
		}
	}

	/*=========================维护邻接关系===============================*/
	void Mesh::AddFace2Neighbor(const iGameFaceHandle& _fh)
	{
		iGameFace f = faces_[_fh];
		size_t n = f.size();
		for (int i = 0; i < n; ++i) {
			NeighborFhOfVertex_[f.vh(i)].insert(_fh);
		}
		for (int i = 0; i < n; ++i) {
			NeighborFhOfEdge_[f.eh(i)].insert(_fh);
		}
	}

	void Mesh::AddEdge2Neighbor(const iGameEdgeHandle& _eh)
	{
		iGameEdge e = edges_[_eh];
		NeighborEhOfVertex_[e.vh1()].insert(_eh);
		NeighborEhOfVertex_[e.vh2()].insert(_eh);

	}

	void Mesh::DeleteFace2Neighbor(const iGameFaceHandle& _fh){
		iGameFace f = faces_[_fh];
		size_t n = f.size();
		for (int i = 0; i < n; ++i) {
			NeighborFhOfVertex_[f.vh(i)].erase(_fh);
		}
		for (int i = 0; i < n; ++i) {
			NeighborFhOfEdge_[f.eh(i)].erase(_fh);
		}
	}
	void Mesh::DeleteEdge2Neighbor(const iGameEdgeHandle& _eh) {
		iGameEdge e = edges_[_eh];
		NeighborEhOfVertex_[e.vh1()].erase(_eh);
		NeighborEhOfVertex_[e.vh2()].erase(_eh);
	}

	Mesh & Mesh::operator=(const Mesh & _surfacemesh) {
		vertices_ = _surfacemesh.vertices_;
		edges_ = _surfacemesh.edges_;
		faces_ = _surfacemesh.faces_;
		Vertex2Vh_ = _surfacemesh.Vertex2Vh_;
		Edge2Eh_ = _surfacemesh.Edge2Eh_;
		Face2Fh_ = _surfacemesh.Face2Fh_;
		NeighborEhOfVertex_ = _surfacemesh.NeighborEhOfVertex_;
		NeighborFhOfVertex_ = _surfacemesh.NeighborFhOfVertex_;
		NeighborFhOfEdge_ = _surfacemesh.NeighborFhOfEdge_;
		VertexHandleID_ = _surfacemesh.VertexHandleID_;
		EdgeHandleID_ = _surfacemesh.EdgeHandleID_;
		FaceHandleID_ = _surfacemesh.FaceHandleID_;
		return *this;

	}

	/*=========================读写元素===============================*/
	// 读取ID为i的顶点
	iGameVertex& Mesh::vertices(iGameVertexHandle _vh) {
		assert(vertices_.count(_vh));
		return vertices_[_vh];
	}
	const iGameVertex Mesh::vertices(iGameVertexHandle _vh) const {
		assert(vertices_.count(_vh));
		return vertices_.find(_vh)->second;                // unordered_map 的 [] 操作符不是常量成员函数，无法对常量函数使用   
	}
	// 读取ID为i的边
	iGameEdge& Mesh::edges(iGameEdgeHandle _eh) {
		assert(edges_.count(_eh));
		return edges_[_eh];
	}
	const iGameEdge& Mesh::edges(iGameEdgeHandle _eh) const {
		assert(edges_.count(_eh));
		return edges_.find(_eh)->second;
	}
	// 读取ID为i的面
	iGameFace& Mesh::faces(iGameFaceHandle _fh) {
		assert(faces_.count(_fh));
		return faces_[_fh];
	}
	const iGameFace Mesh::faces(iGameFaceHandle _fh) const {
		assert(faces_.count(_fh));
		return faces_.find(_fh)->second;
	}

	/*====================根据元素得到对应ID=========================*/
	const iGameVertexHandle Mesh::vertexhanle(iGameVertex _vertex) const {
		if (Vertex2Vh_.find(_vertex) != Vertex2Vh_.end()) return Vertex2Vh_.find(_vertex)->second;
		else return iGameVertexHandle(-1);
	}
	const iGameEdgeHandle Mesh::edgehandle(iGameEdge& _edge) const {
		if (Edge2Eh_.find(_edge) != Edge2Eh_.end()) return Edge2Eh_.find(_edge)->second;
		else return iGameEdgeHandle(-1);
	}
	const iGameFaceHandle Mesh::facehandle(iGameFace& _face) const {
		if (Face2Fh_.find(_face) != Face2Fh_.end()) return Face2Fh_.find(_face)->second;
		else return iGameFaceHandle(-1);
	}

	/*======================得到邻接关系============================*/
	// 顶点的邻接点
	// 先找邻接边，再找邻接点
	std::unordered_set<iGameVertexHandle> Mesh::NeighborVh(iGameVertexHandle _vh) {
		std::unordered_set<iGameVertexHandle> neighborvh;
		auto neighboreh = NeighborEh(_vh);
		// 存在邻接边是前提
		if (neighboreh.size()) {
			for (iGameEdgeHandle eh : neighboreh) {
				if (edges_[eh].vh1() != _vh) neighborvh.insert(edges_[eh].vh1());
				if (edges_[eh].vh2() != _vh) neighborvh.insert(edges_[eh].vh2());
			}
		}
		return neighborvh;
	}
	// 顶点的邻接边
	std::unordered_set<iGameEdgeHandle>& Mesh::NeighborEh(iGameVertexHandle _vh) {
		if (NeighborEhOfVertex_.count(_vh)) return NeighborEhOfVertex_[_vh];
		else return empty_ehs;               // 返回一个空的集合
	}
	// 顶点的邻接面
	std::unordered_set<iGameFaceHandle>& Mesh::NeighborFh(iGameVertexHandle _vh) {
		if (NeighborFhOfVertex_.count(_vh)) return NeighborFhOfVertex_[_vh];
		else return empty_fhs;               // 返回一个空的集合
	}
	// 边的邻接边
	// 两个顶点的所有邻接边去除当前边
	std::unordered_set<iGameEdgeHandle> Mesh::NeighborEh(iGameEdgeHandle _eh) {
		assert(edges_.count(_eh));                     // 保证该边handle存在
		std::unordered_set<iGameEdgeHandle> neighboreh;     // 保存输出的结果
		int k = 0;                                     // 遍历两个顶点
		while (k < 2) {                                
			iGameVertexHandle vh = edges_[_eh].vh(k);
			auto vhneighboreh = NeighborEh(vh);        // 得到点的邻接边
			for (iGameEdgeHandle eh : vhneighboreh) {
				if (eh != _eh) neighboreh.insert(eh);
			}
			++k;                  
		}
		
		return neighboreh;
	}
	// 边的邻接面
	std::unordered_set<iGameFaceHandle>& Mesh::NeighborFh(iGameEdgeHandle _eh) {
		if (NeighborFhOfEdge_.count(_eh)) return NeighborFhOfEdge_[_eh];
		else return empty_fhs;               // 返回一个空的集合
	}
	// 面的邻接面
	// 邻接面：有一条相同边
	std::unordered_set<iGameFaceHandle> Mesh::NeighborFh(iGameFaceHandle _fh) {
		assert(faces_.count(_fh));                     // 保证该边handle存在
		std::unordered_set<iGameFaceHandle> neigborface;
		int k = 0;                                     // 遍历两个顶点
		size_t facesize = faces_[_fh].size();
		while (k < facesize) {
			iGameEdgeHandle eh = faces_[_fh].eh(k);
			auto ehneighborfh = NeighborFh(eh);        // 得到点的邻接边
			for (iGameFaceHandle fh : ehneighborfh) {
				if (fh != _fh) neigborface.insert(fh);
			}
			++k;
		}
		return neigborface;
	}
	// 邻接面2：有一个公共顶点
	std::unordered_set<iGameFaceHandle> Mesh::Neighbor2Fh(iGameFaceHandle _fh) {
		assert(faces_.count(_fh));                     // 保证该边handle存在
		std::unordered_set<iGameFaceHandle> neigborface;
		auto v_indices = faces_[_fh].getVertexHandle();
		for (auto& v_idx : v_indices) {
			auto adjF = NeighborFh(v_idx);
			for (iGameFaceHandle fh : adjF) {
				if (fh != _fh) neigborface.insert(fh);
			}
		}
		return neigborface;
	}
}



// SurfaceMesh 定义
namespace MeshKernel {
	void SurfaceMesh::InitMesh(const std::vector<iGameVertex>& _vertices,
		const std::vector<std::vector<iGameVertexHandle>>& _elements) {
		for (auto v : _vertices) {
			auto vh = AddVertex(iGameVertex(v.x(), v.y(), v.z()));

		}
		for (auto f : _elements) {
			AddFace(f);
		}
	}
	SurfaceMesh& SurfaceMesh::operator=(const SurfaceMesh& _surfacemesh) {
		if (this != &_surfacemesh) {
			Mesh::operator=(_surfacemesh);
		}
		return *this;
	}

	bool SurfaceMesh::isOnBoundary(iGameEdgeHandle eh) {
		auto fcnt = NeighborFh(eh).size();
		return fcnt == 1;
	}

	bool SurfaceMesh::isOnBoundary(iGameVertexHandle vh) {
		for (auto eh : NeighborEh(vh)) {
			if (isOnBoundary(eh)) {
				return true;
			}
		}
		return false;
	}

	bool SurfaceMesh::isOnBoundary(iGameFaceHandle fh) {
		auto face = faces(fh);
		for (auto eh : face.getEdgeHandle()) {
			if (isOnBoundary(eh)) {
				return true;
			}
		}
		return false;
	}

	bool SurfaceMesh::hasLoopBoundary() {

		std::unordered_map<EH, bool> visited;
		EH head_eh(-1);

		for (auto& ep : edges_) {
			if (isOnBoundary(ep.first)) {// 若存在唯一的一圈边界, 那么应当一次遍历完该网格所有边界边
				if (visited.count(ep.first)) continue;
				if (head_eh != -1) return false;// 存在多条边界
				head_eh = ep.first;
				EH pre_eh(-1), cur_eh(head_eh);
				do {
					visited[cur_eh] = true;
					int degree = 0;
					EH next_eh(-1);
					for (auto& adjeh : NeighborEh(cur_eh)) {
						if (isOnBoundary(adjeh)) {
							degree++;
							if (adjeh != pre_eh) {
								next_eh = adjeh;
							}
						}
					}
					if (degree != 2) return false;// 非流形
					pre_eh = cur_eh;
					cur_eh = next_eh;
				} while (cur_eh != head_eh);
				
			}
		}

		return head_eh != -1;

	}

	bool SurfaceMesh::isTriangleMesh() {
		for (auto& fp : faces_) {
			if (fp.second.getVertexHandle().size() != 3) return false;
		}
		return true;
	}

	std::vector<VH> SurfaceMesh::getOrderedBoundaryVH() {
		std::vector<VH> res;
		if (!this->hasLoopBoundary()) return res;
		for (auto& vp : vertices_) {
			if (isOnBoundary(vp.first)) {
				
				VH cur = vp.first, next(-1), pre(-1);
				do {
					res.push_back(cur);
					for (auto& adjvh : NeighborVh(cur)) {
						if (isOnBoundary(adjvh) && adjvh != pre) {
							next = adjvh;
							break;
						}
					}
					if (next == -1) {
						std::cerr << "Error in find loop boundary.\n";
						return {};
					}
					pre = cur;
					cur = next;
					
				} while (cur != res[0]);
				break;
			}
		}

		//std::cout << "We find " << res.size() << " vhs in loop boundary.\n";
		return res;
	}

	void SurfaceMesh::updateAllHandles() {
		int vcnt = VertexSize(), fcnt = FaceSize();
		std::vector<MeshKernel::iGameVertex> newVertices;
		std::vector<std::vector<MeshKernel::iGameVertexHandle>> newFaces;
		std::unordered_map<int, int> mp;// old id to new id
		int idx = 0;
		for (auto& fp : allfaces()) {
			auto vhs = fp.second.getVertexHandle();
			for (auto& vh : vhs) {
				if (!mp.count(vh)) {
					mp[vh] = idx++;
					newVertices.push_back(vertices_[vh]);
				}
				vh = iGameVertexHandle(mp[vh]);
			}
			newFaces.push_back(vhs);
		}

		*this = MeshKernel::SurfaceMesh(newVertices, newFaces);
	}

	void SurfaceMesh::eraseComplicatedEdges() {

		std::unordered_map<VH, std::unordered_set<VH>> vh_to_vhs;
		auto oldfps = this->faces_;

		std::unordered_map<FH, bool> visited;
		for (auto& fp : oldfps) {
			if (visited.count(fp.first)) continue;

			std::queue<FH> fhs;
			fhs.push(fp.first);
			visited[fp.first] = true;
			while (!fhs.empty()) {
				auto fh = fhs.front();
				fhs.pop();

				for (auto& adjfh : NeighborFh(fh)) {
					if (visited.count(adjfh)) continue;
					fhs.push(adjfh);
					visited[adjfh] = true;
				}

				auto& face = faces_[fh];
				auto vhs = face.getVertexHandle();
				bool need_readd = false;// 是否需要重新加面
				bool update_vh_to_vhs = true;
				for (int i = 0; i < vhs.size(); ++i) {
					int j = (i + 1) % vhs.size();
					if (vh_to_vhs[vhs[i]].count(vhs[j])) {// 发现重复边
						std::reverse(vhs.begin(), vhs.end());// 尝试使用翻转来修复
						need_readd = true;
						break;
					}
				}
				if (need_readd) {
					this->DeleteFace(fp.first);// 删去旧面
					bool cannot_add = false;
					for (int i = 0; i < vhs.size(); ++i) {
						int j = (i + 1) % vhs.size();
						if (vh_to_vhs[vhs[i]].count(vhs[j])) {// 无法修复，抛弃此面
							cannot_add = true;
							update_vh_to_vhs = false;
							std::cerr << "Find complicated edge cannot repair!\n";
							break;
						}
					}
					if (!cannot_add) this->AddFace(vhs);
				}
				if (update_vh_to_vhs) {// 将半边关系记录下来
					for (int i = 0; i < vhs.size(); ++i) {
						int j = (i + 1) % vhs.size();
						vh_to_vhs[vhs[i]].insert(vhs[j]);
					}
				}
			}

			

		}

		this->updateAllHandles();

	}

	void SurfaceMesh::flipAllFaces() {

		auto faces_origin = faces_;
		for (auto& fp : faces_origin) {
			FH fh(fp.first);
			auto vhs = fp.second.getVertexHandle();
			this->DeleteFace(fh);
			std::reverse(vhs.begin(), vhs.end());
			this->AddFace(vhs);
		}

		this->updateAllHandles();

	}

	void SurfaceMesh::destory() {// 清除所有数据并重置各种 handle

		this->vertices_.clear();
		this->VertexHandleID_ = 0;
		this->Vertex2Vh_.clear();
		this->NeighborEhOfVertex_.clear();
		this->NeighborFhOfVertex_.clear();

		this->edges_.clear();
		this->EdgeHandleID_ = 0;
		this->Edge2Eh_.clear();
		this->NeighborFhOfEdge_.clear();

		this->faces_.clear();
		this->FaceHandleID_ = 0;
		this->Face2Fh_.clear();

	}

	

	void SurfaceMesh::genNormal(iGameFaceHandle fh) {
		auto& face = faces(fh);
		auto vhs = face.getVertexHandle();
		assert(vhs.size() >= 3 && "should be a face not a line");

		auto& v0 = vertices_[vhs[0]];
		int face_normal_cnt = 0;
		Vec face_normal(0, 0, 0);
		for (int i = 2; i < vhs.size(); ++i) {
			auto& v1 = vertices_[vhs[i - 1]];
			auto& v2 = vertices_[vhs[i]];
			Vec vec01 = (v1 - v0).normalized();
			Vec vec02 = (v2 - v0).normalized();
			face_normal += ((vec01.cross(vec02)).normalized());
			face_normal_cnt++;
		}
		face_normal /= face_normal_cnt;
		face_normal.normalize();
		face.setNormal(face_normal.x(), face_normal.y(), face_normal.z());
	}

	void SurfaceMesh::genNormal(iGameVertexHandle vh) {
		auto& v = vertices(vh);
		Eigen::Vector3d N = Eigen::Vector3d::Zero();
		auto adjFH = NeighborFh(vh);
		for (auto& fh : adjFH) {
			auto& face = faces(fh);
			Eigen::Vector3d faceN = Eigen::Vector3d(face.getNormalX(), face.getNormalY(), face.getNormalZ());
			N += faceN;
		}
		N /= adjFH.size();
		N.normalize();
		v.setNormal(N[0], N[1], N[2]);
	}

	void SurfaceMesh::genAllFacesNormal() {
		for (auto& fp : this->faces_) {
			genNormal(fp.first);
		}
	}

	void SurfaceMesh::genAllVerticesNormal() {
		this->genAllFacesNormal();
		for (auto& vp : this->vertices_) {
			genNormal(vp.first);
		}
	}

	size_t SurfaceMesh::getBoundaryVerticesCount() {
		size_t cnt = 0;
		for (auto& vp : vertices_) {
			if (isOnBoundary(vp.first))
				cnt++;
		}
		return cnt;
	}

}


// VolumeMesh 成员函数定义
namespace MeshKernel {
	void VolumeMesh::InitMesh(const std::vector<iGameVertex>& _vertices,
		const std::vector<std::vector<iGameVertexHandle>>& _elements) {
		for (auto v : _vertices) {
			auto vh = AddVertex(iGameVertex(v.x(), v.y(), v.z()));
		}
		for (auto c : _elements) {
			AddCell(c);
		}
	}
	VolumeMesh& VolumeMesh::operator=(const VolumeMesh& _volumemesh) {
		vertices_ = _volumemesh.vertices_;
		edges_ = _volumemesh.edges_;
		faces_ = _volumemesh.faces_;
		cells_ = _volumemesh.cells_;
		Vertex2Vh_ = _volumemesh.Vertex2Vh_;
		Edge2Eh_ = _volumemesh.Edge2Eh_;
		Face2Fh_ = _volumemesh.Face2Fh_;
		Cell2Ch_ = _volumemesh.Cell2Ch_;
		NeighborEhOfVertex_ = _volumemesh.NeighborEhOfVertex_;
		NeighborFhOfVertex_ = _volumemesh.NeighborFhOfVertex_;
		NeighborChOfVertex_ = _volumemesh.NeighborChOfVertex_;
		NeighborFhOfEdge_ = _volumemesh.NeighborFhOfEdge_;
		NeighborChOfEdge_ = _volumemesh.NeighborChOfEdge_;
		NeighborChOfFace_ = _volumemesh.NeighborChOfFace_;
		VertexHandleID_ = _volumemesh.VertexHandleID_;
		EdgeHandleID_ = _volumemesh.EdgeHandleID_;
		FaceHandleID_ = _volumemesh.FaceHandleID_;
		CellHandleID_ = _volumemesh.CellHandleID_;
		surface_faces_ = _volumemesh.surface_faces_;
		return *this;
	}
	iGameCell& VolumeMesh::cells(iGameCellHandle _ch) {
		assert(cells_.count(_ch));
		return cells_[_ch];
	}
	const iGameCell VolumeMesh::cells(iGameCellHandle _ch) const {
		assert(cells_.count(_ch));
		return cells_.find(_ch)->second;
	}
	const iGameCellHandle VolumeMesh::cellhandle(iGameCell& _cell) const {
		if (Cell2Ch_.find(_cell) != Cell2Ch_.end()) return Cell2Ch_.find(_cell)->second;
		else return iGameCellHandle(-1);
	}
	std::unordered_set<iGameCellHandle> VolumeMesh::NeighborCh(iGameVertexHandle _vh) {
		if (NeighborChOfVertex_.count(_vh)) return NeighborChOfVertex_[_vh];
		else return std::unordered_set<iGameCellHandle>();
	}
	std::unordered_set<iGameCellHandle> VolumeMesh::NeighborCh(iGameEdgeHandle _eh) {
		if (NeighborChOfEdge_.count(_eh)) return NeighborChOfEdge_[_eh];
		else return std::unordered_set<iGameCellHandle>();
	}
	std::unordered_set<iGameCellHandle> VolumeMesh::NeighborCh(iGameFaceHandle _fh) {
		if (NeighborChOfFace_.count(_fh)) return NeighborChOfFace_[_fh];
		else return std::unordered_set<iGameCellHandle>();
	}
	std::unordered_set<iGameCellHandle> VolumeMesh::NeighborCh(iGameCellHandle _ch) {
		assert(cells_.count(_ch));                     // 保证该边handle存在
		std::unordered_set<iGameCellHandle> neigborcell;
		int k = 0;                                     // 遍历两个顶点
		size_t facesize = cells_[_ch].faces_size();
		while (k < facesize) {
			iGameFaceHandle fh = cells_[_ch].fh(k);
			auto fhneighborch = NeighborCh(fh);        // 得到点的邻接边
			for (iGameCellHandle ch : fhneighborch) {
				if (ch != _ch) neigborcell.insert(ch);
			}
			++k;
		}
		return neigborcell;
	}

	iGameCellHandle VolumeMesh::AddCell (const std::vector<iGameVertexHandle>& _vhs) {
		/*std::cout << "加入的体中的各个点为 : " << std::endl;
		for (int i = 0; i < _vhs.size(); i++) std::cout << _vhs[i] << " ";
		std::cout << std::endl;*/
		int facecnt = _vhs.size() == 8 ? 6 : 4;
		int edgecnt = _vhs.size() == 8 ? 12 : 6;
		auto vhs = _vhs;
		if (facecnt == 6) {
			// 保证顶点顺序 OK
			/*
			 * ATTENTION: All the metrics in this section are defined with the "irrational" order
			 *
			 *       P7------P6
			 *      / |     / |
			 *    P4------P5  |
			 *    |   |    |  |
			 *    |  P3----|--P2
			 *    | /      | /
			 *    P0------P1
			 *
			 */
			Vex face_center = (vertices_[_vhs[0]] + vertices_[_vhs[1]] + vertices_[_vhs[2]] + vertices_[_vhs[3]]) / 4.0;
			Vex cell_center = (vertices_[_vhs[4]] + vertices_[_vhs[5]] + vertices_[_vhs[6]] + vertices_[_vhs[7]] + face_center * 4.0) / 8.0;
			Vec vec_fc = (cell_center - face_center).normalized();
			Vec normal = ((vertices_[_vhs[1]] - vertices_[_vhs[0]]).cross(vertices_[_vhs[2]] - vertices_[_vhs[0]])).normalized();
			if (normal.dot(vec_fc) < 0) {
				std::reverse(vhs.begin(), vhs.begin() + 4);
				std::reverse(vhs.begin() + 4, vhs.end());
			}
		}
		std::vector<iGameFaceHandle> fhs(facecnt,(iGameFaceHandle)0);
		std::vector<std::vector<int>> faceform(facecnt);
		if (facecnt == 6) {
			faceform = { {0,3,2,1},{0,4,7,3},{0,1,5,4},{4,5,6,7},{1,2,6,5},{2,3,7,6}};
		}
		else {
			faceform = { {0,1,3},{1,2,3},{2,0,3},{0,2,1}};
		}
		for (int i = 0; i < facecnt; ++i) {
			std::vector<iGameVertexHandle> facevertices(faceform[i].size());
			for (int j = 0; j < faceform[i].size(); ++j) {
				facevertices[j] = vhs[faceform[i][j]];
			}
			fhs[i] = AddFace(facevertices);
		}
		std::vector<iGameEdgeHandle> ehs(edgecnt);
		////////////////////// Test Begin
		//std::cout << "该体的每一个面有的边的数量 以及 fhs[] 中的数  : " << std::endl;
		//for (int i = 0; i < facecnt; i++) {
		//	std::cout << std::to_string(i) << " : " << fhs[i] << " ";
		//	std::cout<< faces_[fhs[i]].getEN() << std::endl;
		//}
		//std::cout << std::endl;
		////////////////////// Test End
		if (edgecnt == 12) {
			ehs = { faces_[fhs[0]].eh(0),faces_[fhs[0]].eh(1),faces_[fhs[0]].eh(2),faces_[fhs[0]].eh(3),
			faces_[fhs[3]].eh(0),faces_[fhs[3]].eh(1),faces_[fhs[3]].eh(2),faces_[fhs[3]].eh(3),
			faces_[fhs[1]].eh(1), faces_[fhs[1]].eh(3), faces_[fhs[4]].eh(2), faces_[fhs[4]].eh(0) };
		} 
		else { 
			ehs = { faces_[fhs[3]].eh(0),faces_[fhs[3]].eh(1),faces_[fhs[3]].eh(2),
			faces_[fhs[0]].eh(0), faces_[fhs[1]].eh(0), faces_[fhs[2]].eh(0) };
		}
		iGameCell c(vhs, ehs, fhs);
		// 如果该体已经存在，则返回面的handle
		if (Cell2Ch_.count(c)) {
			//std::cout<<"该体已经存在 ." << std::endl;
			return Cell2Ch_[c];
		}
		// 否则
		else {
			iGameCellHandle ch = GenCellHandle();						// 生成一个新的handle
			cells_[ch] = c;											// 建立handle与该面之间的映射
			Cell2Ch_[c] = ch;										// 建立该面与handle之间的映射
			AddCell2Neighbor(ch);									// 将该面添加至面所包含的点和边的相邻面中
			return Cell2Ch_[c];										// 返回该面的handle
		}
	}

	iGameCellHandle VolumeMesh::AddCell(const std::vector< std::vector<iGameVertexHandle>>& _vhs) {
		if (_vhs.size() < 4) return iGameCellHandle(-1);
		std::vector<MeshKernel::iGameVertexHandle> vertices;
		std::vector<MeshKernel::iGameEdgeHandle> edges;
		std::vector<MeshKernel::iGameFaceHandle> faces;
		for (int i = 0; i < _vhs.size(); i++) {
			faces.push_back(AddFace(_vhs[i]));
			for (int j = 0; j < _vhs[i].size(); j++) {
				if (j == 0) {
					auto eh = AddEdge(_vhs[i][_vhs[i].size() - 1], _vhs[i][j]);
					auto it = std::find(edges.begin(), edges.end(), eh);
					// 防止handle重复
					if (it == edges.end())
						edges.push_back(eh);
				} else {
					auto eh = AddEdge(_vhs[i][j], _vhs[i][j - 1]);
					auto it = std::find(edges.begin(), edges.end(), eh);
					// 防止handle重复
					if (it == edges.end())
						edges.push_back(eh);
				}
				auto it = std::find(vertices.begin(), vertices.end(), _vhs[i][j]);
				// 防止handle重复
				if (it == vertices.end())
					vertices.push_back(_vhs[i][j]);
			}
		}
		iGameCell c(vertices, edges, faces);
		iGameCellHandle ch = GenCellHandle();
		cells_[ch] = c;
		return ch;
	}

	iGameVertexHandle VolumeMesh::DeleteVertex(const iGameVertexHandle& _vh) {
		if (!vertices_.count(_vh)) return iGameVertexHandle(-1);        // 如果该点不存在，则返回-1
		else {
			// 删除相邻边元素（删除边时自然删除了相邻面）
			auto ve = NeighborEhOfVertex_[_vh];
			for (iGameEdgeHandle eh : ve) {
				DeleteEdge(eh);
			}
			// 删除顶点元素和以其为根据的所有邻接关系
			Vertex2Vh_.erase(vertices_[_vh]);
			vertices_.erase(_vh);
			NeighborEhOfVertex_.erase(_vh);
			NeighborFhOfVertex_.erase(_vh);
			NeighborChOfVertex_.erase(_vh);
			return _vh;
		}
	}
	iGameEdgeHandle VolumeMesh::DeleteEdge(const iGameEdgeHandle& _eh) {
		if (!edges_.count(_eh)) return iGameEdgeHandle(-1);             // 如果该边不存在，则返回-1
		else {
			// 删除邻接关系
			iGameEdge e(edges_[_eh]);
			for (int i = 0; i < 2; ++i) {
				iGameVertexHandle ev = e.vh(i);
				NeighborEhOfVertex_[ev].erase(_eh);                // 删除点邻接边
			}
			// 删除相邻面元素
			auto ef = NeighborFhOfEdge_[_eh];
			for (iGameFaceHandle fh : ef) {
				DeleteFace(fh);
			}
			// 删除边元素和以其为根据的所有邻接关系
			Edge2Eh_.erase(edges_[_eh]);
			edges_.erase(_eh);
			NeighborFhOfEdge_.erase(_eh);
			NeighborChOfEdge_.erase(_eh);
			return _eh;
		}
	}
	iGameFaceHandle VolumeMesh::DeleteFace(const iGameFaceHandle& _fh) {
		if (!faces_.count(_fh)) return iGameFaceHandle(-1);             // 如果该面不存在，则返回-1
		else {                                                     // 如果该面存在，则返回删除的这个面的handle
			// 删除邻接关系
			iGameFace f(faces_[_fh]);
			for (int i = 0; i < f.size(); ++i) {
				iGameVertexHandle fv = f.vh(i);
				iGameEdgeHandle fe = f.eh(i);
				NeighborFhOfVertex_[fv].erase(_fh);               // 删除点邻接面
				NeighborFhOfEdge_[fe].erase(_fh);                 // 删除边邻接面
			}
			auto fc = NeighborChOfFace_[_fh];
			for (iGameCellHandle ch : fc) {
				DeleteCell(ch);
			}
			// 删除面元素
			Face2Fh_.erase(faces_[_fh]);
			faces_.erase(_fh);
			NeighborChOfFace_.erase(_fh);
			return _fh;
		}
	}
	iGameCellHandle VolumeMesh::DeleteCell(const iGameCellHandle& _ch) {
		if (!cells_.count(_ch)) return iGameCellHandle(-1);             
		else {                                                     
			// 删除邻接关系
			iGameCell c(cells_[_ch]);
			std::vector<int> vsize = { 4,8 };
			std::vector<int> esize = { 6,12 };
			std::vector<int> fsize = { 4,6 };
			int volumeType = c.vertices_size() == 4 ? 0 : 1;
			for (int i = 0; i < vsize[volumeType]; ++i) {
				iGameVertexHandle cv = c.vh(i);
				NeighborChOfVertex_[cv].erase(_ch);                               
			}
			for (int i = 0; i < esize[volumeType]; ++i) {
				iGameEdgeHandle ce = c.eh(i);
				NeighborChOfEdge_[ce].erase(_ch);
			}
			for (int i = 0; i < fsize[volumeType]; ++i) {
				iGameFaceHandle cf = c.fh(i);
				NeighborChOfFace_[cf].erase(_ch);
			}
			// 删除面元素
			Cell2Ch_.erase(cells_[_ch]);
			cells_.erase(_ch);
			return _ch;
		}
	}

	void VolumeMesh::AddCell2Neighbor(const iGameCellHandle& _ch)
	{
		iGameCell c(cells_[_ch]);
		std::vector<int> vsize = { 4,8 };
		std::vector<int> esize = { 6,12 };
		std::vector<int> fsize = { 4,6 };
		int volumeType = c.vertices_size() == 4 ? 0 : 1;
		for (int i = 0; i < vsize[volumeType]; ++i) {
			NeighborChOfVertex_[c.vh(i)].insert(_ch);
		}
		for (int i = 0; i < esize[volumeType]; ++i) {
			NeighborChOfEdge_[c.eh(i)].insert(_ch);
		}
		for (int i = 0; i < fsize[volumeType]; ++i) {
			NeighborChOfFace_[c.fh(i)].insert(_ch);
		}
	}

	void VolumeMesh::DeleteCell2Neighbor(const iGameCellHandle& _ch) {
		if (!cells_.count(_ch)) return;
		iGameCell c(cells_[_ch]);
		for (auto& fh : c.getFaceHandle()) {
			NeighborChOfFace_[fh].erase(_ch);
		}
		for (auto& eh : c.getEdgeHandle()) {
			NeighborChOfEdge_[eh].erase(_ch);
		}
		for (auto& vh : c.getVertexHandle()) {
			NeighborChOfVertex_[vh].erase(_ch);
		}
	}

	bool VolumeMesh::isOnBoundary(iGameCellHandle ch) {
		assert(cells_.count(ch));
		auto fhs = cells_[ch].getFaceHandle();
		for (auto& fh : fhs) {
			if (NeighborCh(fh).size() == 1) {
				return true;
			}
		}
		return false;
	}

	bool VolumeMesh::isOnBoundary(iGameFaceHandle fh) {
		assert(faces_.count(fh));
		return NeighborCh(fh).size() == 1;// 0 isnot on the boundary
	}

	bool VolumeMesh::isOnBoundary(iGameEdgeHandle eh) {
		assert(edges_.count(eh));
		for (auto& fh : NeighborFhOfEdge_[eh]) {
			if (NeighborCh(fh).size() == 1) {
				return true;
			}
		}
		return false;
	}

	bool VolumeMesh::isOnBoundary(iGameVertexHandle vh) {
		assert(vertices_.count(vh));
		for (auto& fh : NeighborFhOfVertex_[vh]) {
			if (NeighborCh(fh).size() == 1) {
				return true;
			}
		}
		return false;
	}

	iGameVertex VolumeMesh::getCellCenter(iGameCellHandle ch) {
		iGameVertex center(0, 0, 0);
		const auto& cell = cells_[ch];
		auto vhs = cell.getVertexHandle();
		for (auto& vh : vhs) {
			center = center + vertices_[vh];
		}
		return center / vhs.size();
	}

	double VolumeMesh::getQuadArea(iGameFaceHandle fh) {
		const auto& face = faces_[fh];
		auto vhs = face.getVertexHandle();
		assert(vhs.size() == 4);
		double quad_area = 0;
		const auto& v1 = vertices_[vhs[0]];
		const auto& v2 = vertices_[vhs[1]];
		const auto& v3 = vertices_[vhs[2]];
		const auto& v4 = vertices_[vhs[3]];
		auto vec13 = v3 - v1;
		auto vec12 = v2 - v1;
		auto vec14 = v4 - v1;
		quad_area += (vec12 % vec13).norm() * 0.5;
		quad_area += (vec14 % vec13).norm() * 0.5;
		quad_area += (vec12 % vec14).norm() * 0.5;
		auto vec23 = v3 - v2;
		auto vec24 = v4 - v2;
		quad_area += (vec23 % vec24).norm() * 0.5;
		return quad_area / 2;
	}

	iGameVertex VolumeMesh::getQuadNormal(iGameFaceHandle fh) {
		iGameVertex N(0, 0, 0);
		const auto& face = faces_[fh];
		auto vhs = face.getVertexHandle();
		auto chs = NeighborCh(fh);
		if (chs.size() != 1) return N;
		assert(vhs.size() == 4);
		auto cell_cenetr = getCellCenter(*chs.begin());
		//double quad_area = 0;
		const auto& v1 = vertices_[vhs[0]];
		const auto& v2 = vertices_[vhs[1]];
		const auto& v3 = vertices_[vhs[2]];
		const auto& v4 = vertices_[vhs[3]];
		auto face_center = (v1 + v2 + v3 + v4) / 4;
		auto out_dir = face_center - cell_cenetr;
		out_dir = out_dir.normalized();
		auto vec13 = v3 - v1;
		auto vec12 = v2 - v1;
		auto vec14 = v4 - v1;
		auto vec23 = v3 - v2;
		auto vec24 = v4 - v2;

		auto N123 = vec12 % vec13;
		double area1 = N123.norm() * 0.5;
		N123 = N123.normalized();
		if (N123 * out_dir < 0) N123 = N123 * -1;

		auto N134 = vec14 % vec13;
		double area2 = N134.norm() * 0.5;
		N134 = N134.normalized();
		if (N134 * out_dir < 0) N134 = N134 * -1;

		auto N124 = vec12 % vec14;
		double area3 = N124.norm() * 0.5;
		N124 = N124.normalized();
		if (N124 * out_dir < 0) N124 = N124 * -1;

		auto N234 = vec23 % vec24;
		double area4 = N234.norm() * 0.5;
		N234 = N234.normalized();
		if (N234 * out_dir < 0) N234 = N234 * -1;
		N = (N123 * area1 + N134 * area2 + N124 * area3 + N234 * area4);
		N = N.normalized();
		return N;
	}

	void VolumeMesh::updateAllHandles() {

		std::vector<MeshKernel::iGameVertex> newVertices;
		std::vector<std::vector<MeshKernel::iGameVertexHandle>> newCells;
		std::unordered_map<int, int> mp;// old id to new id
		int idx = 0;
		for (auto& cp : cells_) {
			auto vhs = cp.second.getVertexHandle();
			for (auto& vh : vhs) {
				if (!mp.count(vh)) {
					mp[vh] = idx++;
					newVertices.push_back(vertices_[vh]);
				}
				vh = iGameVertexHandle(mp[vh]);
			}
			newCells.push_back(vhs);
		}

		*this = MeshKernel::VolumeMesh(newVertices, newCells);

	}

	void VolumeMesh::generateAllBoundaryVertexNormals() {

		generateAllBoundaryFaceNormals();
		for (auto& vp : vertices_) {
			if (!isOnBoundary(vp.first)) continue;
			const auto& fhs = NeighborFhOfVertex_[vp.first];
			int bd_fcnt = 0;
			Vec N(0, 0, 0);
			for (auto& fh : fhs) {
				if (!isOnBoundary(fh)) continue;
				auto& f = faces_[fh];
				N += f.getNormal();
				bd_fcnt++;
			}
			if (bd_fcnt == 0) continue;
			N.normalize();
			vp.second.setNormal(N.x(), N.y(), N.z());
		}

	}

	void VolumeMesh::generateAllBoundaryFaceNormals() {

		for (auto& fp : faces_) {
			if (!isOnBoundary(fp.first)) continue;
			const auto& vhs = fp.second.getVertexHandle();
			CH ch = *(NeighborChOfFace_[fp.first].begin());
			Vex center_c = getCellCenter(ch);
			Vex center_f = getFaceCenter(fp.first);
			Vec dir_out = (center_f - center_c).normalized();
			Vec N(0, 0, 0);
			for (int i = 2; i < vhs.size(); ++i) {
				auto& v0 = vertices_[vhs[0]];
				auto& v1 = vertices_[vhs[i - 1]];
				auto& v2 = vertices_[vhs[i - 2]];
				Vec vec01 = (v1 - v0).normalized();
				Vec vec02 = (v2 - v0).normalized();
				Vec nor = (vec01.cross(vec02)).normalized();
				if (nor.dot(dir_out) < 0) {// 保证法向朝外
					nor *= -1;
				}
				N += nor;
			}
			N.normalized();
			fp.second.setNormal(N[0], N[1], N[2]);
		}

	}

	void VolumeMesh::destory() {

		this->vertices_.clear();
		this->VertexHandleID_ = 0;
		this->Vertex2Vh_.clear();
		this->NeighborEhOfVertex_.clear();
		this->NeighborFhOfVertex_.clear();
		this->NeighborChOfVertex_.clear();

		this->edges_.clear();
		this->EdgeHandleID_ = 0;
		this->Edge2Eh_.clear();
		this->NeighborFhOfEdge_.clear();
		this->NeighborChOfEdge_.clear();

		this->faces_.clear();
		this->FaceHandleID_ = 0;
		this->Face2Fh_.clear();
		this->NeighborChOfFace_.clear();

		this->cells_.clear();
		this->CellHandleID_ = 0;
		this->Cell2Ch_.clear();
		

	}

}

// Tetrahedron 成员函数定义
namespace MeshKernel {

	void TetMesh::InitHedra(const std::vector<iGameVertex>& _vertices, std::vector<std::vector<std::vector<iGameVertexHandle>>>& _cells,
		const std::vector<std::vector<iGameVertexHandle>>& _surface_faces) {
		for (auto v : _vertices) {
			auto vh = AddVertex(iGameVertex(v.x(), v.y(), v.z()));
		}
		for (auto ce : _cells) {
			AddCell(ce);
		}
		surface_faces_ = _surface_faces;
	}

	void TetMesh::genNormal(iGameFaceHandle fh) {
		auto& face = faces(fh);
		auto vex = face.getVertexHandle();
		std::vector<Eigen::Vector3d> pos(3, Eigen::Vector3d::Zero());
		for (int i = 0; i < 3; ++i) {
			auto& v = vertices(vex[i]);
			pos[i][0] = v.x();
			pos[i][1] = v.y();
			pos[i][2] = v.z();
		}
		Eigen::Vector3d N = (pos[1] - pos[0]).cross(pos[2] - pos[0]);
		N.normalize();
		face.setNormal(N[0], N[1], N[2]);
	}

	void TetMesh::genNormal(iGameVertexHandle vh) {
		auto& v = vertices(vh);
		Eigen::Vector3d N = Eigen::Vector3d::Zero();
		auto adjFH = NeighborFh(vh);
		for (auto& fh : adjFH) {
			auto& face = faces(fh);
			Eigen::Vector3d faceN = Eigen::Vector3d(face.getNormalX(), face.getNormalY(), face.getNormalZ());
			N += faceN;
		}
		N /= adjFH.size();
		N.normalize();
		v.setNormal(N[0], N[1], N[2]);
	}

	void TetMesh::genAllFacesNormal() {
		for (auto& fp : this->faces_) {
			if (isOnBoundary(fp.first)) {
				genNormal(fp.first);
			}
		}
	}

	void TetMesh::genAllVerticesNormal() {
		this->genAllFacesNormal();
		for (auto& vp : this->vertices_) {
			if (isOnBoundary(vp.first)) {
				genNormal(vp.first);
			}
		}
	}

}

namespace MeshKernel {
	
}


