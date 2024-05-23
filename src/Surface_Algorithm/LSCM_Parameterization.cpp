#include "LSCM_Parameterization.h"
#include <queue>

bool LSCM_Parameterization::Execute(){
    std::function<MeshKernel::iGameVertexHandle(MeshKernel::iGameVertexHandle)> find_furthest_vertex
    = [&](MeshKernel::iGameVertexHandle vh){
        std::vector<int>dist(mesh.VertexSize());
        std::fill(dist.begin(),dist.end(),mesh.VertexSize() + 10);
        std::queue<MeshKernel::iGameVertexHandle>q;
        q.push(vh);
        dist[vh]=0;
        while(!q.empty()){
            MeshKernel::iGameVertexHandle now_v = q.front();
            q.pop();
            for(auto i : mesh.NeighborVh(now_v)){
                if(dist[i] > dist[now_v] + 1 ){
                    dist[i] = dist[now_v] + 1;
                    q.push(i);
                }
            }
        }
        MeshKernel::iGameVertexHandle maxx = vh;
        for(int i=0;i<mesh.VertexSize();i++){
            if(dist[maxx] <mesh.VertexSize() && dist[maxx] < dist[i]){
                maxx = MeshKernel::iGameVertexHandle(i);
            }
        }
        return maxx;
    };
    int fixedPoint1 = find_furthest_vertex(MeshKernel::iGameVertexHandle(0));
    int fixedPoint2 = find_furthest_vertex(MeshKernel::iGameVertexHandle(fixedPoint1));
    //std:: cout << "two fix vertex : " << fixedPoint1 <<" "<< fixedPoint2 << std :: endl;
    Execute(fixedPoint1,fixedPoint2);
    return true;
}

//固定2个顶点
bool LSCM_Parameterization::FixTwoBoundryPoint(int fixedPoint1, int fixedPoint2) {
	auto& vertices = mesh.allvertices();
	std::vector<MeshKernel::iGameVertexHandle> boundaryHandles;
	bool haveBoundary = false;

	for (auto& vp : vertices) {
		auto& vh = vp.first;
		auto& v = mesh.vertices(vh);
		if (mesh.isOnBoundary(vh)) {// 边界固定
			boundaryHandles.push_back(vh);
			haveBoundary = true;
		}
	}
	// 选取任意较远的两个点
	//F_P = boundaryHandles[0].idx();
	//S_P = boundaryHandles[boundaryHandles.size() / 2].idx();
	F_P = fixedPoint1;
	S_P = fixedPoint2;

	return haveBoundary;
}

//每个3D空间三角形对应的2D平面三角形的三条边，按面编号存储
void LSCM_Parameterization::CalculateEdgeVectors() {
	area.resize(f_count, 0.0);
	edgeVectors.resize(f_count);

	for (int k = 0; k < f_count; ++k) {
		edgeVectors[k].resize(3);//每行为3列
	}
	//遍历所有的面
	for (auto f : mesh.allfaces()) {
		std::vector<int> v(3, 0);
		int i = 0;
		auto fh = f.first;
		auto face = f.second;

		for (auto vh : face.getVertexHandle()) {
			v[i] = vh;   //保存该面顶点编号
			i++;
		}
		//计算3条边的边长
		double a, b, c;
		a = mesh.vertices(MeshKernel::iGameVertexHandle(v[0])).dist(mesh.vertices(MeshKernel::iGameVertexHandle(v[1])));
		b = mesh.vertices(MeshKernel::iGameVertexHandle(v[1])).dist(mesh.vertices(MeshKernel::iGameVertexHandle(v[2])));
		c = mesh.vertices(MeshKernel::iGameVertexHandle(v[2])).dist(mesh.vertices(MeshKernel::iGameVertexHandle(v[0])));
		double angle = acos((a * a + c * c - b * b) / (2 * a * c));
		double l = (a + b + c) / 2.0;
		area[fh.idx()] = sqrt(l * (l - a) * (l - b) * (l - c));   //三角形面积（海伦公式）

		//以U1为原点，构建在2D平面的三角形坐标
		Eigen::Vector2d U1, U2, U3;
		U1 << 0, 0;
		U2 << a, 0;
		U3 << c * cos(angle), c* sin(angle);

		//构建2D平面3条边
		std::vector< Eigen::Vector2d> triVectors(3);
		Eigen::Vector2d e1, e2, e3;
		e1 = U3 - U2;   //e1为U1的对边
		e2 = U1 - U3;   //e2为U2的对边
		e3 = U2 - U1;   //e3为U3的对边
		triVectors[0] = e1;
		triVectors[1] = e2;
		triVectors[2] = e3;

		//外层vector按面编号分组，内层vector存放该面的3条2D平面的边
		edgeVectors[fh.idx()] = triVectors;
	}
}

//构建A和BM矩阵
void LSCM_Parameterization::build_A2_BM() {
	A2.resize(2 * f_count, 2 * v_count - 4);
	int A2_face_offset = f_count;
	int A2_vertices_offset = v_count - 2;

	BM.resize(2 * f_count, 4);
	int BM_face_offset = f_count;
	int BM_vertices_offset = 2;
	//利用三元组进行稀疏矩阵插值
	std::vector<Eigen::Triplet<double>> tripletlist_A2;
	std::vector<Eigen::Triplet<double>> tripletlist_BM;
	//遍历每个三角形面
	for (auto f : mesh.allfaces()) {
		auto fh = f.first;
		auto face = f.second;

		//triEdgeVectors存放对应面序号的2D平面三角形的三条边
		std::vector<Eigen::Vector2d> triEdgeVectors = edgeVectors[fh.idx()];
		double triarea = area[fh.idx()];
		double denominator = sqrt(2 * triarea);
		int v_index[3];
		int i = 0;

		//遍历该面的所有点，逆时针0，1，2，获得点的编号
		for (auto vh : face.getVertexHandle()) {
			v_index[i] = vh.idx();
			i++;
		}

		//遍历2D平面上的该面的3条边
		for (i = 0; i < 3; i++) {
			Eigen::Vector2d edgeVector = triEdgeVectors[i];
			//x坐标（实部）
			double dx = edgeVector[0];
			int f_idx = fh.idx();
			int v_idx = VertexMapping[v_index[i]];
			//如果不是固定的点，就插入A矩阵的对应位置
			if (v_index[i] != F_P && v_index[i] != S_P) {

				tripletlist_A2.push_back(Eigen::Triplet<double>(f_idx, v_idx, dx / denominator));
				tripletlist_A2.push_back(Eigen::Triplet<double>(f_idx + A2_face_offset, v_idx + A2_vertices_offset, dx / denominator));
			}
			//如果是固定的点，就插入BM矩阵的对应位置
			else if (v_index[i] == F_P) {
				tripletlist_BM.push_back(Eigen::Triplet<double>(f_idx, 0, -dx / denominator));
				tripletlist_BM.push_back(Eigen::Triplet<double>(f_idx + BM_face_offset, BM_vertices_offset, dx / denominator));
			}
			else if (v_index[i] == S_P) {
				tripletlist_BM.push_back(Eigen::Triplet<double>(f_idx, 1, -dx / denominator));
				tripletlist_BM.push_back(Eigen::Triplet<double>(f_idx + BM_face_offset, 1 + BM_vertices_offset, dx / denominator));
			}
			//y坐标（虚部）
			double dy = edgeVector[1];
			if (v_index[i] != F_P && v_index[i] != S_P) {
				tripletlist_A2.push_back(Eigen::Triplet<double>(f_idx, v_idx + A2_vertices_offset, -dy / denominator));
				tripletlist_A2.push_back(Eigen::Triplet<double>(f_idx + A2_face_offset, v_idx, dy / denominator));
			}
			else if (v_index[i] == F_P) {
				tripletlist_BM.push_back(Eigen::Triplet<double>(f_idx, BM_vertices_offset, dy / denominator));
				tripletlist_BM.push_back(Eigen::Triplet<double>(f_idx + BM_face_offset, 0, -dy / denominator));
			}
			else if (v_index[i] == S_P) {
				tripletlist_BM.push_back(Eigen::Triplet<double>(f_idx, 1 + BM_vertices_offset, dy / denominator));
				tripletlist_BM.push_back(Eigen::Triplet<double>(f_idx + BM_face_offset, 1, -dy / denominator));
			}
		}
	}

	//三元组插入稀疏矩阵
	A2.setFromTriplets(tripletlist_A2.begin(), tripletlist_A2.end());
	BM.setFromTriplets(tripletlist_BM.begin(), tripletlist_BM.end());
}

// 矩阵转化为稀疏矩阵
void LSCM_Parameterization::buildSparseMatrix(Eigen::SparseMatrix<double>& A1_sparse, Eigen::MatrixXd A, int A_rows, int A_cols) {
	std::vector<Eigen::Triplet<double>> tripletlist;  //构建类型为三元组的vector
	for (int i = 0; i < A_rows; i++)
	{
		for (int j = 0; j < A_cols; j++)
		{
			if (std::fabs(A(i, j)) > 1e-6)
			{
				//按Triplet方式填充，速度快
				tripletlist.push_back(Eigen::Triplet<double>(i, j, A(i, j)));
				// 直接插入速度慢
				//A1_sparse.insert(i, j) = A(i, j);
			}
		}
	}
	A1_sparse.setFromTriplets(tripletlist.begin(), tripletlist.end());

}

//将顶点序号与坐标相对应
void LSCM_Parameterization::MapBackUV() {
	int VertexMapping_len = VertexMapping.size();
	int P_1 = antiVertexMapping[VertexMapping_len - 2];
	int P_2 = antiVertexMapping[VertexMapping_len - 1];
	New_U.resize(v_count);
	New_V.resize(v_count);
	New_U[P_1] = U[0];
	New_V[P_1] = U[2];
	New_U[P_2] = U[1];
	New_V[P_2] = U[3];
	int v_offset = v_count - 2;
	for (int i = 0; i < VertexMapping_len - 2; i++) {
		//例如Q在NewP中的序号为i，坐标为Q_Y，Q_Y;但它在实际空间的顶点的序号为antiVertexMapping[i],需要将时实际空间中的序号与新坐标相对应.
		int index = antiVertexMapping[i];
		double P_X = NewP(i);
		double P_Y = NewP(v_offset + i);
		New_U[index] = P_X;
		New_V[index] = P_Y;
	}
}


bool LSCM_Parameterization::Execute(int fixedPoint1, int fixedPoint2) {
	f_count = mesh.FaceSize();
	v_count = mesh.VertexSize();

	auto flag = FixTwoBoundryPoint(fixedPoint1, fixedPoint2);  //固定2个顶点
	if (flag == false) return false;
	CalculateEdgeVectors(); //每个3D空间三角形对应的2D平面三角形的三条边，按面编号存储

	VertexMapping.resize(v_count, 0.0);
	antiVertexMapping.resize(v_count, 0.0);
	//把两个固定的顶点顺序放到矩阵后面，用来做顶点顺序映射，分别定义两个保持正向和反向映射
	VertexMapping[F_P] = v_count - 2;
	VertexMapping[S_P] = v_count - 1;
	antiVertexMapping[v_count - 2] = F_P;
	antiVertexMapping[v_count - 1] = S_P;

	int itrCount = 0;
	//遍历所有的点,构建三角形顶点顺序正向映射与反向映射
	for (auto v : mesh.allvertices()) {
		auto vh = v.first;
		if (vh.idx() == F_P || vh.idx() == S_P) {
			continue;
		}
		VertexMapping[vh.idx()] = itrCount;
		antiVertexMapping[itrCount] = vh.idx();
		itrCount++;
	}

	build_A2_BM();  //构建A和BM矩阵

	//手动设置固定的2个点在平面的位置
	U[0] = 0.0;
	U[1] = 1.0;
	U[2] = 0.0;
	U[3] = 0.0;

	Eigen::MatrixXd bPart = BM * U;
	//std::cout << bPart << std::endl;
	B_Sparse.resize(2 * f_count, 1);
	buildSparseMatrix(B_Sparse, bPart, 2 * f_count, 1);   //将bPart转化成稀疏矩阵，方便进行运算
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > Solver_sparse;
	Solver_sparse.compute(A2);
	NewP.resize(2 * v_count - 4, 1);
	NewP = Solver_sparse.solve(B_Sparse);
	//std::cout << NewP << std::endl;

	MapBackUV();   //将顶点序号与坐标相对应

	for (auto& v : mesh.allvertices()) {
		auto& vh = v.first;
		auto& vex = mesh.vertices(vh);
		vex.setPosition(New_U[vh.idx()], New_V[vh.idx()], 0);
		//std::cout << New_U[vh.idx()] << " " << New_V[vh.idx()] << std::endl;
	}

	return true;
}