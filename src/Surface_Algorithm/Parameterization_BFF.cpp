#include "Parameterization_BFF.h"

#define PI acos(-1)

void bff::Execute() {
	printf("boundary first flattening is starting\n");
	Init();
	Laplace();
	InitCurvature();
	Parameterization();
	printf("boundary first flattening success\n");
}

MeshKernel::EdgeHandle bff::NextEdge(MeshKernel::VertexHandle vh) {
	for (auto& eh : mesh.NeighborEh(vh)) {
		if (mesh.isOnBoundary(eh)) {
			if (!vis[eh]) {
				return eh;
			}
		}
	}
	return MeshKernel::EdgeHandle(mesh.EdgeSize() + 10);
	assert("something wrong" && 1 == 0);
}

std::vector<MeshKernel::EdgeHandle> bff::AroundEdge(MeshKernel::VertexHandle vh) {
	std::set<int> s;
	std::map<int, MeshKernel::EdgeHandle> vtoe;
	std::vector<MeshKernel::EdgeHandle> res;
	int start = -1;
	MeshKernel::VertexHandle cur, nxt;
	for (auto& eh : mesh.NeighborEh(vh)) {
		auto vh2 = mesh.NeighborVhFromEdge(vh, eh);
		if (start == -1 || !mesh.isOnBoundary(MeshKernel::VertexHandle(start))) start = vh2;
		s.insert(vh2);
		vtoe[vh2] = eh;
	}
	cur = MeshKernel::VertexHandle(start);
	while (1) {
		res.push_back(vtoe[cur]);
		s.erase(cur);
		bool ok = 0;
		for (auto& eh : mesh.NeighborEh(cur)) {
			if (s.count(mesh.NeighborVhFromEdge(cur, eh))) {
				nxt = mesh.NeighborVhFromEdge(cur, eh);
				ok = 1;
				break;
			}
		}
		if (!ok) break;
		cur = nxt;
		if (cur == MeshKernel::VertexHandle(start)) break;
	}
	return res;
}

void bff::Parameterization() {
	// 第一步 计算共形因子u
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
	Eigen::SparseMatrix<double> Aii_sparse = Aii.sparseView();
	solver.compute(Aii_sparse);
	//Eigen::MatrixXd ui = Aii.inverse() * -K;
	Eigen::MatrixXd ui = solver.solve(-K);
	Eigen::MatrixXd ub = Eigen::MatrixXd::Zero(bn, 1);
	Eigen::MatrixXd u(n, 1);
	u << ui,
		ub;
	// for (int i = 0; i < in; i++) std::cout << ui(i, 0) << std::endl;
	// 第二步 计算新边界点的角度k_tilde
	Eigen::MatrixXd h = -Abi * ui;
	Eigen::MatrixXd k_tilde = k - h;
	// 第三步：计算边界gamma
	MeshKernel::VertexHandle s = MeshKernel::VertexHandle(0), t = MeshKernel::VertexHandle(0);
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			s = vp.first;
			break;
		}
	}
	MeshKernel::VertexHandle s2 = s; // record
	

	// 微调长度使得边界闭合
	Eigen::MatrixXd N = Eigen::MatrixXd::Zero(bn, bn);
	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(2, bn);
	Eigen::MatrixXd len = Eigen::MatrixXd::Zero(bn, 1), len_tilde = Eigen::MatrixXd::Zero(bn, 1);



	double angle = 0;
	while (1) {
		auto eh = NextEdge(s);
		if ((int)eh > mesh.EdgeSize()) break;
		t = mesh.NeighborVhFromEdge(s, eh);
		vis[eh] = 1;
		angle += k_tilde(to_new[s] - in, 0);
		T(0, to_new[s] - in) = cos(angle), T(1, to_new[s] - in) = sin(angle);
		len(to_new[s] - in, 0) = exp(0.5 * (u(to_new[s], 0) + u(to_new[t], 0))) *
			(mesh.vertices(mesh.edges(eh).vh1()) - mesh.vertices(mesh.edges(eh).vh2())).norm();
		N(to_new[s] - in, to_new[s] - in) = 1 / len(to_new[s] - in, 0);
		s = t;
	}
	len_tilde = len - (N.inverse() * T.transpose() * (T * N.inverse() * T.transpose()).inverse() * T) * len;
	vis.clear();
	angle = 0;
	s = s2;
	std::map<int, MeshKernel::Vertex> gamma;
	gamma[s] = MeshKernel::Vertex(0, 0, 0);
	while (1) {
		auto eh = NextEdge(s);
		if ((int)eh > mesh.EdgeSize()) break;
		t = mesh.NeighborVhFromEdge(s, eh);
		vis[eh] = 1;
		angle += k_tilde(to_new[s] - in, 0);
		gamma[t] = gamma[s] + MeshKernel::Vertex(T(0, to_new[s] - in), T(1, to_new[s] - in), 0) * len_tilde(to_new[s] - in, 0);
		s = t;
	}

	//std::cout << "bn:" << bn << " " << gamma.size() << std::endl;
	//std::cout << " total angle:" << angle*180/PI << std::endl;
	Eigen::MatrixXd gamma_re = Eigen::MatrixXd::Zero(bn, 1);
	Eigen::MatrixXd gamma_im = Eigen::MatrixXd::Zero(bn, 1);
	for (auto x : gamma) {
		gamma_re(to_new[x.first] - in, 0) = x.second.x();
		gamma_im(to_new[x.first] - in, 0) = x.second.y();
	}
	/*
	for (auto x : gamma) {
		std::cout << x.second.x() << " " << x.second.y() << " " << x.second.z() << std::endl;
	}*/

	// 第四步：延拓gamma为f
	//Eigen::MatrixXd ai = -Aii.inverse() * Aib * gamma_re;
	//Eigen::MatrixXd bi = -Aii.inverse() * Aib * gamma_im;
	Eigen::MatrixXd ai = solver.solve(-Aib * gamma_re);
	Eigen::MatrixXd bi = solver.solve(-Aib * gamma_im);
	Eigen::MatrixXd a(n, 1);
	a << ai,
		gamma_re;
	Eigen::MatrixXd b(n, 1);
	b << bi,
		gamma_im;
	for (int i = 0; i < n; i++) {
		auto vh = MeshKernel::VertexHandle(to_old[i]);
		mesh.vertices(vh).setPosition(a(i, 0), b(i, 0), 0.0);
	}
}

void bff::InitCurvature() {
	K = Eigen::MatrixXd::Zero(in, 1);
	k = Eigen::MatrixXd::Zero(bn, 1);
	for (auto& vp : mesh.allvertices()) {
		double defect = 0;
		bool first = 1;
		MeshKernel::Vertex last;
		if (!mesh.isOnBoundary(vp.first)) {
			for (auto& eh : AroundEdge(vp.first)) {
				MeshKernel::Vertex cur;
				auto vh2 = mesh.NeighborVhFromEdge(vp.first, eh);
				cur = mesh.vertices(vh2) - vp.second;
				if (!first) {
					defect += acos(cur * last / (cur.norm() * last.norm()));
				}
				last = cur;
				first = 0;
			}
			auto vh2 = mesh.NeighborVhFromEdge(vp.first, *AroundEdge(vp.first).begin());
			auto cur = mesh.vertices(vh2) - vp.second;
			defect += acos(cur * last / (cur.norm() * last.norm()));
			defect = 2 * PI - defect;
			K(to_new[vp.first], 0) = defect;

		}
		else {
			for (auto& eh : AroundEdge(vp.first)) {
				MeshKernel::Vertex cur;
				auto vh2 = mesh.NeighborVhFromEdge(vp.first, eh);
				cur = mesh.vertices(vh2) - vp.second;
				if (!first) {
					defect += acos(cur * last / (cur.norm() * last.norm()));
				}
				last = cur;
				first = 0;
			}
			defect = PI - defect;
			k(to_new[vp.first] - in, 0) = defect;
		}
	}
}

void bff::Init() {
	n = mesh.VertexSize();
	for (auto& vp : mesh.allvertices()) {
		if (!mesh.isOnBoundary(vp.first)) {
			points.push_back(vp.first);
			to_new[vp.first] = points.size() - 1;
			to_old[points.size() - 1] = vp.first;
		}
	}
	in = points.size();
	bn = n - in;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			points.push_back(vp.first);
			to_new[vp.first] = points.size() - 1;
			to_old[points.size() - 1] = vp.first;
		}
	}
}

inline Eigen::MatrixXd submatrix(Eigen::MatrixXd const& A, int x, int y, int h, int w) {
	Eigen::MatrixXd B(h, w);
	for (int i = 0; i < h; i++) for (int j = 0; j < w; j++) {
		B(i, j) = A(x + i, y + j);
	}
	return B;
}

double bff::CotanAngles(MeshKernel::EdgeHandle eh) {
	double cot = 0;
	for (auto& fh : mesh.NeighborFh(eh)) {
		for (auto& vh : mesh.faces(fh).getVertexHandle()) {
			auto vh1 = mesh.edges(eh).vh1(), vh2 = mesh.edges(eh).vh2();
			if (vh != vh1 && vh != vh2) {
				auto v = mesh.vertices(vh), v1 = mesh.vertices(vh1), v2 = mesh.vertices(vh2);
				cot += 1.0 / tan(acos((v - v1) * (v - v2) / (v - v1).norm() / (v - v2).norm()));
			}
		}
	}
	return cot;
}

void bff::Laplace() {
	A = Eigen::MatrixXd::Zero(n, n);
	for (auto& ep : mesh.alledges()) {
		auto& e = ep.second;
		auto eh = ep.first;
		auto vh1 = to_new[e.vh1()], vh2 = to_new[e.vh2()];
		A(vh1, vh2) = A(vh2, vh1) = -0.5 * CotanAngles(eh);
		A(vh1, vh1) += -A(vh1, vh2);
		A(vh2, vh2) += -A(vh1, vh2);
	}
	Aii = submatrix(A, 0, 0, in, in);
	Aib = submatrix(A, 0, in, in, bn);
	Abi = submatrix(A, in, 0, bn, in);
	Abb = submatrix(A, in, in, bn, bn);
}