#include "Math_PCA.h"


void Math_PCA::get_principal_components(std::vector<Eigen::VectorXd>& vectors, 
	std::vector<Eigen::VectorXd>& result) {
	
	if (vectors.empty()) return;
	result.clear();

	int n = vectors.size(), m = vectors.front().size();

	Eigen::MatrixXd X(m, n), C;// ���������
	for (int i = 0; i < n; ++i) {// ��ʼ�����ݾ���
		X.col(i) = vectors[i];
	}

	Eigen::MatrixXd vec, val;
	Eigen::MatrixXd pc;
	
	// ���ֵ��
	Eigen::MatrixXd meanval = X.colwise().mean();// ����ÿһά�Ⱦ�ֵ
	Eigen::RowVectorXd meanvecRow = meanval;
	X.rowwise() -= meanvecRow;// ������ֵ��Ϊ0

	// ����Э����: ����Э�������C = XTX / n-1;
	C = X.adjoint() * X;
	C = C.array() / (X.cols() - 1);


	// ��������ֵ����������: ʹ��selfadjont���ն��������㷨ȥ���㣬�����ò�����vec��val������������
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(C);

	vec = eig.eigenvectors();
	val = eig.eigenvalues();

	// ������ʧ��: ȷ������ά��
	int dim = 0;
	double sum = 0;
	for (int i = val.rows() - 1; i >= 0; --i) {
		sum += val(i, 0);
		dim = i;
		if (sum / val.sum() >= 0.95) break;
	}
	dim = val.rows() - dim;

	Eigen::MatrixXd res = X * vec.rightCols(dim);

	result.resize(dim);
	for (int i = 0; i < dim; ++i) {
		result[i] = res.col(i);
	}

	//std::cout << "[PCA]: Input dim = " << n << ", output dim = " << dim << std::endl;

}

void Math_PCA::get_principal_components(std::vector<Eigen::Vector3d>& vectors, 
	std::vector<Eigen::Vector3d>& result) {

	if (vectors.empty()) return;
	result.clear();

	std::vector<Eigen::VectorXd> input_para;
	std::vector<Eigen::VectorXd> output_res;
	int n = vectors.size();
	int m = vectors.front().size();
	for (int i = 0; i < n; ++i) {
		Eigen::VectorXd vxd;
		vxd = vectors[i];
		input_para.emplace_back(vxd);
	}

	get_principal_components(input_para, output_res);
	
	int dim = output_res.size();
	for (int i = 0; i < dim; ++i) {
		auto& vxd_res = output_res[i];
		Eigen::Vector3d v3d(vxd_res.x(), vxd_res.y(), vxd_res.z());
		result.emplace_back(v3d);
	}

}