#pragma once
#include <set>
#include <omp.h>
#include "Kernel/Mesh.h"
#include "Tools/meshmath.h"
#include "Streamline.h"

class QuadMesher {
public:
	// ���룺�����������������
	QuadMesher(MeshKernel::SurfaceMesh _trimesh, std::vector<EH> _feature_ehs): trimesh(_trimesh), feature_ehs(_feature_ehs){};
	
	bool meshing(MeshKernel::SurfaceMesh& quadmesh);

	void get_streamlines(std::vector<Streamline>& res, std::vector<SingularPoint>& points) {// ��������
		res = streamlines;
		points = singular_points;
	}

	void get_cross_field(std::vector<std::vector<Vec>>& _cross_field) {
		_cross_field = crossfield;
	}

	void get_streamlines_bad(std::vector<Streamline>& res) {
		res = bad_streamlines;
	}

private:

#pragma region data

	MeshKernel::SurfaceMesh trimesh;// input parameter
	std::vector<EH> feature_ehs;

	std::vector<SingularPoint> singular_points;// ���е������
	std::vector<Streamline> streamlines;// ������Ч������
	std::vector<std::vector<Vec>> crossfield;// ���ж���ı����Ϣ
	std::vector<EH> boundary_ehs;
	std::unordered_map<EH, int> singular_neighbor_ehs;// eh �� �ñ��������id ��ӳ��
	std::unordered_map<FH, int> singular_neighbor_fhs;// fh �� �����������id ��ӳ��

	std::vector<Streamline> bad_streamlines;// ��Ч�����ߣ�������

	// ���������� handle �����Ӧ����� id ��ӳ��
	std::unordered_map<VH, int> vh_to_spid;
	std::unordered_map<EH, int> eh_to_spid;
	std::unordered_map<FH, int> fh_to_spid;

#pragma endregion

#pragma region pipeline
	bool preprocessing();

	void calculate_cross_field();// �����ܳ� by Su
	void generate_initial_streamlines();// ��������
	void simplify_streamlines();// ���ߵļ�
	void quad_meshing();

#pragma endregion

	void init_singular_neighbor_handles();

	/* �·����������� */
	void lengthen_initial_streamlines(const MarchNode&, Streamline&, std::unordered_map<FH, bool>& visited, int bdy_face_cnt = 0);// ��ʹ���ϸ���Ƭ�ı����Ϣ

	void merge_neighbor_edgehandle();// ��ǰ��Ԥ���������ڵ��յ�ߺϲ�

	bool get_first_point_streamline(Vec v, EH eh, FH fh, std::vector<Vec>& march_dir, std::vector<double>& weight_vh1);// ����Ƿ��д������ v �������� eh ������
	std::vector<Vec> get_cross_dirs(const Vec&, EH eh, FH fh);// ���ظõ��ڱ��ϲ�ֵ�õ��ı��
	std::vector<Vec> get_cross_dirs(VH vh, FH fh);// ���ظõ��ڸ���������Ƭ�ϵı��
	std::vector<Vec> get_cross_dirs(const Vec&);// ͨ�������һ����ܳ����������ĸ������ı��
	bool get_intersection_line2(const Vex& source, const Vec& dir, EH eh, Vex& intersection);// ���ض�ά�������߶� eh �Ľ���, true ��ʾ�ཻ, false ��ʾδ�ཻ
	Vec get_most_similar_dir(const Vec& ref_dir, const std::vector<Vec>& dirs);// ���� dirs ���� ref_dir ������ӽ��ķ�������
	std::unordered_set<EH> get_smallest_bounding_ehs(const EH& eh);// ���������ߵ���С����εı�
	std::vector<std::pair<EH, FH>> get_samllest_bounding_ehs_fhs(const EH& eh);// ���ذ�Χ�����ߵ���С����εı��Լ���
	void get_streamline_itrp(const Streamline& sl0, const Streamline& sl1, Streamline& sl_new);// �� sl0 �� sl1 ��ֵ�õ�������, ע��: sl0 ���� sl1 ͬ�� 
	void get_among_ehs(const EH& beg, const EH& end, std::set<EH>& res);// �������� ��beg �� ��end ֮�����С�����߽�߼���(�� beg �� end)
	void generate_initial_streamlines_v(const SingularPoint& point);// ��������
	void generate_initial_streamlines_e(const SingularPoint& point);// ��������
	void generate_initial_streamlines_f(const SingularPoint& point);// ��������
	
};