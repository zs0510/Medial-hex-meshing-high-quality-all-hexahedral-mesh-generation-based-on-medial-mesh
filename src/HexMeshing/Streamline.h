#pragma once
#include "Kernel/Mesh.h"

struct MarchNode {
	FH pre_fh;// 上条交线所在的面
	EH eh;// 交点所在的边
	Vec dir;// 到达交点时的方向向量
	Vex pos;// 交点的坐标
};

struct SingularPoint {
	Vec pos;
	int handle = -1, mode = 0, id = -1;// mode: 0 -> vh, 1 -> eh, 2 -> fh
};

struct Streamline {
	int beg_sid = -1, end_sid = -1;// 奇异点 id
	bool end_on_boundary = false;// true 表示流线终止于边界
	int end_edgehandle = -1;
	bool end_on_fault = false;// true 表示流线无效
	std::vector<Vex> points;// 点坐标
	std::vector<FH> fhs;// 经过的面片
	std::vector<double> length;
	std::vector<double> ratio;
	std::vector<int> opposite_slid;
	void init_length_ratio() {
		length = std::vector<double>(points.size(), 0);
		ratio = std::vector<double>(points.size(), 0);
		for (int i = 1; i < points.size(); ++i) {
			length[i] = length[i - 1] + (points[i] - points[i - 1]).norm();
		}
		for (int i = 1; i < points.size(); ++i) {
			ratio[i] = length[i] / length.back();
		}
	}

	void reverse() {
		std::swap(beg_sid, end_sid);
		std::reverse(points.begin(), points.end());
		this->init_length_ratio();
	}

	Vex get_point(double r) {// 输入参数 r in [0, 1]，返回使用该参数插值得到的点
		if (length.empty()) init_length_ratio();
		Vex res(0, 0, 0);
		if (r <= 0) return points[0];
		if (r >= 1) return points.back();
		int nrid = std::lower_bound(ratio.begin(), ratio.end(), r) - ratio.begin();
		//std::cout << "Input r: " << r << ", nr: " << ratio[nrid] << std::endl;
		double diffp = r - ratio[nrid - 1], diffn = ratio[nrid] - r;
		res = points[nrid] * diffp + points[nrid - 1] * diffn;
		res /= (diffp + diffn);
		/*std::cout << "P0: (" << points[nrid - 1].x() << ", " << points[nrid - 1].y() << "), "
			<< "P1: (" << points[nrid].x() << ", " << points[nrid].y() << "), "
			<< "PI: (" << res.x() << ", " << res.y() << ").\n";*/
		return res;
	}


};



