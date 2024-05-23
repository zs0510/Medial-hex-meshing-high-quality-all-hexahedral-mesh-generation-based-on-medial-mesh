#pragma once
#include "Tools/TinyVector.h"
#include <vector>

// Author: ZS
// Date: 2022/4/25

class ColorBar {
public:
	ColorBar(){
		initColorBar();
	};

	Vector3d getColor(double c);// c should in [0, 1]

private:
	void initColorBar();
	std::vector<std::pair<double, Vector3d>> colorbar;
};

void ColorBar::initColorBar() {

	std::vector<Vector3d> tmp_colors = {
		Vector3d(6, 48, 254), Vector3d(10, 97, 255), Vector3d(12, 147, 255), Vector3d(17, 195, 255), Vector3d(18, 246, 255),
		Vector3d(19, 245, 205), Vector3d(18, 245, 153), Vector3d(17, 246, 103), Vector3d(16, 246, 50), Vector3d(15, 247, 1),
		Vector3d(63, 249, 0), Vector3d(108, 250, 0), Vector3d(156, 253, 0), Vector3d(205, 255, 0), Vector3d(253, 254, 1),
		Vector3d(249, 207, 0), Vector3d(245, 158, 0), Vector3d(243, 108, 0), Vector3d(239, 59, 0), Vector3d(237, 9, 0)
	};
	tmp_colors = {
		Vector3d(0, 0, 255), Vector3d(255, 255, 255), Vector3d(255, 0, 0)
	};
	int sz = tmp_colors.size();
	double step = 1.f / (sz - 1), value = 0.f;
	colorbar.clear();
	colorbar.reserve(sz);
	for (int i = 0; i < sz; ++i) {
		tmp_colors[i] /= 255.f;
		colorbar.emplace_back(value, tmp_colors[i]);
		value += step;
	}

}

Vector3d ColorBar::getColor(double value) {

	if (value < 0.f) value = 0.f;
	if (value > 1.f) value = 1.f;
	Vector3d color(0, 0, 0);
	for (int i = 1; i < colorbar.size(); ++i) {
		if (value <= colorbar[i].first || i == colorbar.size() - 1) {
			double pre_value = colorbar[i - 1].first;
			double next_value = colorbar[i].first;
			Vector3d pre_color = colorbar[i - 1].second;
			Vector3d next_color = colorbar[i].second;
			color = pre_color * (next_value - value) + next_color * (value - pre_value);
			color /= (next_value - pre_value);
			break;
		}
	}
	return color;


}