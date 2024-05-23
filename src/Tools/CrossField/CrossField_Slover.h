#pragma once
#include "basic_types.h"
#include "cross_field.h"
#include "basic_utils.h"
#include "Tools/TinyVector.h"


class CrossField_Slover {
public:
	CrossField_Slover(){};

	void calc_cross_field(std::vector<Vec3>& points, std::vector<ID3>& triangles, std::vector<ID2>& feature_edges, 
		std::vector<ID3>& singularities, std::vector<std::array<double, 9>>& global_triangle_dir) {
		std::vector<Vec3> edge_dir;
		BuildBackgroundMeshAndGuidingField(points,
			triangles,
			feature_edges,
			edge_dir,
			global_triangle_dir,
			singularities);

		std::cout << "#E dir = " << edge_dir.size() << ", #T dir = " << global_triangle_dir.size() << std::endl;
	}

private:
	//std::vector<Vec3> points;
	//std::vector<ID3> triangles;
	//std::vector<ID2> lines;// 特征边，可为空
	//std::vector<Vec3> edge_dir;//边上的标架分量之一
	//std::vector<ID3> singularities;
	//std::vector<std::array<double, 9>> global_triangle_dir;//三角形的三个点的标架分量之一

};