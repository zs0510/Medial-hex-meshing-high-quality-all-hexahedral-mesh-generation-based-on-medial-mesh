#pragma once
#include <vector>
#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Core>

// Principal Components Analysis
// Created by Zhang Sheng at 2022/7/25

namespace Math_PCA {

	void get_principal_components(std::vector<Eigen::VectorXd>& vectors, std::vector<Eigen::VectorXd>& result);

	void get_principal_components(std::vector<Eigen::Vector3d>& vectors, std::vector<Eigen::Vector3d>& result);

};