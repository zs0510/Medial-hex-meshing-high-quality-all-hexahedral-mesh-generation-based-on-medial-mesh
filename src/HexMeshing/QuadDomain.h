#pragma once
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <vector>
#include <QtWidgets/qfiledialog.h>
#include <QtCore/QtCore>
#include <QtWidgets/QLayout>
#include <QtWidgets/QMessageBox>
#include <Eigen/Dense>
#include <TinyVector.h>
#include "Streamline.h"
#include "Project_KD_Tree.h"
#include "Skeletal Mesh/SkeletalMesh.h"
#include <QtWidgets/qinputdialog.h>

struct QuadDomain {
	std::vector<int> singularpoints;
	std::vector<int> streamlines;
	std::vector<int> directions;
};

void read_quad_info(std::vector<Streamline>& streamlines, std::vector<SingularPoint>& points, std::vector<QuadDomain>& quad_domains);// 输出是三维的点