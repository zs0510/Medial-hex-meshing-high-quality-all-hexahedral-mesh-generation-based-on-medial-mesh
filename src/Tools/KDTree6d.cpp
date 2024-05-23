#include "KDTree6d.h"

KDTree6d::KDTree6d(MeshKernel::SurfaceMesh trimesh, unsigned int max_faces, unsigned int max_depth): mesh(trimesh) {

	mesh.genAllVerticesNormal();
	auto fcnt = mesh.FaceSize();
	auto& faces = mesh.allfaces();

	root = new Node();
	root->faces = new Triangles();
	root->fhs = new TrianglesFH();
	root->faces->reserve(fcnt);
	root->fhs->reserve(fcnt);

	mesh.initBBox();
	diagonal_length = (mesh.BBoxMax - mesh.BBoxMin).norm();
	weight_normal = 0.25;// 还在调试中
	double weight_n = diagonal_length * weight_normal;

	for (auto& fp : faces) {
		const auto& f_vhs = fp.second.getVertexHandle();
		vector<Vector6d> triangle;
		for (auto& f_vh : f_vhs) {
			auto& f_v = mesh.vertices(f_vh);
			triangle.emplace_back(Vector6d(f_v.x(), f_v.y(), f_v.z(), f_v.getNormalX() * weight_n, f_v.getNormalY() * weight_n, f_v.getNormalZ() * weight_n));
			//triangle.emplace_back(Vector6d(f_v.x(), f_v.y(), f_v.z(), f_v.getNormalX(), f_v.getNormalY(), f_v.getNormalZ()));
		}
		root->faces->push_back(triangle);
		root->fhs->push_back(fp.first);
	}
	build_recurse(root, max_faces, max_depth);

}

unsigned int KDTree6d::build_recurse(Node* node, unsigned int max_faces, unsigned int depth) {

	if (depth == 0 || (node->faces->size() <= max_faces))
		return depth;

	Vec6 bbox_max(DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
	Vec6 bbox_min(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX);
	for (auto& face : *node->faces) {
		for (auto& v : face) {
			bbox_max[0] = max(bbox_max[0], v.x()), bbox_max[1] = max(bbox_max[1], v.y()), bbox_max[2] = max(bbox_max[2], v.z());
			bbox_min[0] = min(bbox_min[0], v.x()), bbox_min[1] = min(bbox_min[1], v.y()), bbox_min[2] = min(bbox_min[2], v.z());
		}
	}

	// split longest side of bounding box
	Vec6 bb = bbox_max - bbox_min;
	Vec6 bb_center = (bbox_max + bbox_min) / 2;
	double length = bb[0];
	int axis = 0;
	if (bb[1] > length) {
		length = bb[1];
		axis = 1;
	}
	if (bb[2] > length) {
		length = bb[2];
		axis = 2;
	}
	if (bb[3] > length) {
		length = bb[3];
		axis = 3;
	}
	if (bb[4] > length) {
		length = bb[4];
		axis = 4;
	}
	if (bb[5] > length) {
		length = bb[5];
		axis = 5;
	}
	//// split in the middle
	//double split = bb_center[axis];

	// find split position as median
	std::vector<double> axis_pos;
	axis_pos.reserve(node->faces->size() * 3);
	for (auto& face : *(node->faces)) {
		for (auto& v : face) {
			axis_pos.push_back(v[axis]);
		}
	}
	std::sort(axis_pos.begin(), axis_pos.end());
	double split = axis_pos[axis_pos.size() / 2];

	// create children
	auto* left = new Node();
	left->faces = new Triangles();
	left->faces->reserve(node->faces->size() / 2);
	left->fhs = new TrianglesFH();
	left->fhs->reserve(node->faces->size() / 2);

	auto* right = new Node();
	right->faces = new Triangles();
	right->faces->reserve(node->faces->size() / 2);
	right->fhs = new TrianglesFH();
	right->fhs->reserve(node->faces->size() / 2);

	for (int i = 0; i < (*node->faces).size(); ++i) {
		auto& face = (*node->faces)[i];
		auto& fh = (*node->fhs)[i];
		bool l = false, r = false;
		for (auto& v : face) {
			switch (axis) {
			case 0:
				if (v[0] <= split) l = true;
				else r = true;
				break;
			case 1:
				if (v[1] <= split) l = true;
				else r = true;
				break;
			case 2:
				if (v[2] <= split) l = true;
				else r = true;
				break;
			case 3:
				if (v[3] <= split) l = true;
				else r = true;
				break;
			case 4:
				if (v[4] <= split) l = true;
				else r = true;
				break;
			case 5:
				if (v[5] <= split) l = true;
				else r = true;
				break;
			}
		}
		if (l) {
			left->faces->push_back(face);
			left->fhs->push_back(fh);
		} 
		if (r) {
			right->faces->push_back(face);
			right->fhs->push_back(fh);
		} 
	}

	// stop here?
	if (left->faces->size() == node->faces->size() || right->faces->size() == node->faces->size()) {
		node->faces->shrink_to_fit();// compact my memory
		delete left;// delete new nodes
		delete right;
		return depth;// return tree depth
	}
	else {
		// free my memory
		delete node->faces;
		node->faces = nullptr;
		delete node->fhs;
		node->fhs = nullptr;
		// store internal data
		node->axis = axis;
		node->split = split;
		node->left_child = left;
		node->right_child = right;
		// recurse to childen
		int depthLeft = build_recurse(node->left_child, max_faces, depth - 1);
		int depthRight = build_recurse(node->right_child, max_faces, depth - 1);

		return std::min(depthLeft, depthRight);
	}

}

KDTree6d::NearestNeighbor KDTree6d::nearest(Vector6d& v) {

	NearestNeighbor data;
	data.dist = DOUBLE_MAX;
	data.tests = 0;
	nearest_recurse(root, v, data);
	return data;

}

void KDTree6d::nearest_recurse(Node* node, Vector6d& v, NearestNeighbor& data) {


	if (!node->left_child) {// 叶子结点
		double d;
		Vector6d u;
		/*for (auto& face : *node->faces) {
			d = dist_point_triangle(v, face, u);
			data.tests++;
			if (d < data.dist) {
				data.dist = d;
				data.face = face;
				data.nearest = u;
			}
		}*/
		for (int i = 0; i < (*node->faces).size(); ++i) {
			auto& face = (*node->faces)[i];
			auto& fh = (*node->fhs)[i];
			d = dist_point_triangle(v, face, u);
			data.tests++;
			if (d < data.dist) {
				data.dist = d;
				data.face = face;
				data.nearest = u;
				data.fh = fh;
			}
		}
	}
	else {
		double dist = v[node->axis] - node->split;
		if (dist <= 0.0) {
			nearest_recurse(node->left_child, v, data);
			if (fabs(dist) < data.dist)
				nearest_recurse(node->right_child, v, data);
		}
		else {
			nearest_recurse(node->right_child, v, data);
			if (fabs(dist) < data.dist)
				nearest_recurse(node->left_child, v, data);
		}
	}

}

double KDTree6d::dist_point_triangle(Vector6d& v, vector<Vector6d>& face, Vector6d& nearest_vertex) {

	double dist_min = DOUBLE_MAX;
	for (auto& f_v : face) {
		double dist = (f_v - v).Length();
		if (dist < dist_min) {
			dist_min = dist;
			nearest_vertex = f_v;
		}
	}

	return dist_min;

}