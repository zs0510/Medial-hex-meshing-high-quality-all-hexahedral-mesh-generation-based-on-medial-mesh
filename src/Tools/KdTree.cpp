#include "KdTree.h"
#include <cmath>

using namespace MeshKernel;
using namespace std;

namespace MeshKernel {

	iGameKdTree::iGameKdTree(SurfaceMesh _mesh, unsigned int max_faces, unsigned int max_depth) : mesh(_mesh) {
		root = new Node();
		root->faces = new Triangles();
		auto fcnt = mesh.FaceSize();
		auto& faces = mesh.allfaces();
		root->faces->reserve(fcnt);
		for (auto& fp : faces) {
			//auto& vhs = fp.second.getVertexHandle();
			root->faces->push_back(fp.first);
		}
		build_recurse(root, max_faces, max_depth);
		//std::cout << "Build KD Tree success.\n";
	}

	

	unsigned int iGameKdTree::build_recurse(Node* node, unsigned int max_faces, unsigned int depth) {

		if (depth == 0 || (node->faces->size() <= max_faces))
			return depth;

		Vec3 bbox_max(DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN);
		Vec3 bbox_min(DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX);
		for (auto& fh : *node->faces) {
			auto& face = mesh.faces(fh);
			auto vhs = face.getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = mesh.vertices(vh);
				bbox_max[0] = max(bbox_max[0], v.x()), bbox_max[1] = max(bbox_max[1], v.y()), bbox_max[2] = max(bbox_max[2], v.z());
				bbox_min[0] = min(bbox_min[0], v.x()), bbox_min[1] = min(bbox_min[1], v.y()), bbox_min[2] = min(bbox_min[2], v.z());
			}
		}

		// split longest side of bounding box
		Vec3 bb = bbox_max - bbox_min;
		Vec3 bb_center = (bbox_max + bbox_min) / 2;
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
		//// split in the middle
		//double split = bb_center[axis];

		// find split position as median
		std::vector<double> axis_pos;
		axis_pos.reserve(node->faces->size() * 3);
		for (auto& fh : *(node->faces)) {
			auto& face = mesh.faces(fh);
			auto vhs = face.getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = mesh.vertices(vh);
				axis_pos.push_back(v[axis]);
			}
		}
		std::sort(axis_pos.begin(), axis_pos.end());
		double split = axis_pos[axis_pos.size() / 2];

		// create children
		auto* left = new Node();
		left->faces = new Triangles();
		left->faces->reserve(node->faces->size() / 2);
		auto* right = new Node();
		right->faces = new Triangles();
		right->faces->reserve(node->faces->size() / 2);

		for (auto& fh : *node->faces) {
			auto& face = mesh.faces(fh);
			bool l = false, r = false;
			auto vhs = face.getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = mesh.vertices(vh);
				switch (axis) {
					case 0:
						if (v.x() <= split) l = true;
						else r = true;
						break;
					case 1:
						if (v.y() <= split) l = true;
						else r = true;
						break;
					case 2:
						if (v.z() <= split) l = true;
						else r = true;
						break;
				}
			}
			if (l) left->faces->push_back(fh);
			if (r) right->faces->push_back(fh);
		}

		// stop here?
		if (left->faces->size() == node->faces->size() || right->faces->size() == node->faces->size()) {
			node->faces->shrink_to_fit();// compact my memory
			delete left;// delete new nodes
			delete right;
			return depth;// return tree depth
		} else {
			// free my memory
			delete node->faces;
			node->faces = nullptr;
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

	iGameKdTree::NearestNeighbor iGameKdTree::nearest(iGameVertex& v) {
		NearestNeighbor data;
		data.dist = DOUBLE_MAX;
		data.tests = 0;
		nearest_recurse(root, v, data);
		return data;
	}

	void iGameKdTree::nearest_recurse(Node* node, iGameVertex& v, NearestNeighbor& data) {

		if (!node->left_child) {// р╤вс╫А╣Ц
			double d;
			iGameVertex u;
			for (auto& fh : *node->faces) {
				d = dist_point_triangle(v, fh, u);
				data.tests++;
				if (d < data.dist) {
					data.dist = d;
					data.fh = fh;
					data.nearest = u;
				}
			}
		} else {
			double dist = v[node->axis] - node->split;
			if (dist <= 0.0) {
				nearest_recurse(node->left_child, v, data);
				if (fabs(dist) < data.dist)
					nearest_recurse(node->right_child, v, data);
			} else {
				nearest_recurse(node->right_child, v, data);
				if (fabs(dist) < data.dist)
					nearest_recurse(node->left_child, v, data);
			}
		}
	}

	double iGameKdTree::dist_point_triangle(iGameVertex& v, FH& fh, iGameVertex& nearest_vertex) {
		auto& face = mesh.faces(fh);
		auto vhs = face.getVertexHandle();
		auto& v0 = mesh.vertices(vhs[0]), & v1 = mesh.vertices(vhs[1]), & v2 = mesh.vertices(vhs[2]);
		iGameVertex vec01 = v1 - v0;
		iGameVertex vec02 = v2 - v0;
		iGameVertex n = vec01 % vec02;
		double d = n.norm2();
		// Check if the triangle is degenerated -> measure dist to line segments
		if (fabs(d) < 1E-6F) {
			iGameVertex q, qq;// nearest_vertex
			double d, dd = DOUBLE_MAX;
			dd = dist_point_line_segment(v, v0, v1, qq);
			d = dist_point_line_segment(v, v1, v2, q);
			if (d < dd) {
				dd = d;
				qq = q;
			}
			d = dist_point_line_segment(v, v2, v0, q);
			if (d < dd) {
				dd = d;
				qq = q;
			}
			nearest_vertex = qq;
			return dd;
		}

		double inv_d = 1.f / d;
		iGameVertex vec12 = v2 - v1;
		iGameVertex v0v = v - v0;
		iGameVertex t = v0v % n;
		double a = (t * vec02) * (-inv_d);
		double b = (t * vec01) * (inv_d);
		double s01, s02, s12;
		// TODO
		if (a < 0.f) {
			s02 = vec02 * v0v / vec02.norm2();
			if (s02 < 0.f) {
				s01 = vec01 * v0v / vec01.norm2();
				if (s01 <= 0.f) v0v = v0;
				else if (s01 >= 1.f) v0v = v1;
				else {
					vec01 *= s01;
					v0v = v0 + vec01;
				}
			} else if (s02 > 1.f) {
				s12 = vec12 * (v - v1) / vec12.norm2();
				if (s12 >= 1.f) v0v = v2;
				else if (s12 <= 0.f) v0v = v1;
				else {
					vec12 *= s12;
					v0v = v1 + vec12;
				}
			} else {
				vec02 *= s02;
				v0v = v0 + vec02;
			}
		} else if (b < 0.f) {
			s01 = vec01 * v0v / vec01.norm2();
			if (s01 < 0.f) {
				s02 = vec02 * (v - v1) / vec02.norm2();
				if (s02 <= 0.f) v0v = v0;
				else if (s02 >= 1.f) v0v = v2;
				else {
					vec02 *= s02;
					v0v = v0 + vec02;
				}
			} else if (s01 > 1.f) {
				s12 = vec12 * (v - v1) / vec12.norm2();
				if (s12 >= 1.f) v0v = v2;
				else if (s12 <= 0.f) v0v = v1;
				else {
					vec12 *= s12;
					v0v = v1 + vec12;
				}
			} else {
				vec01 *= s01;
				v0v = v0 + vec01;
			}
		} else if (a + b > 1.f) {
			s12 = vec12 * (v - v1) / vec12.norm2();
			if (s12 >= 1.f) {
				s02 = vec02 * v0v / vec02.norm2();
				if (s02 <= 0.f) v0v = v0;
				else if (s02 >= 1.f) v0v = v2;
				else {
					vec02 *= s02;
					v0v = v0 + vec02;
				}
			} else if (s12 <= 0.f) {
				s01 = vec01 * v0v / vec01.norm2();
				if (s01 <= 0.0) v0v = v0;
				else if (s01 >= 1.f) v0v = v1;
				else {
					vec01 *= s01;
					v0v = v0 + vec01;
				}
			} else {
				vec12 *= s12;
				v0v = v1 + vec12;
			}
		} else {
			n *= ((n * v0v) * inv_d);
			v0v = v - n;
		}
		nearest_vertex = v0v;
		v0v -= v;
		return v0v.norm();
	}

	double iGameKdTree::dist_point_line_segment(iGameVertex& v, iGameVertex& v0, iGameVertex& v1, iGameVertex& nearest_vertex) {
		iGameVertex vec1 = v - v0;
		iGameVertex vec2 = v1 - v0;
		double t = vec2 * vec2;
		iGameVertex min_v = v0;
		if (t > 1E-6F) {
			t = vec1 * vec2 / t;
			if (t > 1.f) {
				min_v = v1;
				vec1 = v - min_v;
			} else if (t > 0.f) {
				min_v = v0 + vec2 * t;
				vec1 = v - min_v;
			}
		}
		nearest_vertex = min_v;
		return vec1.norm();
	}


};

