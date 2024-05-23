#include "BVH.h"

void BVH_Tree::buildBVH_Tree(std::vector<BVH_Face>& faces) {// 请确保每个面都是三角形

	std::vector<BVH_Node*> total_nodes;
	total_nodes.reserve(faces.size());

	for (auto& face : faces) {

		BVH_Node* node = new BVH_Node();

		node->face = face;

		// 面积
		Vector3d E1 = face.vertices[1] - face.vertices[0];
		Vector3d E2 = face.vertices[2] - face.vertices[0];
		Vector3d C = E1.Cross(E2);
		node->area = C.Length() * 0.5f;

		// 包围盒
		Vector3d _min = face.vertices[0], _max = face.vertices[0];
		for (int i = 1; i < face.vertices.size(); ++i) {
			_min[0] = std::min(_min[0], face.vertices[i][0]);
			_min[1] = std::min(_min[1], face.vertices[i][1]);
			_min[2] = std::min(_min[2], face.vertices[i][2]);
			_max[0] = std::max(_max[0], face.vertices[i][0]);
			_max[1] = std::max(_max[1], face.vertices[i][1]);
			_max[2] = std::max(_max[2], face.vertices[i][2]);
		}
		node->box = Bound3(_min, _max);

		total_nodes.push_back(node);

	}

	std::sort(total_nodes.begin(), total_nodes.end(), [](BVH_Node* node1, BVH_Node* node2) {
		auto& box1 = node1->box;
		auto& box2 = node2->box;
		if (box1.centroid[0] != box2.centroid[0]) return box1.centroid[0] < box2.centroid[0];
		if (box1.centroid[1] != box2.centroid[1]) return box1.centroid[1] < box2.centroid[1];
		return box1.centroid[2] < box2.centroid[2];
		});

	leaves = total_nodes.size();
	//std::cout << "Init leaves success. leaves = " << leaves << std::endl;

	root = new BVH_Node();
	root = recursiveBuild(total_nodes, 2);
	auto& pos_min = root->box.pos_min;
	auto& pos_max = root->box.pos_max;
	/*std::cout << "Init BVH Tree success. height = " << height << ". Bounding Box: (" << pos_min[0] << ", " << pos_min[1] << ", " << pos_min[2] << ") --> " <<
		"(" << pos_max[0] << ", " << pos_max[1] << ", " << pos_max[2] << ")\n";*/

}

void BVH_Tree::buildBVH_Tree(MeshKernel::SurfaceMesh& mesh) {

	if (mesh.fsize() == 0) return;

	std::vector<BVH_Node*> total_nodes;
	
	for (auto& fp : mesh.allfaces()) {
		auto& fh = fp.first;

		BVH_Node* node = new BVH_Node();
		auto& face = mesh.faces(fh);
		auto vhs = face.getVertexHandle();
		std::vector<Vector3d> _vertices;
		for (auto& vh : vhs) {
			auto& v = mesh.vertices(vh);
			_vertices.push_back(Vector3d(v.x(), v.y(), v.z()));
		}

		// 叶子结点的面
		BVH_Face& bvh_face = node->face;
		bvh_face.vertices = _vertices;
		bvh_face.index = fh;

		// 面积
		Vector3d E1 = _vertices[1] - _vertices[0];
		Vector3d E2 = _vertices[2] - _vertices[0];
		Vector3d C = E1.Cross(E2);
		node->area = C.Length() * 0.5f;
		
		// 包围盒
		Vector3d _min = _vertices[0], _max = _vertices[0];
		for (int i = 1; i < _vertices.size(); ++i) {
			_min[0] = std::min(_min[0], _vertices[i][0]);
			_min[1] = std::min(_min[1], _vertices[i][1]);
			_min[2] = std::min(_min[2], _vertices[i][2]);
			_max[0] = std::max(_max[0], _vertices[i][0]);
			_max[1] = std::max(_max[1], _vertices[i][1]);
			_max[2] = std::max(_max[2], _vertices[i][2]);
		}
		node->box = Bound3(_min, _max);

		total_nodes.push_back(node);

	}

	std::sort(total_nodes.begin(), total_nodes.end(), [](BVH_Node* node1, BVH_Node* node2) {
		auto& box1 = node1->box;
		auto& box2 = node2->box;
		if (box1.centroid[0] != box2.centroid[0]) return box1.centroid[0] < box2.centroid[0];
		if (box1.centroid[1] != box2.centroid[1]) return box1.centroid[1] < box2.centroid[1];
		return box1.centroid[2] < box2.centroid[2];
		});

	leaves = total_nodes.size();
	//std::cout << "Init leaves success. leaves = " << leaves << std::endl;
	
	root = new BVH_Node();
	root = recursiveBuild(total_nodes, 2);
	auto& pos_min = root->box.pos_min;
	auto& pos_max = root->box.pos_max;
	/*std::cout << "Init BVH Tree success. height = " << height << ". Bounding Box: (" << pos_min[0] << ", " << pos_min[1] << ", " << pos_min[2] << ") --> " <<
		"(" << pos_max[0] << ", " << pos_max[1] << ", " << pos_max[2] << ")\n";*/

}

void BVH_Tree::buildBVH_Tree(SkeletalMesh& mesh) {

	if (mesh.fsize() == 0) return;

	std::vector<BVH_Node*> total_nodes;

	for (auto& fp : mesh.allfaces()) {
		auto& fh = fp.first;

		BVH_Node* node = new BVH_Node();
		auto& face = mesh.faces(fh);
		auto vhs = face.getVertexHandle();
		std::vector<Vector3d> _vertices;
		for (auto& vh : vhs) {
			auto& v = mesh.vertices(vh);
			_vertices.push_back(Vector3d(v.x(), v.y(), v.z()));
		}

		// 叶子结点的面
		BVH_Face& bvh_face = node->face;
		bvh_face.vertices = _vertices;
		bvh_face.index = fh;

		// 面积
		Vector3d E1 = _vertices[1] - _vertices[0];
		Vector3d E2 = _vertices[2] - _vertices[0];
		Vector3d C = E1.Cross(E2);
		node->area = C.Length() * 0.5f;

		// 包围盒
		Vector3d _min = _vertices[0], _max = _vertices[0];
		for (int i = 1; i < _vertices.size(); ++i) {
			_min[0] = std::min(_min[0], _vertices[i][0]);
			_min[1] = std::min(_min[1], _vertices[i][1]);
			_min[2] = std::min(_min[2], _vertices[i][2]);
			_max[0] = std::max(_max[0], _vertices[i][0]);
			_max[1] = std::max(_max[1], _vertices[i][1]);
			_max[2] = std::max(_max[2], _vertices[i][2]);
		}
		node->box = Bound3(_min, _max);

		total_nodes.push_back(node);

	}

	std::sort(total_nodes.begin(), total_nodes.end(), [](BVH_Node* node1, BVH_Node* node2) {
		auto& box1 = node1->box;
		auto& box2 = node2->box;
		if (box1.centroid[0] != box2.centroid[0]) return box1.centroid[0] < box2.centroid[0];
		if (box1.centroid[1] != box2.centroid[1]) return box1.centroid[1] < box2.centroid[1];
		return box1.centroid[2] < box2.centroid[2];
		});

	leaves = total_nodes.size();
	//std::cout << "Init leaves success. leaves = " << leaves << std::endl;

	root = new BVH_Node();
	root = recursiveBuild(total_nodes, 2);
	auto& pos_min = root->box.pos_min;
	auto& pos_max = root->box.pos_max;
	/*std::cout << "Init BVH Tree success. height = " << height << ". Bounding Box: (" << pos_min[0] << ", " << pos_min[1] << ", " << pos_min[2] << ") --> " <<
		"(" << pos_max[0] << ", " << pos_max[1] << ", " << pos_max[2] << ")\n";*/

}

void BVH_Tree::buildBVH_Tree(MeshKernel::TetMesh& mesh) {

	if (mesh.fsize() == 0) return;

	std::vector<BVH_Node*> total_nodes;

	for (auto& fp : mesh.allfaces()) {
		auto& fh = fp.first;
		if (!mesh.isOnBoundary(fh)) continue;

		BVH_Node* node = new BVH_Node();
		auto& face = mesh.faces(fh);
		auto vhs = face.getVertexHandle();
		std::vector<Vector3d> _vertices;
		for (auto& vh : vhs) {
			auto& v = mesh.vertices(vh);
			_vertices.push_back(Vector3d(v.x(), v.y(), v.z()));
		}

		// 叶子结点的面
		BVH_Face& bvh_face = node->face;
		bvh_face.vertices = _vertices;
		bvh_face.index = fh;

		// 面积
		Vector3d E1 = _vertices[1] - _vertices[0];
		Vector3d E2 = _vertices[2] - _vertices[0];
		Vector3d C = E1.Cross(E2);
		node->area = C.Length() * 0.5f;
		if (vhs.size() == 4) {
			Vector3d E3 = _vertices[3] - _vertices[0];
			C = E2.Cross(E3);
			node->area += C.Length() * 0.5f;
		}


		// 包围盒
		Vector3d _min = _vertices[0], _max = _vertices[0];
		for (int i = 1; i < _vertices.size(); ++i) {
			_min[0] = std::min(_min[0], _vertices[i][0]);
			_min[1] = std::min(_min[1], _vertices[i][1]);
			_min[2] = std::min(_min[2], _vertices[i][2]);
			_max[0] = std::max(_max[0], _vertices[i][0]);
			_max[1] = std::max(_max[1], _vertices[i][1]);
			_max[2] = std::max(_max[2], _vertices[i][2]);
		}
		node->box = Bound3(_min, _max);

		total_nodes.push_back(node);

	}

	std::sort(total_nodes.begin(), total_nodes.end(), [](BVH_Node* node1, BVH_Node* node2) {
		auto& box1 = node1->box;
		auto& box2 = node2->box;
		if (box1.centroid[0] != box2.centroid[0]) return box1.centroid[0] < box2.centroid[0];
		if (box1.centroid[1] != box2.centroid[1]) return box1.centroid[1] < box2.centroid[1];
		return box1.centroid[2] < box2.centroid[2];
		});

	leaves = total_nodes.size();
	std::cout << "Init leaves success. leaves = " << leaves << std::endl;

	root = new BVH_Node();
	root = recursiveBuild(total_nodes, 2);
	auto& pos_min = root->box.pos_min;
	auto& pos_max = root->box.pos_max;
	std::cout << "Init BVH Tree success. height = " << height << ". Bounding Box: (" << pos_min[0] << ", " << pos_min[1] << ", " << pos_min[2] << ") --> " <<
		"(" << pos_max[0] << ", " << pos_max[1] << ", " << pos_max[2] << ")\n";

}

void BVH_Tree::buildBVH_Tree(MeshKernel::VolumeMesh& mesh) {

	if (mesh.fsize() == 0) return;

	std::vector<BVH_Node*> total_nodes;

	for (auto& fp : mesh.allfaces()) {
		auto& fh = fp.first;
		if (!mesh.isOnBoundary(fh)) continue;

		BVH_Node* node = new BVH_Node();
		auto& face = mesh.faces(fh);
		auto vhs = face.getVertexHandle();
		std::vector<Vector3d> _vertices;
		for (auto& vh : vhs) {
			auto& v = mesh.vertices(vh);
			_vertices.push_back(Vector3d(v.x(), v.y(), v.z()));
		}

		// 叶子结点的面
		BVH_Face& bvh_face = node->face;
		bvh_face.vertices = _vertices;
		bvh_face.index = fh;

		// 面积
		Vector3d E1 = _vertices[1] - _vertices[0];
		Vector3d E2 = _vertices[2] - _vertices[0];
		Vector3d C = E1.Cross(E2);
		node->area = C.Length() * 0.5f;
		if (vhs.size() == 4) {
			Vector3d E3 = _vertices[3] - _vertices[0];
			C = E2.Cross(E3);
			node->area += C.Length() * 0.5f;
		}


		// 包围盒
		Vector3d _min = _vertices[0], _max = _vertices[0];
		for (int i = 1; i < _vertices.size(); ++i) {
			_min[0] = std::min(_min[0], _vertices[i][0]);
			_min[1] = std::min(_min[1], _vertices[i][1]);
			_min[2] = std::min(_min[2], _vertices[i][2]);
			_max[0] = std::max(_max[0], _vertices[i][0]);
			_max[1] = std::max(_max[1], _vertices[i][1]);
			_max[2] = std::max(_max[2], _vertices[i][2]);
		}
		node->box = Bound3(_min, _max);

		total_nodes.push_back(node);

	}
	std::sort(total_nodes.begin(), total_nodes.end(), [](BVH_Node* node1, BVH_Node* node2) {
		auto& box1 = node1->box;
		auto& box2 = node2->box;
		if (box1.centroid[0] != box2.centroid[0]) return box1.centroid[0] < box2.centroid[0];
		if (box1.centroid[1] != box2.centroid[1]) return box1.centroid[1] < box2.centroid[1];
		return box1.centroid[2] < box2.centroid[2];
		});

	leaves = total_nodes.size();
	//std::cout << "Init leaves success. leaves = " << leaves << std::endl;

	root = new BVH_Node();
	root = recursiveBuild(total_nodes, 2);
	auto& pos_min = root->box.pos_min;
	auto& pos_max = root->box.pos_max;
	std::cout << "Init BVH Tree success. height = " << height << ". Bounding Box: (" << pos_min[0] << ", " << pos_min[1] << ", " << pos_min[2] << ") --> " <<
		"(" << pos_max[0] << ", " << pos_max[1] << ", " << pos_max[2] << ")\n";

}

BVH_Node* BVH_Tree::recursiveBuild(std::vector<BVH_Node*>& nodes, int h) {

	int sz = nodes.size();
	BVH_Node* res = new BVH_Node();

	if (sz == 1) {
		res = nodes[0];
		if (h > height) height = h;
	} else if (sz == 2) {
		res->left = nodes[0];
		res->right = nodes[1];
		res->box = Union(nodes[0]->box, nodes[1]->box);
		res->area = nodes[0]->area + nodes[1]->area;
		if (++h > height) height = h;
	} else {
		

		auto beginning = nodes.begin();
		auto middling = nodes.begin() + (sz / 2);
		auto ending = nodes.end();
		std::vector<BVH_Node*> left_nodes(beginning, middling);
		std::vector<BVH_Node*> right_nodes(middling, ending);
		assert(sz == left_nodes.size() + right_nodes.size());
		res->left = recursiveBuild(left_nodes, h + 1);
		res->right = recursiveBuild(right_nodes, h + 1);
		res->box = Union(res->left->box, res->right->box);
		res->area = res->left->area + res->right->area;
	}

	return res;

}

void BVH_Tree::destoryBVH_Tree() {

	std::queue<BVH_Node*> que;
	if (root) que.push(root);
	while (!que.empty()) {
		auto* node = que.front();
		que.pop();
		if (node->left) que.push(node->left);
		if (node->right) que.push(node->right);
		delete node;
		node = nullptr;
	}
	
}

Bound3 BVH_Tree::Union(const Bound3& b1, const Bound3& b2) {

	double x_min = std::min(b1.pos_min[0], b2.pos_min[0]);
	double y_min = std::min(b1.pos_min[1], b2.pos_min[1]);
	double z_min = std::min(b1.pos_min[2], b2.pos_min[2]);
	double x_max = std::max(b1.pos_max[0], b2.pos_max[0]);
	double y_max = std::max(b1.pos_max[1], b2.pos_max[1]);
	double z_max = std::max(b1.pos_max[2], b2.pos_max[2]);
	Bound3 box(Vector3d(x_min, y_min, z_min), Vector3d(x_max, y_max, z_max));
	return box;

}

bool BVH_Tree::checkIntersect(const Ray& ray, const Bound3& box) {
	// 若射线与包围盒相交返回 true, 否则返回 false
	double times[3][2];
	times[0][0] = (box.pos_min[0] - ray.pos[0]) * ray.invDir[0];
	times[0][1] = (box.pos_max[0] - ray.pos[0]) * ray.invDir[0];
	times[1][0] = (box.pos_min[1] - ray.pos[1]) * ray.invDir[1];
	times[1][1] = (box.pos_max[1] - ray.pos[1]) * ray.invDir[1];
	times[2][0] = (box.pos_min[2] - ray.pos[2]) * ray.invDir[2];
	times[2][1] = (box.pos_max[2] - ray.pos[2]) * ray.invDir[2];

	// 检查有无负向
	for (int i = 0; i < 3; ++i) {
		if (ray.isNeg[i]) {
			std::swap(times[i][0], times[i][1]);
		}
	}

	double time_enter, time_exit;

	time_enter = std::max(times[0][0], std::max(times[1][0], times[2][0])); //坐标全部进去才算进去
	time_exit = std::min(times[0][1], std::min(times[1][1], times[2][1])); //某一个坐标出去就算出去

	return (time_enter <= time_exit && time_exit >= 0);

}

Intersection BVH_Tree::getIntersection(const Ray& ray, const std::vector<Vector3d>& vertices) {

	// 检查射线是否与三角形相交, 并返回交点
	Intersection isec;

	if (vertices.size() != 3) {
		std::cerr << "Vertices size != 3\n";
		return isec;
	}

	// Refer to https://www.cnblogs.com/graphics/archive/2010/08/09/1795348.html
	double t, u, v;
	const Vector3d& dir = ray.dir;
	Vector3d E1 = vertices[1] - vertices[0];
	Vector3d E2 = vertices[2] - vertices[0];
	// 叉乘
	Vector3d P = dir.Cross(E2);
	// 点乘
	double det = P.Dot(E1);
	Vector3d T;
	if (det > 0) {
		T = ray.pos - vertices[0];
	} else {
		T = vertices[0] - ray.pos;
		det = -det;
	}
	// 未相交
	if (det < 0)
		return isec;

	// Calculate u and make sure u <= 1
	u = T.Dot(P);
	if (u < 0.0f || u > det)
		return isec;

	// Q
	Vector3d Q = T.Cross(E1);

	// Calculate v and make sure u + v <= 1
	v = dir.Dot(Q);
	if (v < 0.0f || u + v > det)
		return isec;

	//// Calculate t, scale parameters, ray intersects triangle
	//t = E2.x() * Q.x() + E2.y() * Q.y() + E2.z() * Q.z();
	float fInvDet = 1.0f / det;
	//t *= fInvDet;
	u *= fInvDet;
	v *= fInvDet;
	isec.pos = u * E1 + v * E2 + vertices[0];
	Vector3d dir1 = isec.pos - ray.pos;
	dir1.Normalize();
	if (dir1.Dot(dir) >= 0.f) {
		isec.happened = true;
		isec.distance = (isec.pos - ray.pos).Length();
		//intersection_count++;
		isec.intersection_count = 1;
	}

	return isec;

}

Intersection BVH_Tree::getIntersection(const Ray& ray, const BVH_Face& face) {

	// 检查射线是否与面相交, 并返回交点
	Intersection isec;
	if (face.vertices.size() == 3) {
		isec = getIntersection(ray, face.vertices);
	} else if (face.vertices.size() == 4) {
		std::vector<Vector3d> triangle0 = { face.vertices[0], face.vertices[1], face.vertices[2] };
		isec = getIntersection(ray, triangle0);
		if (!isec.happened) {
			std::vector<Vector3d> triangle1 = { face.vertices[0], face.vertices[3], face.vertices[2] };
			isec = getIntersection(ray, triangle1);
		}
	}
	isec.index = face.index;
	return isec;

}

Intersection BVH_Tree::getIntersection(const Ray& ray, BVH_Node* node) {

	// 检查射线是否与结点相交, 并返回交点
	Intersection isec;

	// 未相交
	if (node == nullptr || !checkIntersect(ray, node->box)) return isec;

	// 叶子结点
	if (node->left == nullptr && node->right == nullptr) {
		isec = getIntersection(ray, node->face);
		return isec;
	}

	Intersection isec_left, isec_right;
	if (node->left) isec_left = getIntersection(ray, node->left);
	if (node->right) isec_right = getIntersection(ray, node->right);

	isec = (isec_left.distance < isec_right.distance) ? isec_left : isec_right;
	isec.intersection_count = isec_left.intersection_count + isec_right.intersection_count;

	return isec;

}