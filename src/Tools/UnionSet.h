#pragma once
#include <vector>
#include <algorithm>

using namespace std;

class UnionSet {
public:
	UnionSet(int num_of_element);

private:
	int n;
	vector<int> parent, sz;

	int find(int a);
	int get_sz(int a);// 返回 a 所在连通分量的单元数量

	bool merge(int a, int b);

};