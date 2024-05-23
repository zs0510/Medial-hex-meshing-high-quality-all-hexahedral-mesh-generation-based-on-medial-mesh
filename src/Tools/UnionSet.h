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
	int get_sz(int a);// ���� a ������ͨ�����ĵ�Ԫ����

	bool merge(int a, int b);

};