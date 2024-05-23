#include "UnionSet.h"

UnionSet::UnionSet(int num_of_element) {
	n = num_of_element;
	sz.resize(n, 1);
	parent.resize(n);
	for (int i = 0; i < n; ++i) {
		parent[i] = i;
	}

}

int UnionSet::find(int a) {
	if (parent[a] == a) {
		return a;
	}
	return parent[a] = find(parent[a]);
}

bool UnionSet::merge(int a, int b) {
	int pa = find(a);
	int pb = find(b);
	if (pa == pb) return false;
	parent[pb] = parent[b] = pa;
	sz[pa] += sz[pb];
	return true;
}

int UnionSet::get_sz(int a) {
	int pa = find(a);
	return sz[pa];
}