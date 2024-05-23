#pragma once
#include "Handle.h"
#include "Vertex.h"
#include <iostream>


namespace MeshKernel { 
	// 网格中的边
	class iGameEdge {
	public:
		/*=========================构造函数=============================*/
		iGameEdge() :vertices_(2, (iGameVertexHandle)-1) {}
		iGameEdge(const iGameEdge& _e) {
			vertices_ = _e.vertices_;
		}
		iGameEdge(const iGameVertexHandle& _vh1, const iGameVertexHandle& _vh2) :
			vertices_{ _vh1,_vh2 } {
		}

		/*=======================基本操作符==========================*/
		inline iGameEdge& operator=(const iGameEdge& _e) {
			vertices_ = _e.vertices_;
			weight = _e.weight;
			length = _e.length;
			return *this;
		}
		inline bool operator==(const iGameEdge& _e) const {
			return ((vh1() == _e.vh1() && vh2() == _e.vh2())
				|| (vh1() == _e.vh2() && vh2() == _e.vh1()));
		}

		/*=========================读写元素==============================*/
		const iGameVertexHandle vh1() const { return vertices_[0]; }
		const iGameVertexHandle vh2() const { return vertices_[1]; }
		iGameVertexHandle& vh1() { return vertices_[0]; }
		iGameVertexHandle& vh2() { return vertices_[1]; }

		const iGameVertexHandle vh(int k) const { 
			assert(k <= 1 && k >= 0);
			return vertices_[k]; 
		}
		iGameVertexHandle& vh(int k) {
			assert(k <= 1 && k >= 0);
			return vertices_[k];
		}

		/*=================对当前vertices_进行排序=======================*/
		inline void SortVertices() {
			if (vertices_[0] > vertices_[1]) {
				std::swap(vertices_[0], vertices_[1]);
			}
		}

		inline void setWeight(double w) {
			weight = w;
		}

		inline double getWeight() {
			return weight;
		}

		inline void setLength(double l) {
			length = l;
		}

		inline double getLength() {
			return length;
		}

	private:
		// 保存该边两个顶点的vertexhandle
		std::vector<iGameVertexHandle> vertices_;
		double weight = 0.f;
		double length = 0.f;

	};

}

/*======================特化iGameEdge的哈希映射============================*/
// To do: 修改Hash映射
namespace std
{
	template<> struct hash<MeshKernel::iGameEdge>
	{
		size_t operator()(const MeshKernel::iGameEdge& _e)const
		{
			if (_e.vh1() <= _e.vh2()) {
				return hash<int>()(_e.vh1()) ^ hash<int>()(_e.vh2());
			}
			else {
				return hash<int>()(_e.vh2()) ^ hash<int>()(_e.vh1());
			}
			
		}
	};
}