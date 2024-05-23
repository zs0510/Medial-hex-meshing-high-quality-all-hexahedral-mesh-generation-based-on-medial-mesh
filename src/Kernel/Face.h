#pragma once
#include "Edge.h"
#include <algorithm>
typedef MeshKernel::iGameVertex iGameVector;

namespace MeshKernel{
	// 网格中的面片
	class iGameFace {
	public: 
		/*=========================构造函数=============================*/
		iGameFace() :n_(0) {}
		iGameFace(const std::vector<iGameVertexHandle>& _vertices, const std::vector<iGameEdgeHandle>& _edges)  {
			assert(_edges.size() == _vertices.size());
			assert(_edges.size() >= 3);
			n_ = _vertices.size();
			edges_.assign(_edges.begin(), _edges.end());
			vertices_.assign(_vertices.begin(), _vertices.end());
		}
		iGameFace(const iGameFace& _f) {
			*this = _f;
		}
		/*=======================基本操作符==========================*/
		inline iGameFace& operator=(const iGameFace& _f) {
			n_ = _f.n_;
			normal = _f.normal;
			edges_.assign(_f.edges_.begin(), _f.edges_.end());
			vertices_.assign(_f.vertices_.begin(), _f.vertices_.end());
			return *this;
		}
		bool operator==(const iGameFace& _f) const {
			return (n_== _f.n_) && isSameEdges(_f);
		}

		/*=======================读写元素==========================*/
		const iGameVertexHandle vh(int k) const { return vertices_[k]; }            // 常量对象只支持读
		const iGameEdgeHandle eh(int k) const { return edges_[k]; }
		iGameVertexHandle& vh(int k) { return vertices_[k]; }		               // 非常量对象支持读写
		iGameEdgeHandle& eh(int k) { return edges_[k]; }
		
		const size_t size() const { return n_; }

		/*================得到有序的顶点handle排列==========================*/
		std::vector<iGameVertexHandle> getSortedVertexHandle() const {
			std::vector<iGameVertexHandle> sortedvertexhandle = vertices_;
			std::sort(sortedvertexhandle.begin(), sortedvertexhandle.end());
			return sortedvertexhandle;
		}

		/*================顶点handle排列==========================*/
		std::vector<iGameVertexHandle> getVertexHandle() const {
			return vertices_;
		}
		
		/*================边handle排列==========================*/
		std::vector<iGameEdgeHandle> getEdgeHandle() const {
			return edges_;
		}

		void setNormal(double x, double y, double z) {
			normal = iGameVector(x, y, z);
		}

		iGameVector getNormal() {
			return normal;
		}

		inline double getNormalX() {
			return normal.x();
		}

		inline double getNormalY() {
			return normal.y();
		}

		inline double getNormalZ() {
			return normal.z();
		}


		
	private:
		// 判断输入面的边是否与当前面全部相同
		bool isSameEdges(const iGameFace& _f) const {
			std::unordered_set<iGameEdgeHandle> ehset;
			for (const auto& e1 : _f.edges_) ehset.insert(e1);
			for (const auto& e2 : edges_) ehset.insert(e2);
			return ehset.size() == n_;
		}
	private:
		// 保存该面所有边的handle
		std::vector<iGameEdgeHandle> edges_;
		// 保存该面所有点的handle
		std::vector<iGameVertexHandle> vertices_;
		// 保存该面的边数
		size_t n_;
		// 法向量
		iGameVector normal;
	};	
}

/*======================特化iGameFace的哈希映射============================*/
// To do: 修改Hash映射
namespace std
{
	template<> struct hash <MeshKernel::iGameFace>
	{
		size_t operator()(const MeshKernel::iGameFace& _f)const
		{
			size_t res = 0;
			assert(_f.size() >= 3);
			auto fv = _f.getSortedVertexHandle();
			for (int i = 0; i < fv.size(); ++i) {
				res ^= hash<int>()(fv[i]);
			}
			return res;
		}
	};
}