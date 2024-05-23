#pragma once
#include"Face.h"

namespace MeshKernel { 
	class iGameCell {
	public:
		/*=========================构造函数=============================*/
		iGameCell() :n_(0) {}
		iGameCell(const std::vector<iGameVertexHandle>& _vertices, const std::vector<iGameEdgeHandle>& _edges,
			 const std::vector<iGameFaceHandle>& _faces){
			//assert(_faces.size() + _vertices.size() - _edges.size() == 2);     //是否满足欧拉公式
			assert(_faces.size() >= 4 && _vertices.size() >= 4);
			n_ = _faces.size();
			vertices_.assign(_vertices.begin(), _vertices.end());
			edges_.assign(_edges.begin(), _edges.end());
			faces_.assign(_faces.begin(), _faces.end());
		}
		iGameCell(const iGameCell& _c) {
			*this = _c;
		}
		/*=======================基本操作符==========================*/
		inline iGameCell& operator=(const iGameCell& _c) {
			n_ = _c.n_;
			vertices_.assign(_c.vertices_.begin(), _c.vertices_.end());
			edges_.assign(_c.edges_.begin(), _c.edges_.end());
			faces_.assign(_c.faces_.begin(), _c.faces_.end());
			return *this;
		}
		bool operator==(const iGameCell& _c) const {
			return (n_ == _c.n_)  && isSameFaces(_c);
		}
		/*=======================读写元素==========================*/
		const iGameVertexHandle vh(int k) const { return vertices_[k]; }          // 常量对象只支持读
		const iGameEdgeHandle eh(int k) const { return edges_[k]; }
		const iGameFaceHandle fh(int k) const { return faces_[k]; }
		
		iGameVertexHandle& vh(int k) { return vertices_[k]; }                     // 非常量对象支持读写
		iGameEdgeHandle& eh(int k) { return edges_[k]; }
		iGameFaceHandle& fh(int k) { return faces_[k]; }
		
		size_t vertices_size() const { return n_; }
		size_t faces_size() const { return faces_.size(); }
		size_t edges_size() const { return edges_.size(); }
		/*================得到有序的顶点handle排列==========================*/
		std::vector<iGameVertexHandle> getSortedVertexHandle() const {
			std::vector<iGameVertexHandle> sortedvertexhandle = vertices_;
			std::sort(sortedvertexhandle.begin(), sortedvertexhandle.end());
			return sortedvertexhandle;
		}
		std::vector<iGameVertexHandle> getVertexHandle() {
			std::vector<iGameVertexHandle> res = vertices_;
			return res;
		}
		std::vector<iGameVertexHandle> getVertexHandle() const {
			std::vector<iGameVertexHandle> res = vertices_;
			return res;
		}
		std::vector<iGameFaceHandle> getFaceHandle() {
			std::vector<iGameFaceHandle> res = faces_;
			return res;
		}
		std::vector<iGameFaceHandle> getFaceHandle() const {
			std::vector<iGameFaceHandle> res = faces_;
			return res;
		}
		std::vector<iGameEdgeHandle> getEdgeHandle() {
			std::vector<iGameEdgeHandle> res = edges_;
			return res;
		}
		std::vector<iGameEdgeHandle> getEdgeHandle() const {
			std::vector<iGameEdgeHandle> res = edges_;
			return res;
		}
		int getN() const {// 得到体上的全部点的个数
			return vertices_.size();
		}
		
	private:
		// 判断输入体的面是否与当前体全部相同
		bool isSameFaces(const iGameCell& _c) const {
			std::unordered_set<iGameFaceHandle> fhset;
			for (const auto& f1 : _c.faces_) fhset.insert(f1);
			for (const auto& f2 : faces_) fhset.insert(f2);
			return fhset.size() == n_;
		}
	private:
		// 保存该体所有点的handle
		std::vector<iGameVertexHandle> vertices_;
		// 保存该体所有边的handle
		std::vector<iGameEdgeHandle> edges_;
		// 保存该体所有面的handle
		std::vector<iGameFaceHandle> faces_;
		// 保存该体的面数
		size_t n_;                     
	};

}

/*======================特化iGameCell的哈希映射============================*/
// To do: 修改Hash映射
namespace std
{
	template<> struct hash<MeshKernel::iGameCell>
	{
		size_t operator()(const MeshKernel::iGameCell& _c)const
		{
			size_t res = 0;
			assert(_c.vertices_size() >= 4);
			auto cv = _c.getSortedVertexHandle();
			for (int i = 0; i < cv.size(); ++i) {
				res ^= hash<int>()(cv[i]);
			}
			return res;
		}
	};
}