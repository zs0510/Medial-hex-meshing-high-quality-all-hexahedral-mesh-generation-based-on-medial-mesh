#pragma once
#include "Edge.h"
#include <algorithm>
typedef MeshKernel::iGameVertex iGameVector;

namespace MeshKernel{
	// �����е���Ƭ
	class iGameFace {
	public: 
		/*=========================���캯��=============================*/
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
		/*=======================����������==========================*/
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

		/*=======================��дԪ��==========================*/
		const iGameVertexHandle vh(int k) const { return vertices_[k]; }            // ��������ֻ֧�ֶ�
		const iGameEdgeHandle eh(int k) const { return edges_[k]; }
		iGameVertexHandle& vh(int k) { return vertices_[k]; }		               // �ǳ�������֧�ֶ�д
		iGameEdgeHandle& eh(int k) { return edges_[k]; }
		
		const size_t size() const { return n_; }

		/*================�õ�����Ķ���handle����==========================*/
		std::vector<iGameVertexHandle> getSortedVertexHandle() const {
			std::vector<iGameVertexHandle> sortedvertexhandle = vertices_;
			std::sort(sortedvertexhandle.begin(), sortedvertexhandle.end());
			return sortedvertexhandle;
		}

		/*================����handle����==========================*/
		std::vector<iGameVertexHandle> getVertexHandle() const {
			return vertices_;
		}
		
		/*================��handle����==========================*/
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
		// �ж�������ı��Ƿ��뵱ǰ��ȫ����ͬ
		bool isSameEdges(const iGameFace& _f) const {
			std::unordered_set<iGameEdgeHandle> ehset;
			for (const auto& e1 : _f.edges_) ehset.insert(e1);
			for (const auto& e2 : edges_) ehset.insert(e2);
			return ehset.size() == n_;
		}
	private:
		// ����������бߵ�handle
		std::vector<iGameEdgeHandle> edges_;
		// ����������е��handle
		std::vector<iGameVertexHandle> vertices_;
		// �������ı���
		size_t n_;
		// ������
		iGameVector normal;
	};	
}

/*======================�ػ�iGameFace�Ĺ�ϣӳ��============================*/
// To do: �޸�Hashӳ��
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