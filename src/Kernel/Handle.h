#pragma once
#include <iostream>

namespace MeshKernel { 
	class iGameHandle {
	public:
		// 显式构造，不允许 "iGameHandle h = 1;" 的语法进行隐式转换
		explicit iGameHandle(int _idx) : idx_(_idx) {};
		// 赋值构造
		iGameHandle& operator=(const iGameHandle& _h) {
			idx_ = _h.idx_;
			return *this;
		}

		// 基类handle允许用int进行修改
		iGameHandle& operator=(int _idx) {
			idx_ = _idx;
			return *this;
		}

		//判定handle是否存在
		inline bool is_valid() const { return idx_ != -1; }
		/*===========================================handle的比较操作===========================================*/
		inline bool operator<(const iGameHandle& _h) const { return (this->idx_ < _h.idx_); }

		inline bool operator<(int _idx) const { return idx_ < _idx; }

		inline bool operator>(const iGameHandle& _h) const { return (this->idx_ > _h.idx_); }

		inline bool operator>(int _idx) const { return idx_ > _idx; }

		inline bool operator==(const iGameHandle& _h) const { return _h.idx_ == this->idx_; }

		inline bool operator!=(const iGameHandle& _h) const { return _h.idx_ != this->idx_; }

		/*===========================================修改与重置===========================================*/

		inline const int& idx() const { return idx_; }      //取元素

		void idx(const int& _idx) { idx_ = _idx; }      //修改元素

		inline operator int() const { return idx_; } //隐式转换为int

		void reset() { idx_ = -1; }  //初始化

	private:
		int idx_;
	};

	//显式构造能够避免其他类型的handle赋值给当前类型的handle
	class iGameVertexHandle : public iGameHandle { public: explicit iGameVertexHandle(int _idx = -1) : iGameHandle(_idx) {} };
	class iGameEdgeHandle : public iGameHandle { public: explicit iGameEdgeHandle(int _idx = -1) : iGameHandle(_idx) {} };
	class iGameFaceHandle : public iGameHandle { public: explicit iGameFaceHandle(int _idx = -1) : iGameHandle(_idx) {} };
	class iGameCellHandle : public iGameHandle { public: explicit iGameCellHandle(int _idx = -1) : iGameHandle(_idx) {} };
}

/*======================特化各种iGameHandle的哈希映射============================*/
namespace std
{
	template<>
	struct hash<MeshKernel::iGameVertexHandle>
	{
		size_t operator()(const MeshKernel::iGameVertexHandle& h)const
		{
			return hash<int>()(h.idx());
		}
	};
	template<>
	struct hash<MeshKernel::iGameEdgeHandle>
	{
		size_t operator()(const MeshKernel::iGameEdgeHandle& h)const
		{
			return hash<int>()(h.idx());
		}
	};
	template<>
	struct hash<MeshKernel::iGameFaceHandle>
	{
		size_t operator()(const MeshKernel::iGameFaceHandle& h)const
		{
			return hash<int>()(h.idx());
		}
	};
	template<>
	struct hash<MeshKernel::iGameCellHandle>
	{
		size_t operator()(const MeshKernel::iGameCellHandle& h)const
		{
			return hash<int>()(h.idx());
		}
	};

}






// note:
// 派生类如果想调用基类已被覆盖的操作符
// 例如 "iGameHandle& operator=(int _idx)"
//void test() {
//	iGameVertexHandle vh(1);
//	vh = 2;                      错误，不允许隐式转换
//	vh.iGameHandle::operator= (2);    正确，调用基类操作符
//}