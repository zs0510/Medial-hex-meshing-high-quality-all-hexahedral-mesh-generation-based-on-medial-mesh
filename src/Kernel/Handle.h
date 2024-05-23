#pragma once
#include <iostream>

namespace MeshKernel { 
	class iGameHandle {
	public:
		// ��ʽ���죬������ "iGameHandle h = 1;" ���﷨������ʽת��
		explicit iGameHandle(int _idx) : idx_(_idx) {};
		// ��ֵ����
		iGameHandle& operator=(const iGameHandle& _h) {
			idx_ = _h.idx_;
			return *this;
		}

		// ����handle������int�����޸�
		iGameHandle& operator=(int _idx) {
			idx_ = _idx;
			return *this;
		}

		//�ж�handle�Ƿ����
		inline bool is_valid() const { return idx_ != -1; }
		/*===========================================handle�ıȽϲ���===========================================*/
		inline bool operator<(const iGameHandle& _h) const { return (this->idx_ < _h.idx_); }

		inline bool operator<(int _idx) const { return idx_ < _idx; }

		inline bool operator>(const iGameHandle& _h) const { return (this->idx_ > _h.idx_); }

		inline bool operator>(int _idx) const { return idx_ > _idx; }

		inline bool operator==(const iGameHandle& _h) const { return _h.idx_ == this->idx_; }

		inline bool operator!=(const iGameHandle& _h) const { return _h.idx_ != this->idx_; }

		/*===========================================�޸�������===========================================*/

		inline const int& idx() const { return idx_; }      //ȡԪ��

		void idx(const int& _idx) { idx_ = _idx; }      //�޸�Ԫ��

		inline operator int() const { return idx_; } //��ʽת��Ϊint

		void reset() { idx_ = -1; }  //��ʼ��

	private:
		int idx_;
	};

	//��ʽ�����ܹ������������͵�handle��ֵ����ǰ���͵�handle
	class iGameVertexHandle : public iGameHandle { public: explicit iGameVertexHandle(int _idx = -1) : iGameHandle(_idx) {} };
	class iGameEdgeHandle : public iGameHandle { public: explicit iGameEdgeHandle(int _idx = -1) : iGameHandle(_idx) {} };
	class iGameFaceHandle : public iGameHandle { public: explicit iGameFaceHandle(int _idx = -1) : iGameHandle(_idx) {} };
	class iGameCellHandle : public iGameHandle { public: explicit iGameCellHandle(int _idx = -1) : iGameHandle(_idx) {} };
}

/*======================�ػ�����iGameHandle�Ĺ�ϣӳ��============================*/
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
// �������������û����ѱ����ǵĲ�����
// ���� "iGameHandle& operator=(int _idx)"
//void test() {
//	iGameVertexHandle vh(1);
//	vh = 2;                      ���󣬲�������ʽת��
//	vh.iGameHandle::operator= (2);    ��ȷ�����û��������
//}