#pragma once
#include <vector>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cassert>
#include "Kernel.h"
#include "Handle.h"
#include "Tools/TinyVector.h"

namespace MeshKernel{  
	// 网格中的顶点
	class iGameVertex {
	public:
		iGameVertex(){ }
		iGameVertex(const iGameVertex& _v) :position(_v.position) { }
		iGameVertex(double x, double y, double z){
			position = Vector3d(x, y, z);
		}

		/*===========读写单一分量元素=================*/
		// 常量成员只能读元素
		const double x() const { return position[0]; }
		const double y() const { return position[1]; }
		const double z() const { return position[2]; }

		// 普通成员可以读写元素
		double& x() { return position[0]; }
		double& y() { return position[1]; }
		double& z() { return position[2]; }

		/*============基本运算=================*/
		inline iGameVertex operator+(const iGameVertex& _rhs) const {                             //加法
			return iGameVertex(x() + _rhs.x(), y() + _rhs.y(), z() + _rhs.z());
		}
		inline iGameVertex operator-(const iGameVertex& _rhs) const {                             //减法
			return iGameVertex(x() - _rhs.x(), y() - _rhs.y(), z() - _rhs.z());
		}
		inline iGameVertex operator*(double k) const {                                        //数乘
			return iGameVertex(x() * k, y() * k, z() * k);
		}
		inline iGameVertex operator/(double k) const {                                        //数除
			return iGameVertex(x() / k, y() / k, z() / k);
		}
		inline double operator*(const iGameVertex& _rhs) const {                              //点乘
			return double(x() * _rhs.x() + y() * _rhs.y() + z() * _rhs.z());
		}
		inline double dot(const iGameVertex& _rhs) const {
			return double(x() * _rhs.x() + y() * _rhs.y() + z() * _rhs.z());
		}
		inline iGameVertex operator%(const iGameVertex& _rhs) const {                              //叉乘
			return iGameVertex(y() * _rhs.z() - z() * _rhs.y(),
				z() * _rhs.x() - x() * _rhs.z(),
				x() * _rhs.y() - y() * _rhs.x());
		}
		inline iGameVertex cross(const iGameVertex& _rhs) const  {
			return iGameVertex(y() * _rhs.z() - z() * _rhs.y(),
				z() * _rhs.x() - x() * _rhs.z(),
				x() * _rhs.y() - y() * _rhs.x());
		}

		inline iGameVertex& operator+=(const iGameVertex& _rhs) {
			position[0] += _rhs.x();
			position[1] += _rhs.y();
			position[2] += _rhs.z();
			return *this;
		}
		inline iGameVertex& operator-=(const iGameVertex& _rhs) {
			position[0] -= _rhs.x();
			position[1] -= _rhs.y();
			position[2] -= _rhs.z();
			return *this;
		}
		inline double& operator[](int i) {
			assert(i >= 0 && i < 3);
			return position[i];
		}
		inline iGameVertex& operator*=(double k) {
			position[0] *= k;
			position[1] *= k;
			position[2] *= k;
			return *this;
		}
		inline iGameVertex& operator/=(double k) {
			assert(k != 0.f);
			position[0] /= k;
			position[1] /= k;
			position[2] /= k;
			return *this;
		}
		inline double norm() const {                                                      //模长
			return sqrt(x() * x() + y() * y() + z() * z());
		}
		inline double norm2() const {                                                      //模长方
			return x() * x() + y() * y() + z() * z();
		}
		inline iGameVertex normalized() const {                                                 //返回单位化的自己(自己并不单位化)
			auto m = this->norm();
			if (m < 1e-10) {
				std::cerr << "Divide zero!!!\n";
				return *this;
			}
			return iGameVertex(x() / m, y() / m, z() / m);
		}
		inline iGameVertex& normalize() {                                                 //将自己单位化
			*this = this->normalized();
			return *this;
		}

		inline double dist(const iGameVertex& _rhs) const {                                   //计算俩个向量的距离
			return (*this - _rhs).norm();
		}
		inline void setPosition(double x, double y, double z) {
			position = Vector3d(x, y, z);
		}
		inline void setPosition(Vector3d P) {
			position = P;
		}
		inline void setPosition(std::vector<double> Pos) {
			assert(Pos.size() == 3);
			position = Vector3d(Pos[0], Pos[1], Pos[2]);
		}
		inline void setX(double x) {
			position[0] = x;
		}
		inline void setY(double y) {
			position[1] = y;
		}
		inline void setZ(double z) {
			position[2] = z;
		}
		inline void setNormal(double x, double y, double z) {
			normal = Vector3d(x, y, z);
		}
		inline void setNormal(Vector3d N) {
			normal = N;
		}
		inline void setNormal(std::vector<double> N) {
			if (N.size() != 3) return;
			normal = Vector3d(N[0], N[1], N[2]);
		}
		inline Vector3d getNormal() {
			return normal;
		}
		inline double getNormalX() {
			return normal[0];
		}
		inline double getNormalY() {
			return normal[1];
		}
		inline double getNormalZ() {
			return normal[2];
		}
		
		/*===========比较运算=================*/
		inline bool operator==(const iGameVertex& _rhs) const {
			return (x() == _rhs.x() && y() == _rhs.y() && z() == _rhs.z());
		}
		/*===========赋值运算=================*/
		//返回值为引用，支持连等
		iGameVertex& operator=(const iGameVertex& _v) {
			normal = _v.normal;
			position = _v.position; 
			return *this;
		}
	private:
		Vector3d position;
		Vector3d normal;
	};
}

/*======================特化V3f的哈希映射============================*/
//冲突时会调用相等函数
namespace std
{
	template<> struct hash<MeshKernel::iGameVertex>
	{
		size_t operator()(const MeshKernel::iGameVertex& v)const
		{
			return hash<double>()(v.x()) ^ hash<double>()(v.y()) ^ hash<double>()(v.z());
		}
	};
}