#pragma once
#include "Kernel/Mesh.h"

using namespace std;

namespace HexAddCell {

	void add_cell_2fh(MeshKernel::VolumeMesh& mesh, FH fh1, FH fh2);// ͨ���ִ��������������

	void add_cell_3fh(MeshKernel::VolumeMesh& mesh, FH fh1, FH fh2, FH fh3);// ͨ���ִ��������������, ����������Ծ�����

	void add_cell_eh_fh(MeshKernel::VolumeMesh& mesh, EH eh, FH fh);// ͨ���ִ��һ��һ��������

	void add_cell_fh(MeshKernel::VolumeMesh& mesh, FH fh);// ͨ���ִ��һ����������

};
