#pragma once

namespace MeshKernel {

	// 基本元素在网格中的唯一ID 
	class iGameHandle;
	class iGameVertexHandle;                  // 顶点的ID 
	class iGameEdgeHandle;                    // 边的ID
	class iGameFaceHandle;                    // 面单元的ID
	class iGameCellHandle;                    // 体单元的ID

	// 网格基本元素
	class iGameVertex;                        // 顶点
	class iGameEdge;                          // 边
	class iGameFace;                          // 面
	class iGameCell;                          // 体

	// 不同类型的网格
	class Mesh;                          // 基础的网格概念
	class SurfaceMesh;                   // 表面网格
	class TriMesh;                       // 三角形网格
	class QuadMesh;                      // 四边形网格
	class VolumeMesh;                    // 体网格
	class TetMesh;                       // 四面体网格
	class HexMesh;                       // 六面体网格


}




