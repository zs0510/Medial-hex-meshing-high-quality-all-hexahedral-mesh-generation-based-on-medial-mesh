#pragma once
#include "Cell.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <queue>

#ifndef IGMAE_KERNEL_SIMPLIFY
#define IGMAE_KERNEL_SIMPLIFY
typedef MeshKernel::iGameVertex Vex;
typedef MeshKernel::iGameVertex Vec;
typedef MeshKernel::iGameVertexHandle VH;
typedef MeshKernel::iGameEdgeHandle EH;
typedef MeshKernel::iGameFaceHandle FH;
typedef MeshKernel::iGameCellHandle CH;
#endif

namespace MeshKernel { 
	// ����ĸ���
	class Mesh {
	public:
		iGameVertex BBoxMin, BBoxMax;
		

		void initBBox();
		inline bool isValid(iGameVertexHandle _vh) { return vertices_.count(_vh); }
		inline bool isValid(iGameEdgeHandle _eh) { return edges_.count(_eh); }
		inline bool isValid(iGameFaceHandle _fh) { return faces_.count(_fh); }


		/*=========================��дԪ��===============================*/
		// ��ȡIDΪi�Ķ���
		iGameVertex& vertices(iGameVertexHandle _vh);
		const iGameVertex vertices(iGameVertexHandle _vh) const;        // unordered_map �� [] ���������ǳ�����Ա�������޷��Գ�������ʹ��   
		// ��ȡIDΪi�ı�
		iGameEdge& edges(iGameEdgeHandle _eh);
		const iGameEdge& edges(iGameEdgeHandle _eh) const;
		// ��ȡIDΪi����
		iGameFace& faces(iGameFaceHandle _fh);
		const iGameFace faces(iGameFaceHandle _fh) const;

		size_t vsize() const { return vertices_.size(); }
		size_t esize() const { return edges_.size(); }
		size_t fsize() const { return faces_.size(); }
		const std::unordered_map<iGameVertexHandle, iGameVertex>& allvertices() const { return vertices_; }// ֻ�ɶ�
		const std::unordered_map<iGameEdgeHandle, iGameEdge>& alledges() const { return edges_; }
		const std::unordered_map<iGameFaceHandle, iGameFace>& allfaces() const { return faces_; }
		std::unordered_map<iGameVertexHandle, iGameVertex>& allvertices() { return vertices_; }// �ɶ���д
		std::unordered_map<iGameEdgeHandle, iGameEdge>& alledges() { return edges_; }
		std::unordered_map<iGameFaceHandle, iGameFace>& allfaces() { return faces_; }

		/*====================����Ԫ�صõ���ӦID=========================*/
		const iGameVertexHandle vertexhanle(iGameVertex _vertex) const;
		const iGameEdgeHandle edgehandle(iGameEdge& _edge) const;
		const iGameFaceHandle facehandle(iGameFace& _face) const;
		/*======================�õ��ڽӹ�ϵ============================*/
		// ������ڽӵ�
		std::unordered_set<iGameVertexHandle> NeighborVh(iGameVertexHandle _vh);
		// ������ڽӱ�
		std::unordered_set<iGameEdgeHandle>& NeighborEh(iGameVertexHandle _vh);
		// ������ڽ���
		std::unordered_set<iGameFaceHandle>& NeighborFh(iGameVertexHandle _vh);
		// �ߵ��ڽӱ�
		std::unordered_set<iGameEdgeHandle> NeighborEh(iGameEdgeHandle _eh);
		// �ߵ��ڽ���
		std::unordered_set<iGameFaceHandle>& NeighborFh(iGameEdgeHandle _eh);
		// ����ڽ���
		std::unordered_set<iGameFaceHandle> NeighborFh(iGameFaceHandle _fh);// share common edge
		std::unordered_set<iGameFaceHandle> Neighbor2Fh(iGameFaceHandle _fh);// share common vertex
		/*=========================���Ԫ��=============================*/
		iGameVertexHandle AddVertex(const iGameVertex& _v);
		iGameEdgeHandle AddEdge(const iGameVertexHandle& _vh1, const iGameVertexHandle& _vh2);
		iGameFaceHandle AddFace(const std::vector<iGameVertexHandle>& _vhs);

		/*=========================ɾ��Ԫ��=============================*/
		// ɾ���ͼ�Ԫ��ʱɾ��һ���ڽӵĸ߼�Ԫ��
		// ɾ���߼�Ԫ��ʱ�����Ա�����Ԫ��ʹ�õĵͼ�Ԫ��
		iGameVertexHandle DeleteVertex(iGameVertexHandle _vh);
		iGameEdgeHandle DeleteEdge(iGameEdgeHandle _eh);
		iGameFaceHandle DeleteFace(iGameFaceHandle _fh);


		/*=========================����ΨһID========================*/
		iGameVertexHandle GenVertexHandle() { return (iGameVertexHandle)VertexHandleID_++; }
		iGameEdgeHandle GenEdgeHandle() { return (iGameEdgeHandle)EdgeHandleID_++; }
		iGameFaceHandle GenFaceHandle() { return (iGameFaceHandle)FaceHandleID_++; }

		size_t VertexSize() { return vertices_.size(); }
		size_t EdgeSize() {return edges_.size(); }
		size_t FaceSize() { return faces_.size(); }

		/*=========================һЩ���õĻ�������========================*/

		iGameVertex getFaceCenter(iGameFaceHandle fh);
		iGameVertex getEdgeMidpoint(iGameEdgeHandle eh);

		void genAllEdgesLength();
		void genLength(iGameEdgeHandle);
		double getLength(iGameEdgeHandle);

		Vec getFaceNormal(iGameFaceHandle);

		// ������������Ƿ��ڽ�
		bool isConnected(iGameFaceHandle fh1, iGameFaceHandle fh2);
		bool isConnected(iGameEdgeHandle eh1, iGameEdgeHandle eh2);
		bool isConnected(iGameVertexHandle vh1, iGameVertexHandle vh2);

		// �����������֮���eh
		iGameEdgeHandle getEdgeHandle(iGameVertexHandle vh1, iGameVertexHandle vh2);


	protected:
		/*============��һ����ʹ�õ�ID=========*/
		int VertexHandleID_ = 0;
		int EdgeHandleID_ = 0;
		int FaceHandleID_ = 0;
		/*============����������Ϊ���õĿշ���=========*/
		std::unordered_set<iGameVertexHandle> empty_vhs;
		std::unordered_set<iGameEdgeHandle> empty_ehs;
		std::unordered_set<iGameFaceHandle> empty_fhs;

		/*=============�޸��ڽӹ�ϵ============*/
		// ��������������������ĵ�ͱߵ���������
		void AddFace2Neighbor(const iGameFaceHandle& _fh);
		// ���ñ���������������ĵ����ڱ���
		void AddEdge2Neighbor(const iGameEdgeHandle& _eh);
		// Todo: Delete Neighbor
		void DeleteFace2Neighbor(const iGameFaceHandle& _fh);
		void DeleteEdge2Neighbor(const iGameEdgeHandle& _eh);
	protected:
		// handle��Ԫ�صĶ�Ӧ
		std::unordered_map<iGameVertexHandle, iGameVertex> vertices_;
		std::unordered_map<iGameEdgeHandle, iGameEdge> edges_;
		std::unordered_map<iGameFaceHandle, iGameFace> faces_;

		// Ԫ�ص�handle�Ķ�Ӧ
		std::unordered_map<iGameVertex, iGameVertexHandle> Vertex2Vh_;
		std::unordered_map<iGameEdge, iGameEdgeHandle> Edge2Eh_;
		std::unordered_map<iGameFace, iGameFaceHandle> Face2Fh_;

		// �ڽӹ�ϵ
		std::unordered_map<iGameVertexHandle, std::unordered_set<iGameEdgeHandle>> NeighborEhOfVertex_;          //����ڽӱ�
		std::unordered_map<iGameVertexHandle, std::unordered_set<iGameFaceHandle>> NeighborFhOfVertex_;          //����ڽ���
		std::unordered_map<iGameEdgeHandle, std::unordered_set<iGameFaceHandle>> NeighborFhOfEdge_;              //�ߵ��ڽ���
	protected:
		virtual void InitMesh(const std::vector<iGameVertex>& _vertices,
			const std::vector<std::vector<iGameVertexHandle>>& _elements) = 0;                              // ������д
		Mesh& operator=(const Mesh& _mesh);
	};

	// ��������
	class SurfaceMesh : public Mesh {
	public:
		/*==========================���캯��==============================*/
		SurfaceMesh() {};
		SurfaceMesh(const std::vector<iGameVertex>& _vertices, const std::vector<std::vector<iGameVertexHandle>>& _faces) {
			InitMesh(_vertices, _faces);
		}
		//SurfaceMesh(const SurfaceMesh& _surfacemesh);
		/*=============��ʼ������=============*/
		virtual void InitMesh(const std::vector<iGameVertex>& _vertices,
			const std::vector<std::vector<iGameVertexHandle>>& _elements) override;
		SurfaceMesh& operator=(const SurfaceMesh& _surfacemesh);

		bool isOnBoundary(iGameEdgeHandle);
		bool isOnBoundary(iGameVertexHandle);
		bool isOnBoundary(iGameFaceHandle);

		size_t valence(iGameVertexHandle vh) { return NeighborEhOfVertex_[vh].size(); }

		bool hasLoopBoundary();

		void genNormal(iGameFaceHandle);
		void genNormal(iGameVertexHandle);
		void genAllFacesNormal();
		void genAllVerticesNormal();// ��Ȼ������������ķ�����

		void updateAllHandles();
		void eraseComplicatedEdges();// ʹ��̰�Ĳ���ɾȥ�ظ��ߣ�ʹ�����������ݽṹ
		void flipAllFaces();// ��ת���е���
		std::vector<VH> getOrderedBoundaryVH();

		bool isTriangleMesh();
		size_t getBoundaryVerticesCount();

		void destory();// ����������ݲ����ø��� handle
	};

	// ������
	class VolumeMesh : public Mesh {
	public:
		/*==========================���캯��==============================*/
		VolumeMesh() {};
		VolumeMesh(const std::vector<iGameVertex>& _vertices, const std::vector<std::vector<iGameVertexHandle>>& _cells) {
			InitMesh(_vertices, _cells);
		}
		
		/*=============��ʼ������=============*/
		virtual void InitMesh(const std::vector<iGameVertex>& _vertices,
			const std::vector<std::vector<iGameVertexHandle>>& _elements) override;
		VolumeMesh& operator=(const VolumeMesh& _volumemesh);

	public:
		iGameCell& cells(iGameCellHandle _ch);
		const iGameCell cells(iGameCellHandle _ch) const;
		inline bool isValid(iGameVertexHandle _vh) { return vertices_.count(_vh); }
		inline bool isValid(iGameEdgeHandle _eh) { return edges_.count(_eh); }
		inline bool isValid(iGameFaceHandle _fh) { return faces_.count(_fh); }
		inline bool isValid(iGameCellHandle _ch) { return cells_.count(_ch); }

		size_t csize() const { return cells_.size(); }
		size_t CellSize() const { return cells_.size(); }
		const std::unordered_map <iGameCellHandle, iGameCell> & allcells() const { return cells_; }
		
		const iGameCellHandle cellhandle(iGameCell& _cell) const;

		std::unordered_set<iGameCellHandle> NeighborCh(iGameVertexHandle _vh);// ���ڽ���
		std::unordered_set<iGameCellHandle> NeighborCh(iGameEdgeHandle _eh);// ���ڽ���
		std::unordered_set<iGameCellHandle> NeighborCh(iGameFaceHandle _fh);// ���ڽ���
		std::unordered_set<iGameCellHandle> NeighborCh(iGameCellHandle _ch);// ���ڽ���

		iGameCellHandle AddCell(const std::vector<iGameVertexHandle>& _vhs);// ����
		iGameCellHandle AddCell(const std::vector< std::vector<iGameVertexHandle>>& _vhs);

		iGameVertexHandle DeleteVertex(const iGameVertexHandle& _vh);// ɾ��
		iGameEdgeHandle DeleteEdge(const iGameEdgeHandle& _eh);// ɾ��
		iGameFaceHandle DeleteFace(const iGameFaceHandle& _fh);// ɾ��
		iGameCellHandle DeleteCell(const iGameCellHandle& _ch);// ɾ��

		iGameCellHandle GenCellHandle() { return (iGameCellHandle)CellHandleID_++; }

		void updateAllHandles();
		void destory();

		bool isOnBoundary(iGameCellHandle ch);
		bool isOnBoundary(iGameFaceHandle fh);
		bool isOnBoundary(iGameEdgeHandle eh);
		bool isOnBoundary(iGameVertexHandle vh);

		void generateAllBoundaryFaceNormals();// �������ɶ���εķ�����
		void generateAllBoundaryVertexNormals();

		
		iGameVertex getCellCenter(iGameCellHandle ch);

		iGameVertex getQuadNormal(iGameFaceHandle fh);
		double getQuadArea(iGameFaceHandle fh);

		std::vector<std::vector<iGameVertexHandle>> surface_faces_;// LEE
	protected:
		int CellHandleID_ = 0;
		std::unordered_set<iGameCellHandle> empty_chs;
		
		void AddCell2Neighbor(const iGameCellHandle& _ch);
		void DeleteCell2Neighbor(const iGameCellHandle& _ch);

		std::unordered_map<iGameCellHandle, iGameCell> cells_;
		std::unordered_map<iGameCell, iGameCellHandle> Cell2Ch_;

		std::unordered_map<iGameVertexHandle, std::unordered_set<iGameCellHandle>> NeighborChOfVertex_;          //����ڽ���
		std::unordered_map<iGameEdgeHandle, std::unordered_set<iGameCellHandle>> NeighborChOfEdge_;              //�ߵ��ڽ���
		std::unordered_map<iGameFaceHandle, std::unordered_set<iGameCellHandle>> NeighborChOfFace_;              //����ڽ���
	};

	class TriMesh : SurfaceMesh{};
	class QuadMesh : SurfaceMesh{};

	class TetMesh : public VolumeMesh {
	public:
		TetMesh(const std::vector<iGameVertex>& _vertices, const std::vector<std::vector<iGameVertexHandle>>& _cells) {
			InitMesh(_vertices, _cells);
			/*InitHedra(_vertices, _cells, _surface_faces);*/
			printf("init tetrahedron success\n");
		}
		TetMesh(const std::vector<iGameVertex>& _vertices, std::vector<std::vector<std::vector<iGameVertexHandle>>>& _cells,
			const std::vector<std::vector<iGameVertexHandle>>& _surface_faces) {
			InitHedra(_vertices, _cells, _surface_faces);
			printf("init tetrahedron success\n");
		}

		void InitHedra(const std::vector<iGameVertex>& _vertices, std::vector<std::vector<std::vector<iGameVertexHandle>>>& _cells,
			const std::vector<std::vector<iGameVertexHandle>>& _surface_faces);
		void genNormal(iGameFaceHandle);
		void genNormal(iGameVertexHandle);
		void genAllFacesNormal();
		void genAllVerticesNormal();
		
	};


	class HexMesh : VolumeMesh{};
}