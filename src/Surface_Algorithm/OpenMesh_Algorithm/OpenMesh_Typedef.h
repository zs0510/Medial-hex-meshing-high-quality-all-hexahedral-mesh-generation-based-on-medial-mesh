#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/Property.hh>

// iGame-Lab's OpenMesh setting
struct iGameTraits : public OpenMesh::DefaultTraits {
	// change float to double
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
};

typedef OpenMesh::TriMesh_ArrayKernelT<iGameTraits> iGame_OpenMesh_TriMesh;
typedef OpenMesh::TriMesh_ArrayKernelT<iGameTraits> OM_TriMesh;
typedef iGame_OpenMesh_TriMesh::Point OM_Vex;
typedef iGame_OpenMesh_TriMesh::Point OM_Vec;
typedef iGame_OpenMesh_TriMesh::Point OM_Pos;
typedef iGame_OpenMesh_TriMesh::VertexHandle OM_VH;
typedef iGame_OpenMesh_TriMesh::EdgeHandle OM_EH;
typedef iGame_OpenMesh_TriMesh::FaceHandle OM_FH;
typedef iGame_OpenMesh_TriMesh::HalfedgeHandle OM_HEH;