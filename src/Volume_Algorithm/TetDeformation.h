//
// Created by teitoku on 2021-11-24.
//

#ifndef IGAMEVIEW_TETDEFORMATION_H
#define IGAMEVIEW_TETDEFORMATION_H

#endif //IGAMEVIEW_TETDEFORMATION_H
#include <set>
#include <Kernel/Mesh.h>

#ifndef IGMAE_KERNEL_SIMPLIFY
#define IGMAE_KERNEL_SIMPLIFY
typedef MeshKernel::iGameVertex Vex;
typedef MeshKernel::iGameVertex Vec;
typedef MeshKernel::iGameVertexHandle VH;
typedef MeshKernel::iGameEdgeHandle EH;
typedef MeshKernel::iGameFaceHandle FH;
typedef MeshKernel::iGameCellHandle CH;
#endif

class TetDeformation{
private:
    MeshKernel::TetMesh * mesh;
public:
    TetDeformation(){}
    TetDeformation(MeshKernel::TetMesh * mesh);
    void Execute(std::vector<std::pair<MeshKernel::iGameVertexHandle,MeshKernel::iGameVertex> > &point_fixed,
                 std::vector<std::pair<MeshKernel::iGameVertexHandle,MeshKernel::iGameVertex> > &point_moved );
};