#include <set>
#include <Kernel/Mesh.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <omp.h>

class Deformation{
private:
    MeshKernel::SurfaceMesh * mesh;
    std::vector<double> weights;
    bool use_cotweight;
public:
    Deformation(){}
    Deformation(MeshKernel::SurfaceMesh * mesh, bool _use_cotweight = true);
    void Execute(std::vector<int> &point_fixed,
                 std::vector<std::pair<MeshKernel::iGameVertexHandle, MeshKernel::iGameVertex> > &point_moved );
    void initCotWeight();
};