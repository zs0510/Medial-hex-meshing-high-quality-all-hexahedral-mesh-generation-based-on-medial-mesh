//
// Created by teitoku on 2021-11-24.
//



#include "TetDeformation.h"
TetDeformation::TetDeformation(MeshKernel::TetMesh * mesh){
    this->mesh=mesh;
}

void TetDeformation::Execute(std::vector<std::pair<MeshKernel::iGameVertexHandle,MeshKernel::iGameVertex> > &point_fixed,
                          std::vector<std::pair<MeshKernel::iGameVertexHandle,MeshKernel::iGameVertex> > &point_moved ){
    std::cout <<"start tet deformation fix points number :" <<point_fixed.size() << "  move points number  " << point_moved.size()<< std::endl;
    std::set<int> handles_f;
    for(auto i:point_fixed){
        handles_f.insert(i.first.idx());
    }
    std::vector<int> handles_m;
    std::vector<MeshKernel::iGameVertex> handles_m_pos;
    std::unordered_map<MeshKernel::iGameVertexHandle,MeshKernel::iGameVertex> pos_mesh_ref;
    for(auto i : mesh->allvertices()){
        pos_mesh_ref[i.first]=i.second;
    }

    for(int i=0;i<point_moved.size();i++){
        handles_m.push_back(point_moved[i].first.idx());
        handles_m_pos.push_back(point_moved[i].second);
        mesh->vertices(point_moved[i].first)=point_moved[i].second;
    }

    std::set<int> handles = handles_f;
    handles.insert(handles_m.begin(), handles_m.end());
    int nv = mesh->VertexSize();
    std::vector<Eigen::Matrix3d> Lts;
    std::vector<Eigen::Matrix3d> A;
    std::vector<Eigen::Matrix3d> A_inverse;
    A.resize( mesh->CellSize());
    A_inverse.resize( mesh->CellSize());

    for(int i=0;i< mesh->CellSize();i++) {
        Vex p0 =  pos_mesh_ref[mesh->cells((CH(i))).vh(0)];
        Vex p1 =  pos_mesh_ref[mesh->cells((CH(i))).vh(1)];
        Vex p2 =  pos_mesh_ref[mesh->cells((CH(i))).vh(2)];
        Vex p3 =  pos_mesh_ref[mesh->cells((CH(i))).vh(3)];

        Vex v1 = p0 - p1;
        Vex v2 = p0 - p2;
        Vex v3 = p0 - p3;

        Eigen::Vector3d e1(v1.x(), v1.y(), v1.z());
        Eigen::Vector3d e2(v2.x(), v2.y(), v2.z());
        Eigen::Vector3d e3(v3.x(), v3.y(), v3.z());
        A[i].col(0)=e1;
        A[i].col(1)=e2;
        A[i].col(2)=e3;
        A_inverse[i]=A[i].inverse();
    }

    Lts.resize(mesh->CellSize());
    for (int iter = 0; iter < 5; iter++) {
        //local
        for(int i=0;i< mesh->CellSize();i++){
            Vex v1p = mesh->vertices(mesh->cells((CH(i))).vh(0))-
                    mesh->vertices(mesh->cells((CH(i))).vh(1));
            Vex v2p = mesh->vertices(mesh->cells((CH(i))).vh(0))-
                    mesh->vertices(mesh->cells((CH(i))).vh(2));
            Vex v3p = mesh->vertices(mesh->cells((CH(i))).vh(0))-
                    mesh->vertices(mesh->cells((CH(i))).vh(3));
            Eigen::Vector3d e1p(v1p.x(),v1p.y(),v1p.z());
            Eigen::Vector3d e2p(v2p.x(),v2p.y(),v2p.z());
            Eigen::Vector3d e3p(v3p.x(),v3p.y(),v3p.z());

            Eigen::Matrix3d B;

            B.col(0)=e1p;
            B.col(1)=e2p;
            B.col(2)=e3p;

            Eigen::Matrix3d J = B*A_inverse[i] ;

            Eigen::JacobiSVD<Eigen::Matrix3d> svd(J, Eigen::ComputeFullU| Eigen::ComputeFullV);
            Eigen::Matrix3d U = svd.matrixU();
            Eigen::Matrix3d V = svd.matrixV();

            Eigen::Matrix3d R = U * V.transpose();

            if (R.determinant() < 0)
            {
                U(0, 2) *= -1;
                U(1, 2) *= -1;
                U(2, 2) *= -1;
                R = U * V.transpose();
            }
            Lts[i] = R;
        }
        //global
        std::vector<Eigen::Triplet<double> > trivec;
        Eigen::MatrixX3d b;
        b.resize(nv, 3);
        b.setZero();

        for(auto i:point_fixed){
            trivec.emplace_back(i.first.idx(),i.first.idx(),1);
            b(i.first.idx(),0)=i.second.x();
            b(i.first.idx(),1)=i.second.y();
            b(i.first.idx(),2)=i.second.z();
        }

        for(int i=0;i<point_moved.size();i++){
            trivec.emplace_back(point_moved[i].first.idx(),point_moved[i].first.idx(),1);
            b(point_moved[i].first.idx(),0)=point_moved[i].second.x();
            b(point_moved[i].first.idx(),1)=point_moved[i].second.y();
            b(point_moved[i].first.idx(),2)=point_moved[i].second.z();
        }

        for(int i=0;i< mesh->CellSize();i++){
            for(int j=0;j<3;j++) {
                double kj = A_inverse[i].coeffRef(0, j);
                double kk = A_inverse[i].coeffRef(1, j);
                double kl = A_inverse[i].coeffRef(2, j);
                double i_factor = kj + kk+ kl;
                double j_factor = -kj;
                double k_factor = -kk;
                double l_factor = - kl;

                double Rtx = Lts[i].coeffRef(0,j);
                double Rty = Lts[i].coeffRef(1,j);
                double Rtz = Lts[i].coeffRef(2,j);
                double prime_i = 2 * i_factor;
                double prime_j = 2 * j_factor;
                double prime_k = 2 * k_factor;
                double prime_l = 2 * l_factor;
                if(!handles.count(mesh->cells((CH(i))).vh(0).idx()) ){
                    trivec.emplace_back(mesh->cells((CH(i))).vh(0), mesh->cells((CH(i))).vh(0), prime_i * i_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(0), mesh->cells((CH(i))).vh(1), prime_i * j_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(0), mesh->cells((CH(i))).vh(2), prime_i * k_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(0), mesh->cells((CH(i))).vh(3), prime_i * l_factor);
                    b(mesh->cells((CH(i))).vh(0).idx(),0) += prime_i * Rtx;
                    b(mesh->cells((CH(i))).vh(0).idx(),1) += prime_i * Rty;
                    b(mesh->cells((CH(i))).vh(0).idx(),2) += prime_i * Rtz;
                }
                if(!handles.count(mesh->cells((CH(i))).vh(1).idx())) {
                    trivec.emplace_back(mesh->cells((CH(i))).vh(1).idx(), mesh->cells((CH(i))).vh(0).idx(), prime_j * i_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(1).idx(), mesh->cells((CH(i))).vh(1).idx(), prime_j * j_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(1).idx(), mesh->cells((CH(i))).vh(2).idx(), prime_j * k_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(1).idx(), mesh->cells((CH(i))).vh(3).idx(), prime_j * l_factor);
                    b(mesh->cells((CH(i))).vh(1).idx(),0) += prime_j * Rtx;
                    b(mesh->cells((CH(i))).vh(1).idx(),1) += prime_j * Rty;
                    b(mesh->cells((CH(i))).vh(1).idx(),2) += prime_j * Rtz;

                }
                if(!handles.count(mesh->cells((CH(i))).vh(2).idx())) {
                    trivec.emplace_back(mesh->cells((CH(i))).vh(2).idx(),mesh->cells((CH(i))).vh(0).idx(), prime_k * i_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(2).idx(),mesh->cells((CH(i))).vh(1).idx(), prime_k * j_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(2).idx(), mesh->cells((CH(i))).vh(2).idx(), prime_k * k_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(2).idx(), mesh->cells((CH(i))).vh(3).idx(), prime_k * l_factor);
                    b(mesh->cells((CH(i))).vh(2).idx(),0) += prime_k * Rtx;
                    b(mesh->cells((CH(i))).vh(2).idx(),1) += prime_k * Rty;
                    b(mesh->cells((CH(i))).vh(2).idx(),2) += prime_k * Rtz;
                }
                if(!handles.count(mesh->cells((CH(i))).vh(3).idx())) {
                    trivec.emplace_back(mesh->cells((CH(i))).vh(3).idx(), mesh->cells((CH(i))).vh(0).idx(), prime_l * i_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(3).idx(), mesh->cells((CH(i))).vh(1).idx(), prime_l * j_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(3).idx(), mesh->cells((CH(i))).vh(2).idx(), prime_l * k_factor);
                    trivec.emplace_back(mesh->cells((CH(i))).vh(3).idx(), mesh->cells((CH(i))).vh(3).idx(), prime_l * l_factor);
                    b(mesh->cells((CH(i))).vh(3).idx(),0) += prime_l * Rtx;
                    b(mesh->cells((CH(i))).vh(3).idx(),1) += prime_l * Rty;
                    b(mesh->cells((CH(i))).vh(3).idx(),2) += prime_l * Rtz;
                }
            }
        }
        Eigen::SparseMatrix<double> smat;
        smat.resize(nv, nv);
        smat.setFromTriplets(trivec.begin(),trivec.end());
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(smat);
        solver.factorize(smat);
        Eigen::MatrixX3d uv;
        uv.resize(nv, 3);
        uv.col(0) = solver.solve(b.col(0));
        uv.col(1) = solver.solve(b.col(1));
        uv.col(2) = solver.solve(b.col(2));

        for(auto j: mesh->allvertices()){
            mesh->vertices(j.first).x() = uv(j.first.idx(), 0);
            mesh->vertices(j.first).y() = uv(j.first.idx(), 1);
            mesh->vertices(j.first).z() = uv(j.first.idx(), 2);
        }
    }
}