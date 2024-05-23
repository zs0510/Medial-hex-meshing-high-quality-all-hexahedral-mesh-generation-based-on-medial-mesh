//
// Created by teitoku on 2021-11-11.
//

#include "Deformation.h"

Deformation::Deformation(MeshKernel::SurfaceMesh* mesh, bool _use_cotweight) {
    this->mesh = mesh;
    this->use_cotweight = _use_cotweight;
    if (use_cotweight) initCotWeight();// 初始化时就将权值计算好
}
void Deformation::Execute(std::vector<int>& point_fixed,
    std::vector<std::pair<MeshKernel::iGameVertexHandle, MeshKernel::iGameVertex> >& point_moved) {
    std::cout << "start surface mesh deformation fix points number :" << point_fixed.size() << "  move points number  " << point_moved.size() << std::endl;
    std::set<int> handles_f = {};
    for (auto i : point_fixed) {
        handles_f.insert(i);
    }
    std::vector<int> handles_m = {};
    std::vector<MeshKernel::iGameVertex> handles_m_pos = {  };

    for (int i = 0; i < point_moved.size(); i++) {
        handles_m.push_back(point_moved[i].first.idx());
        handles_m_pos.push_back(point_moved[i].second);
    }


    std::set<int> handles = handles_f;
    handles.insert(handles_m.begin(), handles_m.end());

    int nv = mesh->VertexSize();
    std::unordered_map<MeshKernel::iGameVertexHandle, MeshKernel::iGameVertex> pos_mesh_ref;
    for (auto i : mesh->allvertices()) {
        pos_mesh_ref[i.first] = i.second;
    }

    std::vector<Eigen::Triplet<double>> trivec;
    for (auto i : mesh->allvertices()) {
        if (handles.count(i.first.idx()) > 0)
        {
            trivec.emplace_back(i.first.idx(), i.first.idx(), 1.);
            continue;
        }
        double weight_sum = 0.;
        for (auto eh : mesh->NeighborEh(i.first)) {
            auto& edge = mesh->edges(eh);
            auto j = (i.first == edge.vh1()) ? edge.vh2() : edge.vh1();
            double weight_ = (use_cotweight) ? weights[eh] : 2;
            weight_sum += weight_;
            trivec.emplace_back(i.first.idx(), j.idx(), -weight_);
        }
       /* for (auto j : mesh->NeighborVh(i.first)) {
            double weight_ = 2;
            weight_sum += weight_;
            trivec.emplace_back(i.first.idx(), j.idx(), -weight_);
        }*/
        trivec.emplace_back(i.first.idx(), i.first.idx(), weight_sum);
    }
    Eigen::SparseMatrix<double> smat;
    smat.resize(nv, nv);
    smat.setFromTriplets(trivec.begin(), trivec.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(smat);

    Eigen::MatrixX3d uv;
    uv.resize(nv, 3);

    std::vector<Eigen::Matrix3d> Lts;
    Lts.resize(nv);
    Eigen::MatrixX3d b;
    b.resize(nv, 3);

    for (int iter = 0; iter < 5; iter++) {

#pragma omp parallel for
        for (int i = 0; i < nv; ++i) {
            MeshKernel::iGameVertexHandle vh(i);
            auto& v = mesh->vertices(vh);
            Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
            for (auto eh : mesh->NeighborEh(vh)) {
                auto& edge = mesh->edges(eh);
                auto j = (vh == edge.vh1()) ? edge.vh2() : edge.vh1();
                double weight_ = (use_cotweight) ? weights[eh] : 2;
                auto e_ = pos_mesh_ref[vh] - pos_mesh_ref[j];
                auto ep_ = v - mesh->vertices(j);
                Eigen::Vector3d ep(ep_.x(), ep_.y(), ep_.z());
                Eigen::Vector3d e(e_.x(), e_.y(), e_.z());
                J += weight_ * (e * ep.transpose());
            }
           /* for (auto j : mesh->NeighborVh(vh))
            {
                auto e_ = pos_mesh_ref[vh] - pos_mesh_ref[j];
                auto ep_ = v - mesh->vertices(j);
                double weight_ = 2;
                Eigen::Vector3d ep(ep_.x(), ep_.y(), ep_.z());
                Eigen::Vector3d e(e_.x(), e_.y(), e_.z());
                J += weight_ * (e * ep.transpose());
            }*/
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Matrix3d U = svd.matrixU();
            Eigen::Matrix3d V = svd.matrixV();
            Eigen::Matrix3d R = V * U.transpose();
            if (R.determinant() < 0)
            {
                U(0, 2) *= -1;
                U(1, 2) *= -1;
                U(2, 2) *= -1;
                R = V * U.transpose();
            }
            Lts[vh.idx()] = R;
        }

        /*for (auto i : mesh->allvertices())
        {
            Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
            for (auto j : mesh->NeighborVh(i.first))
            {
                auto e_ = pos_mesh_ref[i.first] - pos_mesh_ref[j];
                auto ep_ = i.second - mesh->vertices(j);
                double weight_ = 2;
                Eigen::Vector3d ep(ep_.x(), ep_.y(), ep_.z());
                Eigen::Vector3d e(e_.x(), e_.y(), e_.z());
                J += weight_ * (e * ep.transpose());
            }
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::Matrix3d U = svd.matrixU();
            Eigen::Matrix3d V = svd.matrixV();
            Eigen::Matrix3d R = V * U.transpose();
            if (R.determinant() < 0)
            {
                U(0, 2) *= -1;
                U(1, 2) *= -1;
                U(2, 2) *= -1;
                R = V * U.transpose();
            }
            Lts[i.first.idx()] = R;
        }*/
#pragma omp parallel for
        for (int i = 0; i < nv; ++i) {
            MeshKernel::iGameVertexHandle vh(i);
            auto& v = mesh->vertices(vh);
            Eigen::Vector3d b_tmp(0., 0., 0.);
            for (auto eh : mesh->NeighborEh(vh)) {
                auto& edge = mesh->edges(eh);
                auto j = (vh == edge.vh1()) ? edge.vh2() : edge.vh1();
                double weight_ = (use_cotweight) ? weights[eh] : 2;
                weight_ *= 0.5;
                auto ep_ = pos_mesh_ref[vh] - pos_mesh_ref[j];
                Eigen::Vector3d ep(ep_.x(), ep_.y(), ep_.z());
                Eigen::Matrix3d JR = Lts[vh.idx()] + Lts[j.idx()];
                b_tmp += weight_ * (JR * ep);
            }
            /*for (auto j : mesh->NeighborVh(vh))
            {
                auto ep_ = pos_mesh_ref[vh] - pos_mesh_ref[j];
                Eigen::Vector3d ep(ep_.x(), ep_.y(), ep_.z());
                Eigen::Matrix3d JR = Lts[vh.idx()] + Lts[j.idx()];
                double weight_ = 1;
                b_tmp += weight_ * (JR * ep);
            }*/
            b(vh.idx(), 0) = b_tmp[0];
            b(vh.idx(), 1) = b_tmp[1];
            b(vh.idx(), 2) = b_tmp[2];
        }

        /*for (auto i : mesh->allvertices())
        {
            Eigen::Vector3d b_tmp(0., 0., 0.);
            for (auto j : mesh->NeighborVh(i.first))
            {
                auto ep_ = pos_mesh_ref[i.first] - pos_mesh_ref[j];
                Eigen::Vector3d ep(ep_.x(), ep_.y(), ep_.z());
                Eigen::Matrix3d JR = Lts[i.first.idx()] + Lts[j.idx()];
                double weight_ = 1;
                b_tmp += weight_ * (JR * ep);
            }
            b(i.first.idx(), 0) = b_tmp[0];
            b(i.first.idx(), 1) = b_tmp[1];
            b(i.first.idx(), 2) = b_tmp[2];
        }*/
        for (int i : handles_f) {
            MeshKernel::iGameVertexHandle vh(i);
            auto& b_tmp = mesh->vertices(vh);
            /*for (auto j : mesh->allvertices())
                if (j.first.idx() == i)
                    b_tmp = pos_mesh_ref[j.first];*/
            b(i, 0) = b_tmp.x();
            b(i, 1) = b_tmp.y();
            b(i, 2) = b_tmp.z();
        }

        for (int i = 0; i < handles_m.size(); i++)
        {
            auto b_tmp = handles_m_pos[i];
            b(handles_m[i], 0) = b_tmp.x();
            b(handles_m[i], 1) = b_tmp.y();
            b(handles_m[i], 2) = b_tmp.z();
        }
        uv.col(0) = solver.solve(b.col(0));
        uv.col(1) = solver.solve(b.col(1));
        uv.col(2) = solver.solve(b.col(2));
        for (auto i : mesh->allvertices()) {
            mesh->vertices(i.first).x() = uv(i.first.idx(), 0);
            mesh->vertices(i.first).y() = uv(i.first.idx(), 1);
            mesh->vertices(i.first).z() = uv(i.first.idx(), 2);
        }
    }
}

void Deformation::initCotWeight() {

    auto ecnt = mesh->alledges().size();
    weights.resize(ecnt, 0.f);

#pragma omp parallel for
    for (int i = 0; i < ecnt; ++i) {
        MeshKernel::iGameEdgeHandle eh(i);
        auto& edge = mesh->edges(eh);
        auto vh = edge.vh1();
        auto uh = edge.vh2();
        auto& v = mesh->vertices(vh);
        auto& u = mesh->vertices(uh);
        float cot_weight = 0.f;
        for (auto face : mesh->NeighborFh(eh)) {
            int idx = 0;
            auto wh = mesh->faces(face).vh(idx++);
            while (wh == vh || wh == uh) {
                assert(idx < 3);
                wh = mesh->faces(face).vh(idx++);
            }
            auto& w = mesh->vertices(wh);
            Eigen::Vector3f v_pos(v.x(), v.y(), v.z());
            Eigen::Vector3f u_pos(u.x(), u.y(), u.z());
            Eigen::Vector3f w_pos(w.x(), w.y(), w.z());
            Eigen::Vector3f vec1 = (v_pos - w_pos).normalized();
            Eigen::Vector3f vec2 = (u_pos - w_pos).normalized();
            float cos_theta = vec1.dot(vec2);
            float cot_theta = cos_theta / std::sqrt(1 - cos_theta * cos_theta);
            cot_weight += std::max(1E-6F, cot_theta);
        }
        weights[eh] = cot_weight * 0.5f;
#ifdef DEBUG
        std::cout << "weights: " << weights[eh] << std::endl;
#endif // DEBUG
    }
    std::cout << "init cotangent weights success." << std::endl;
}