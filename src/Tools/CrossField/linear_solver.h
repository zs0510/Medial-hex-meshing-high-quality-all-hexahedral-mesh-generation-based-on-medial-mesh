//
// Created by su on 2022/3/13.
//

#ifndef CROSSFIELD__LINEAR_SOLVER_H_
#define CROSSFIELD__LINEAR_SOLVER_H_


#include <vector>

#include <Eigen/Sparse>

enum linearSystemEigenSolver {
  EigenCholeskyLLT,
  EigenCholeskyLDLT,
  EigenSparseLU,
  EigenSparseQR,
  EigenCG,
  EigenCGLeastSquare,
  EigenBiCGSTAB
};

class LinearSystem {
 private:
  Eigen::VectorXd x_;
  Eigen::VectorXd b_;
  Eigen::SparseMatrix<double> a_;
  linearSystemEigenSolver solver_type_;


  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_lu;

  //for lu
  bool factorize_;
  bool preprocess_sparsity_pattern_;

 public:
  LinearSystem();

  void set_solver_type(linearSystemEigenSolver solverName);

  bool IsAllocated() const;
  void Allocate(int nbRows);
  void Clear();
  void ZeroMatrix();

  void ZeroRightHandSide();
  void ZeroSolution();
  int SystemSolve();
  double NormInfRightHandSide() const;
  double NormInfSolution() const;

  void AddToMatrix(int row, int col, const double &val);
  void GetFromMatrix(int row, int col, double &val) const;
  void AddToRightHandSide(int row, const double &val);
  bool AddToRightHandSide(const std::vector<double> &b);
  void GetFromRightHandSide(int row, double &val) const;
  void GetFromSolution(int row, double &val) const;
  void AddToSolution(int row, const double &val);

  bool AddSparseCoefficients(
      const std::vector<std::vector<size_t>> &columns,
      const std::vector<std::vector<double>> &values,
      bool first_time = false);

  bool PreprocessSparsityPattern();
  bool Factorize();

  bool Solve(std::vector<double> &x);
};

bool SolveSparseLinearSystem(
    const std::vector<std::vector<size_t>> &columns,
    const std::vector<std::vector<double>> &values,
    const std::vector<double> &rhs,
    std::vector<double> &x);



#endif //CROSSFIELD__LINEAR_SOLVER_H_
