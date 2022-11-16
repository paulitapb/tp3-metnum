#pragma once


#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMatrix;

typedef Eigen::VectorXd Vector;
typedef Eigen::SparseVector<double> SVector;


typedef Eigen::Triplet<double> T; //to build sparce matrix efficiently

