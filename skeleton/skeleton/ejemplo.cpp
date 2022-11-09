#include <Eigen/Sparse>
#include <vector>
#include <iostream>


/**
https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
**/


typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n);
void saveAsBitmap(const Eigen::VectorXd& x, int n, const char* filename);

 
int main(int argc, char** argv){

std::vector<T> tl;
tl.reserve(12);
tl.push_back( T(0,3,1) );
tl.push_back( T(0,1,1) );
tl.push_back( T(0,2,1) );

SpMat mat(4,4);
mat.setFromTriplets(tl.begin(), tl.end());

std::cout << mat << std::endl;
std::cout << -mat << std::endl;

std::cout << -mat * 2 << std::endl;

std::cout << mat.transpose() << std::endl;

// Acá necesitamos casting para que la suma este definida
std::cout << SpMat(mat.transpose()) + mat << std::endl;

}
