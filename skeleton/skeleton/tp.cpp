#include <Eigen/Sparse>
#include <vector>
#include <iostream>
 
typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n);
void saveAsBitmap(const Eigen::VectorXd& x, int n, const char* filename);
typedef vector< tuple<int, int> > ejes;

 
int main(int argc, char** argv){

std::vector<T> tripletList;
tripletList.reserve(12);
ejes ejemplo;
tl.push_back( tuple<int, int, int>(1,4,1) );
tl.push_back( tuple<int, int, int>(1,3,1) );
tl.push_back( tuple<int, int, int>(1,2,1) );

SparseMatrixType mat(4,4);
mat.setFromTriplets(tl.begin(), tl.end());
cout << mat << endl;
}
