#include <iostream>

#include "types.h"
#include "utils.h"
#include "eigen.h"

using namespace std;

int main()
{
    // Centralidad del autovector del club de karate
    SparseMatrix mat_karate = read_matrix_karate();

    pair<double, Vector> centralidad_karate = power_iteration(mat_karate, 1000, pow(10, -10));

    Eigen::MatrixXd eigvect(1, 34);
    eigvect.row(0) = centralidad_karate.second;

    out_eigvectors(eigvect, "../_outs/centralidad_karate.out");

    // facebook
    SparseMatrix ego_facebook = read_ego_facebook();

    return 0;
}
