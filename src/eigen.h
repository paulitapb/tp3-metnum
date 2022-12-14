#pragma once
#include "types.h"

void elim_gauss(SparseMatrix &A, Vector &v, double epsilon);
Vector backward_sust(SparseMatrix &A, Vector &b);
tuple<Vector, int> jacobi(Vector &x, Vector &b, SparseMatrix &A, int k, double epsilon);
tuple<Vector, int> gauss_seidel(Vector &x, Vector &b, SparseMatrix &A, int k, double epsilon);
