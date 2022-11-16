#include <algorithm>
#include <chrono>
#include <iostream>
#include "eigen.h"

#define CHECKPIVOT 1;
using namespace std;

// Matriz, dinámica, de tamaño arbitrario (X), con elementos del tipo Double(D)
using Eigen::MatrixXd;

// Matriz de numero de filas arbitrario (X) , de 1 columna y con elementos del
// tipo Double(D)
using Eigen::VectorXd;

void elim_gauss(SparseMatrix &A, Vector &v, double epsilon){
	
	for (int i = 0; i < A.outerSize()-1; i++)
	{ // Por cada fila pivot
		double aii = A.coeff(i, i);
		//cout << "aii " << aii <<endl; 
		if (abs(aii) > epsilon)
		{
			for (int j = i + 1; j < A.outerSize(); j++)
			{ //Recorro cada una de las filas de adentro
					
				SVector filaPivot = A.innerVector(i); 
				SVector filaJ = A.innerVector(j); 	

				if (abs(A.coeff(j, i)) < epsilon)
				{
					continue;
				}
				
				double mji = A.coeff(j, i) / aii;
				//cout << "Mji es " << mji << endl;
				v.coeffRef(j) = v.coeffRef(j) - mji*v.coeffRef(i); 
				// It sobre la fila pivot
				SVector::InnerIterator itPivot(filaPivot);
				SVector::InnerIterator it(filaJ);  

				while (itPivot)
				{
					/* if(i >= 1){
						cout << filaJ << endl; 
						cout << filaPivot <<endl; 
						cout << "it " << it.value() << " " << it.index() <<endl; 
						cout << "itPivot " << itPivot.value() << " " << itPivot.index() <<endl;  
						cout << " pivots" <<endl;
					} */
					if (itPivot.index() == it.index())
					{
						A.coeffRef(j, it.index()) = it.value() - mji * itPivot.value();
						++it;
						++itPivot;
					}else if(itPivot.index() > it.index()){
						++it; 
					}
					else if (itPivot.index() < it.index() or it )
					{
						A.coeffRef(j, itPivot.index()) = -mji * itPivot.value();
						++itPivot; 
					}
				}
			}
			
		}else
		{ // TODO //NO CHEQUEA SI LA COLUMNA ES CERO!!!
#ifdef CHECKPIVOT
			bool rompe = false;
			Eigen::Block<Eigen::SparseMatrix<double, 1, int>, -1, 1, false> columnita = A.col(i);
			for (int k = i + 1; k < columnita.innerSize(); k++)
			{
				if (abs(columnita.coeffRef(k, 0)) < epsilon)
				{
					cout << "No se pivotear" << endl;
					rompe = true;
					break;
				}
			}
			if (rompe)
			{
				break;
			}
#endif
		}
	} 
 
}

vector<double> backward_sust(SparseMatrix &A){
	//Falta que tome un vector y que tomo el ultimo valor del vector
	vector<double> res(A.outerSize(), 0);

	for(int i = A.outerSize()-1; i >= 0; i--){
		double suma = 0;
		SVector filai = A.row(i);

		if(filai.size() < 2){
			cout << "Hay una variable libre" <<endl;
			break;
		}

		int k = 0;
		for(int j = filai.size()-1; j >= 0; j--){

			if(i == filai.coeff(j).first) continue;

			suma += filai(j).second * res[filai(j).first];

		}
		double aii = A.coeff(i, i);

		res[i]= (A.coeff(i, A.outerSiz()-1) - suma) / aii;

	}
	return res;
}

double jacobiSum(Vector x, SparseMatrix &A, int i)
{
	double sum = 0;
	for (int j = 0; j < x.size(); j++)
	{
		if (j == i)
			continue;

		sum += A.coeff(i, j) * x(j);
	}
	return sum;
}

Vector jacobi(Vector x, Vector b, SparseMatrix &A, int k, double epsilon)
{
	Vector prev = Vector::Ones(x.size());
	int iter = 0;
	while (abs(prev(0) - x(0)) >= epsilon || iter < k)
	{
		prev = x;
		for (int i = 0; i < x.size(); i++)
		{
			x(i) = (b(i) - jacobiSum(prev, A, i)) / A.coeff(i, i);
		}
		iter++;
	}
	return x;
}

pair<double, Vector> power_iteration(const Matrix &X, unsigned niter, double eps)
{
	Eigen::VectorXd b = Vector::Random(X.cols());
	bool break_crit = false;

	for (int i = 0; i < niter && !break_crit; i++)
	{
		Eigen::VectorXd new_b = X * b;
		new_b = new_b / new_b.norm();
		double cos_angle = new_b.transpose() * b;
		break_crit = (1 - eps) < cos_angle && cos_angle <= 1;
		b = new_b;
	}

	double eigval = b.transpose() * X * b;
	return make_pair(eigval, b);
}

pair<Vector, Matrix> deflation(const Matrix &X, unsigned num, unsigned num_iter, double epsilon)
{
	Matrix A;
	A = X;
	Eigen::VectorXd eigvalues(num);
	Eigen::MatrixXd eigvectors(A.rows(), num);

	double a = 0;
	Eigen::VectorXd v = Vector::Zero(A.rows());
	for (int i = 0; i < num; i++)
	{
		A = A - (a * v * v.transpose());
		tie(a, v) = power_iteration(A, num_iter, epsilon);
		eigvalues[i] = a;
		eigvectors.col(i) = v;
	}
	return make_pair(eigvalues, eigvectors);
}