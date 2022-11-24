#include <algorithm>
#include <chrono>
#include <iostream>

#include "eigen.h"

#define CHECKPIVOT 0;
using namespace std;

// Matriz, dinámica, de tamaño arbitrario (X), con elementos del tipo Double(D)
using Eigen::MatrixXd;

// Matriz de numero de filas arbitrario (X) , de 1 columna y con elementos del
// tipo Double(D)
using Eigen::VectorXd;

void elim_gauss(SparseMatrix &A, Vector &v, double epsilon)
{
	// CAMBIAR COEFFREFF POR ValueRef

	for (int i = 0; i < A.outerSize() - 1; i++)
	{ // Por cada fila pivot
		double aii = A.coeff(i, i);

		if (abs(aii) > epsilon)
		{
			for (int j = i + 1; j < A.outerSize(); j++)
			{ // Recorro cada una de las filas de adentro

				SVector filaPivot = A.innerVector(i);
				SVector filaJ = A.innerVector(j);

				/* if (abs(A.coeff(j, i)) < epsilon)
				{
					cout <<"Salta fila"<<endl;
					continue;
				} */

				double mji = A.coeff(j, i) / aii;
				// cout << "Mji es " << mji << endl;
				v.coeffRef(j) = v.coeffRef(j) - mji * v.coeffRef(i);
				// It sobre la fila pivot
				SVector::InnerIterator itPivot(filaPivot);
				SVector::InnerIterator it(filaJ);

				while (itPivot)
				{
					if (itPivot.index() == it.index())
					{
						A.coeffRef(j, it.index()) = it.value() - mji * itPivot.value();
						++it;
						++itPivot;
					}
					else if (itPivot.index() > it.index())
					{
						++it;
					}
					else if (itPivot.index() < it.index() or it)
					{
						A.coeffRef(j, itPivot.index()) = -mji * itPivot.value();
						++itPivot;
					}
				}
			}
		}
		else
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

Vector backward_sust(SparseMatrix &A, Vector &b)
{
	// Falta que tome un vector y que tomo el ultimo valor del vector

	Vector res = Vector::Zero(A.outerSize());

	for (int i = A.outerSize() - 1; i >= 0; i--)
	{
		double suma = 0;
		SVector filai = A.innerVector(i);

		int k = 0;
		for (SVector::ReverseInnerIterator itFila(filai); itFila; --itFila)
		{
			if (i == itFila.index())
				continue;
			suma += itFila.value() * res(itFila.index());
		}
		double aii = A.coeff(i, i);
		res(i) = (b[i] - suma) / aii;
	}
	return res;
}

double iterativeSum(Vector x, SparseMatrix &A, int i)
{
	// como el x va a tener actualizadas las posiciones hasta i-1 y las siguientes van a ser de la iteracion anterior (GS),
	// la suma es la misma para gauss-seidel y jacobi, nada mas que a uno le paso el vector que estoy actualizando y al otro el previo
	double sum = 0;
	for (int j = 0; j < x.size(); j++)
	{
		if (j == i)
			continue;

		sum += A.coeff(i, j) * x(j);
	}
	return sum;
}

bool areEqual(Vector x1, Vector x2, double e)
{
	if (x1.size() != x2.size())
		return false;
	for (int i = 0; i < x1.size(); i++)
	{
		if (abs(x1(i) - x2(i)) > e)
		{
			return false;
		}
	}
	return true;
}

tuple<Vector, int> jacobi(Vector x, Vector b, SparseMatrix &A, int k, double epsilon)
{
	Vector prev = x + Vector::Ones(x.size());
	int iter = 0;
	while (!areEqual(prev, x, epsilon) && iter < k)
	{
		prev = x;
		for (int i = 0; i < x.size(); i++)
		{
			x(i) = (b(i) - iterativeSum(prev, A, i)) / A.coeff(i, i);
		}
		iter++;
	}
	// cout << iter <<endl;
	return make_tuple(x, iter);
}

tuple<Vector, int> gauss_seidel(Vector x, Vector b, SparseMatrix &A, int k, double epsilon)
{
	// el vector previous no lo necesito mas como vector porque no lo voy a recorrer nunca
	// lo unico que quizas puede ser que pase es que queramos comparar todo el vector previo con el ultimo x
	Vector prev = x + Vector::Ones(x.size());
	int iter = 0;
	while (!areEqual(prev, x, epsilon) && iter < k)
	{
		prev = x;
		for (int i = 0; i < x.size(); i++)
		{
			x(i) = (b(i) - iterativeSum(x, A, i)) / A.coeff(i, i);
		}
		iter++;
	}
	// cout << iter <<endl;
	return make_tuple(x, iter);
}