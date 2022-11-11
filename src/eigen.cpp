#include <algorithm>
#include <chrono>
#include <iostream>
#include "eigen.h"

using namespace std;


// Matriz, dinámica, de tamaño arbitrario (X), con elementos del tipo Double(D)
using Eigen::MatrixXd;

// Matriz de numero de filas arbitrario (X) , de 1 columna y con elementos del 
//tipo Double(D)
using Eigen::VectorXd;

void elim_gauss(SparseMatrix &A, double epsilon){
    for(int i = 0; i< A.outerSize(); i++){ //Por cada fila pivot 
		cout << i <<endl; 
		double aii = A.coeff(i, i);
		if(abs(aii) > epsilon){
			
			for(int j = i+1; j < A.innerSize(); j++){ 
				SparseMatrix::InnerIterator it(A, j); 
				double mji = A.coeff(j, i)/aii;

				//It sobre la fila pivot 
				for(SparseMatrix::InnerIterator itPivot(A, i);itPivot; ++itPivot){ 

					if(itPivot.col() == it.col()){
						A.coeffRef(j, it.col()) = it.value() - mji *itPivot.value(); 
						++it; 
					}//
					else if(itPivot.col() < it.col()){
						A.coeffRef(j, it.col()) = - mji *itPivot.value();  
					}
				}
			}	
		}else{
			if(A.row(i).nonZeros() == A.outerSize()){
				cout << "No puedo pivotear" <<endl; 
			}
		}	
	}
}


/* vector<double> backward_sust2(MatrizRalaCSR &A){
    vector<double> res(A.n(), 0); 
    
    for(int i = A.n()-1; i >= 0; i--){ 
        double suma = 0; 
        vector<pair<int, double>> filai = A.row(i); 
         
        if(filai.size() < 2){
            cout << "Hay una variable libre" <<endl; 
            break;
        }

        int k = 0; 
        for(int j = filai.size()-1; j >= 0; j--){
            
            if(i == filai[j].first) continue; 
            
            suma += filai[j].second * res[filai[j].first];
             
        } 
        double aii = A.dameValor(i, i);
        
        res[i]= (A.dameValor(i, A.m()-1) - suma) / aii; 

    } 
    return res; 
} */

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