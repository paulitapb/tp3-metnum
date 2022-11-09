#include "utils.h"

SparseMatrix read_matrix_karate(){
    double epsilon = pow(10, -5);

    ifstream input("../../datasets/karateclub_matriz.txt");  
     
    std::vector<T> tripletList;
    tripletList.reserve(156);

    for(int i = 0; i < 34; i++){
        for(int j = 0; j < 34; j++){
            double e;
            input >> e; 
            if(abs(e-1) < epsilon ){
                tripletList.push_back(T(i,j,e));                  
            }
        }
    }

    SparseMatrix mat_karate(34, 34);
    mat_karate.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat_karate; 
}

void print_sparce_matrix(SparseMatrix &m){
    cout<< MatrixXd(m) <<endl;
}


SparseMatrix read_ego_facebook(){
    double epsilon = pow(10, -5);

    ifstream input("../../datasets/ego-facebook.edges");  
    int m = 28048; 
    int n = 3436; 

    std::vector<T> tripletList;
    tripletList.reserve(m);

    for(int i = 0; i < m; i++){
        double a, b;
        input >> a >> b; 
        tripletList.push_back(T(a-1,b-1,1));
    }

    SparseMatrix mat_ego_face(n, n);
    mat_ego_face.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat_ego_face; 
}

void out_eigvalues(Vector &eigvals, string path){
    ofstream output(path); 
    for(Vector::iterator it = eigvals.begin(); it != eigvals.end();  it++){
        output << *it << "\n"; 
    }
}

void out_eigvectors(Matrix eigvect, string path){
    ofstream output(path); 
    for (int k = 0; k< eigvect.outerSize(); ++k){ 
        for(Matrix::InnerIterator it(eigvect, k); it; ++it){
            output << it.value()<< "\t"; 
        }
        output <<endl; 
    }
}


