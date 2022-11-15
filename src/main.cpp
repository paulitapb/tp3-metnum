#include <iostream>

#include "types.h"
#include "utils.h"
#include "eigen.h"

using namespace std;
double epsilon = 1e-5;
vector<string> tests = {"test_15_segundos.txt",  
                        "test_30_segundos.txt",
                        "test_aleatorio_desordenado.txt",
                        "test_aleatorio.txt",
                        "test_completo.txt",
                        "test_sin_links.txt",
                        "test_trivial.txt"};

vector<double> p = {0.9, 0.8, 0.76, 0.85, 0.5, 0.64, 0.3}; 

vector<string> res_output = {"NO PASO", "OK"};


vector<double> page_rank(string test_path, double p){  
    
    SparseMatrix W = read_test("test_aleatorio.txt");

    //print_sparce_matrix(W);
    int n = W.outerSize();

    //armamos  D y una identidad 
    SparseMatrix I(n,n);
    SparseMatrix D(n,n);

    for(int i = 0; i < n; i++){
        I.coeffRef(i, i) = 1;
        //I.coeffRef(i, n) = 1;
        double colSum = W.col(i).sum();
        if(abs(colSum) > epsilon){
            D.coeffRef(i, i) = 1/colSum;
        }
    } 

    cout << "Identidad" << endl;
    print_sparce_matrix(I);

    cout << "Diagonal" << endl;
    print_sparce_matrix(D);
    

    //Armamos WD
    SparseMatrix WD(n,n); 
    WD = W*D;

    cout << "WD" << endl;
    print_sparce_matrix(WD);

    // pW
    WD = p* WD;

    cout << "pWD" << endl;
    print_sparce_matrix(WD);
     
    
    // A = I-pWD
    SparseMatrix A(n, n);
    A = I - WD; 
    cout << "A" << endl;
    print_sparce_matrix(A);
    
    // Generamos b
    vector <double> b(n, 1);

    //Aplicamos EG
    elim_gauss(A, b, epsilon);

    cout << "Con EG A vale (wodka):" << endl;


    return {};
    /*
    vector<double> e(w.n(), 1); 

    //Agregamos e para obtener el sistema Ax = e  
    A.agregarColumna(e);

    //Triangulamos el sistema
    elim_gauss2(A);

    //Resuelvo el sistema
    vector<double> ranks = backward_sust2(A);

    //Calculo error
      
    //Vector z
    vector<double> z(w.n());
    for (int i = 0; i < z.size(); i++){
        if(cjs[i] == 0){
            z[i] = 1/w.n();
        }else{
            z[i] = (1-p)/w.n();
        }
    }
    
    //Calculo e * z traspuesta
    vector<double> ez_a;
    vector<int> ez_ja;
    vector<int> ez_ia(w.n()+1, 0);
    MatrizRalaCSR ZE = MatrizRalaCSR(w.n(), w.n(), ez_a, ez_ja, ez_ia);
    
    for(int i = 0; i < z.size(); i++){
        for(int j = 0; j < z.size(); j++){ 
            ZE.asignarValor(i,j, -z[i]); //Multiplico x -1 para despues restar 
        }
    }
    MatrizRalaCSR x = MatrizRalaCSR(w.n()+1, w.n(), ez_a, ez_ja, ez_ia);
    
    
    x.agregarColumna(ranks); 
    // Generamos A2 y calculamos |Ax'- x'| ~= 0
    MatrizRalaCSR A2 = restar_ralas2(WD, ZE); 
    A2 = multiplicar_ralas2(A2,x); 
    A2 = restar_ralas2(A2,x);   

    bool aproxima_cero = true;
    for (int i = 0; i < A2.n(); i++)
    {
        if(fabs(A2.m_real(i)-0) > 0.00001) aproxima_cero = false; 
    }
    cout << "x es solucion aproximada: " <<  res_output[aproxima_cero] <<endl; 

    //Normalizo la solucion
    normalizar_vector(ranks); 
    
    return ranks; */
}

int main()
{
    page_rank("test_aleatorio.txt", 0.9);
   //SparseMatrix W = read_test("test_aleatorio.txt");
   //print_sparce_matrix(W);

   return 0;
}
