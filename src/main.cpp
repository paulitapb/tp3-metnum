#include <iostream>

#include "types.h"
#include "utils.h"
#include "eigen.h"
#include <chrono>

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

vector<double> page_rank_EG(string test_path, double p){
    SparseMatrix W = read_test(test_path);
    int n = W.outerSize();

    // armamos  D y una identidad
    SparseMatrix I(n, n);
    SparseMatrix D(n, n);

    for (int i = 0; i < n; i++){
        I.coeffRef(i, i) = 1;
        double colSum = W.col(i).sum();
        if (abs(colSum) > epsilon)
        {
            D.coeffRef(i, i) = 1 / colSum;
        }
    }

    // Armamos WD
    SparseMatrix WD(n, n);
    WD = W * D;

    // pW
    WD = p * WD;

    // A = I-pWD
    SparseMatrix A(n, n);

    A = I - WD;

    //print_sparce_matrix(A);
    
    Vector e = Vector::Ones(n);

    //Triangulamos el sistema
    elim_gauss(A, e, 1e-5);

    print_sparce_matrix(A);
    cout << e <<endl; 

    //Resuelvo el sistema
    vector<double> ranks = backward_sust(A, e);

    //Normalizo la solucion
    normalizar_vector(ranks);

    print_vector(ranks); 

    return ranks; 
}

pair<bool, double> resultados_tests(string res_path, vector<double> const &v){
    string archivo = "tests/" + res_path;
    string archivo_out = "test_out/" + res_path;  
    ifstream entrada(archivo);
    ofstream salida(archivo_out); 
    double p;  
    entrada >> p;
    salida << p << endl; 
    vector<double> res(v.size(), 0);
    bool paso_test = true; 
    double suma_error = 0; 
    for(int i = 0; i< v.size();i++){
        entrada >> res[i];  
        salida << fixed << setprecision(4) << v[i] <<endl; 

        suma_error += fabs(v[i]-res[i]); 
        if(fabs(v[i]-res[i]) > 0.0001){
            paso_test &= false; 
        }
    }
    return make_pair(paso_test, suma_error/v.size()); 
}


void correr_test_catedra(){
 
    string output_path = "tests/med_tiempo.txt"; 
    //archivo de salida que me da el resultado de los test 
    ofstream archivo_salida(output_path);
    archivo_salida << "Instancia \t \t \t \t  \t \t \t& n \t \t & m \t \t & Tiempo tardado  \t & test_result \t \t & Error Abs \t \t" << endl;
    
    int i = 0;
    for (string test : tests) {
        i++;
        if(i <= 2){
            continue;
        }
        cout << "corriendo test "  << test << " con p " << p[i-1] <<endl;
        
        //Medir tiempos 
        auto inicio = chrono::high_resolution_clock::now();

        //Calculo de rankings   
        vector<double> puntajes_finales = page_rank_EG(test, p[i-1]); 
        
        auto final = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;
        
        string res_path = test + ".out";  
         
        pair<bool, double> res_test = resultados_tests(res_path, puntajes_finales); 

        string archivo = "tests/" + test;
        ifstream entrada(archivo);
        int n,m;
        entrada >> n >> m;
        if(test == "test_aleatorio_desordenado.txt"){
            archivo_salida << "" << test << " \t \t & " << n << " \t \t & "  << m << " \t \t & ";
        }else{
            archivo_salida << "" << test << " \t \t \t \t \t & " << n << " \t \t & "  << m << " \t \t & ";
        }
        archivo_salida << fixed << setprecision(8) << tiempo_ejecucion1.count() / 1000 << " seg \t   & \t" 
        << res_output[res_test.first] << "& \t" << res_test.second <<endl;
    }
    archivo_salida.close();
}


int main(){
    correr_test_catedra(); 

    return 0;
}
