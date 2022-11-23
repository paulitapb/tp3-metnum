#include <iostream>

//#include "types.h"
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

vector<string> tests_implementaciones = {
    "test_supernodo_10.txt",
    "test_supernodo_20.txt",
    "test_supernodo_50.txt",

    "test_antisupernodo_10.txt",
    "test_antisupernodo_20.txt",
    "test_antisupernodo_50.txt"};

vector<double> p = {0.9, 0.8, 0.76, 0.85, 0.5, 0.64, 0.3};

vector<string> res_output = {"NO PASO", "OK"};

void armarMatrizA(string test_path, SparseMatrix &W, SparseMatrix &A, double p)
{
    int n = W.outerSize();

    // armamos  D y una identidad
    SparseMatrix I(n, n);
    SparseMatrix D(n, n);

    for (int i = 0; i < n; i++)
    {
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
    A = I - WD;
}

Vector page_rank_EG(string test_path, double p)
{

    SparseMatrix W = read_test(test_path);
    int n = W.outerSize();

    SparseMatrix A(n, n);

    armarMatrizA(test_path, W, A, p);

    Vector e = Vector::Ones(n);

    // Triangulamos el sistema
    elim_gauss(A, e, 1e-5);

    // Resuelvo el sistema
    Vector ranks = backward_sust(A, e);

    // Normalizo la solucion
    normalizar_vector(ranks);

    return ranks;
}

Vector page_rank_Jacobi(string test_path, double p)
{

    SparseMatrix W = read_test(test_path);
    int n = W.outerSize();

    SparseMatrix A(n, n);

    armarMatrizA(test_path, W, A, p);

    Vector xo = Vector::Ones(n);

    Vector v = jacobi(xo, xo, A, 10000, 1e-5);

    // Triangulamos el sistema
    Vector ranks;

    // Normalizo la solucion
    normalizar_vector(ranks);

    return ranks;
}

Vector page_rank_GS(string test_path, double p)
{

    SparseMatrix W = read_test(test_path);
    int n = W.outerSize();

    SparseMatrix A(n, n);

    armarMatrizA(test_path, W, A, p);

    Vector xo = Vector::Ones(n);

    Vector v = gauss_seidel(xo, xo, A, 10000, 1e-5);

    // Triangulamos el sistema
    Vector ranks;

    // Normalizo la solucion
    normalizar_vector(ranks);

    return ranks;
}

void correr_test_catedra()
{

    string output_path = "tests/med_tiempo.txt";
    // archivo de salida que me da el resultado de los test
    ofstream archivo_salida(output_path);
    archivo_salida << "Instancia \t \t \t \t  \t \t \t& n \t \t & m \t \t & Tiempo tardado  \t & test_result \t \t & Error Abs \t \t" << endl;

    int i = 0;
    for (string test : tests)
    {
        i++;
        /* if(i <= 2){
            continue;
        } */
        cout << "corriendo test " << test << " con p " << p[i - 1] << endl;
        // Medir tiempos
        auto inicio = chrono::high_resolution_clock::now();

        // Calculo de rankings
        Vector puntajes_finales = page_rank_EG("tests/" + test, p[i - 1]);
        // Vector puntajes_finales = page_rank_Jacobi(test, p[i-1]);

        auto final = chrono::high_resolution_clock::now();
        chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;

        string res_path = test + ".out";

        pair<bool, double> res_test = resultados_tests(res_path, puntajes_finales);

        string archivo = "tests/" + test;
        ifstream entrada(archivo);
        int n, m;
        entrada >> n >> m;
        if (test == "test_aleatorio_desordenado.txt")
        {
            archivo_salida << "" << test << " \t \t & " << n << " \t \t & " << m << " \t \t & ";
        }
        else
        {
            archivo_salida << "" << test << " \t \t \t \t \t & " << n << " \t \t & " << m << " \t \t & ";
        }
        archivo_salida << fixed << setprecision(8) << tiempo_ejecucion1.count() / 1000 << " seg \t   & \t"
                       << res_output[res_test.first] << "& \t" << res_test.second << endl;
    }
    archivo_salida.close();
}

void correr_test_catedra_experimentacion(int reps)
{

    int i = 0;
    for (string test : tests)
    {
        i++;
        if (i <= 2)
        {
            continue;
        }

        vector<double> tiempos(reps, 0);

        cout << "corriendo test " << test << " con p " << p[i - 1] << endl;
        for (int j = 0; j < reps; j++)
        {

            // Medir tiempos
            auto inicio = chrono::high_resolution_clock::now();

            // Calculo de rankings
            //Vector puntajes_finales = page_rank_EG("tests/" + test, p[i - 1]);
            Vector puntajes_finales = page_rank_GS("tests/" + test, p[i-1]);

            auto final = chrono::high_resolution_clock::now();
            chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;

            tiempos[j] = tiempo_ejecucion1.count();
        }
        std::string res_out = (string) "tiempos_exp/" + test + "GS" + ".exp";

        write_time_results(res_out, tiempos);
    }
}

void correr_test_nuestros(int reps)
{

    string output_path = "test_nuestros/med_tiempo.txt";
    // archivo de salida que me da el resultado de los test
    ofstream archivo_salida(output_path);
    archivo_salida << "Instancia \t \t \t \t  \t \t \t& n \t \t & m \t \t & Tiempo tardado  \t & test_result \t \t & Error Abs \t \t" << endl;
    double p = 0.75;
    int i = 0;

    for (string test : tests_implementaciones)
    {
        for (int i = 0; i < reps; i++)
        {

            cout << "corriendo test " << test << " con p " << p << endl;

            // Medir tiempos
            auto inicio = chrono::high_resolution_clock::now();

            // Calculo de rankings
            Vector puntajes_finales = page_rank_EG("test_nuestros/" + test, p);

            // Vector puntajes_finales = page_rank_Jacobi("test_nuestros/" + test, p);
            auto final = chrono::high_resolution_clock::now();

            chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;

            string res_path = test + ".out";

            // pair<bool, double> res_test = resultados_tests(res_path, puntajes_finales);

            string archivo = "test_nuestros/" + test;
            ifstream entrada(archivo);
            int n, m;
            entrada >> n >> m;

            if (test == "test_aleatorio_desordenado.txt")
            {
                archivo_salida << "" << test << " \t \t & " << n << " \t \t & " << m << " \t \t & ";
            }
            else
            {
                archivo_salida << "" << test << " \t \t \t \t \t & " << n << " \t \t & " << m << " \t \t & ";
            }

            archivo_salida << fixed << setprecision(8) << tiempo_ejecucion1.count() / 1000 << " seg \t   & \t" << endl;
            //<< res_output[res_test.first] << "& \t" << res_test.second <<endl;
        }
    }
    archivo_salida.close();
}

void correr_test_nuestros_experimentacion(int reps)
{
    double p = 0.75;
    int i = 0;
    for (string test : tests_implementaciones)
    {
        cout << "corriendo test " << test << " con p " << p << endl;
        vector<double> tiempos(reps, 0);

        for (int j = 0; j < reps; j++)
        {
            // Medir tiempos
            auto inicio = chrono::high_resolution_clock::now();

            // Calculo de rankings
            Vector puntajes_finales = page_rank_EG("test_nuestros/" + test, p);
            //Vector puntajes_finales = page_rank_GS("test_nuestros/" + test, p);

            auto final = chrono::high_resolution_clock::now();
            chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;

            tiempos[j] = tiempo_ejecucion1.count();
        }
        std::string res_out = (string) "tiempos_exp/" + test + "EG" + ".exp";

        write_time_results(res_out, tiempos);
    }
}

int main()
{
    //correr_test_catedra();
    correr_test_catedra_experimentacion(1000);
    //correr_test_nuestros_experimentacion(1000);

    return 0;
}
