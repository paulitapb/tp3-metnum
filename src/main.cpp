#include <iostream>

//#include "types.h"
#include "utils.h"

#include "eigen.h"
#include <chrono>

using namespace std;
double epsilon = 1e-10;
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
    "test_supernodo_30.txt",
    "test_supernodo_40.txt",
    "test_supernodo_50.txt",

    "test_antisupernodo_10.txt",
    "test_antisupernodo_20.txt",
    "test_antisupernodo_30.txt",
    "test_antisupernodo_40.txt",
    "test_antisupernodo_50.txt",

    "clique10.txt",
    "clique20.txt",
    "clique30.txt",
    "clique40.txt",
    "clique50.txt",

    "binomial_graph50_1.txt",
    "binomial_graph50_2.txt",
    "binomial_graph50_3.txt",
    "binomial_graph50_4.txt",
    "binomial_graph50_5.txt",
    "binomial_graph50_6.txt",
    "binomial_graph50_8.txt",

    "density50_10.txt",
    "density50_20.txt",
    "density50_30.txt",
    "density50_40.txt",
    "density50_50.txt"};

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
    elim_gauss(A, e, epsilon);

    // Resuelvo el sistema
    Vector ranks = backward_sust(A, e);

    // Normalizo la solucion
    normalizar_vector(ranks);

    return ranks;
}

pair<Vector, int> page_rank_Jacobi(string test_path, double p)
{

    SparseMatrix W = read_test(test_path);
    int n = W.outerSize();

    SparseMatrix A(n, n);

    armarMatrizA(test_path, W, A, p);

    Vector xo = Vector::Ones(n);

    // Triangulamos el sistema
    Vector ranks;
    int iters;
    tie(ranks, iters) = jacobi(xo, xo, A, 10000, epsilon);

    // Normalizo la solucion
    normalizar_vector(ranks);

    return make_pair(ranks, iters);
}

pair<Vector, int> page_rank_GS(string test_path, double p)
{

    SparseMatrix W = read_test(test_path);
    int n = W.outerSize();

    SparseMatrix A(n, n);

    armarMatrizA(test_path, W, A, p);

    Vector xo = Vector::Ones(n);

    // Triangulamos el sistema
    Vector ranks;
    int iters;
    tie(ranks, iters) = gauss_seidel(xo, xo, A, 10000, epsilon);

    // Normalizo la solucion
    normalizar_vector(ranks);

    return make_pair(ranks, iters);
}

void correr_test_catedra()
{

    string output_path = "tests/med_tiempoGS.txt";
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
        Vector puntajes_finales;
        //puntajes_finales = page_rank_EG("tests/" + test, p[i - 1]);
        std::pair<Vector, int> r = page_rank_GS("tests/" + test, p[i - 1]);
        puntajes_finales = r.first;
        // int iters;
        // tie(puntajes_finales,iters) = page_rank_GS("tests/" + test, p[i-1]);

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
            Vector puntajes_finales;
            puntajes_finales = page_rank_EG("tests/" + test, p[i - 1]);
            // int iters;
            // tie(puntajes_finales,iters) = page_rank_Jacobi("tests/" + test, p[i-1]);

            auto final = chrono::high_resolution_clock::now();
            chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;

            tiempos[j] = tiempo_ejecucion1.count();
        }
        std::string res_out = (string) "tiempos_exp/" + test + "EG" + ".exp";

        write_time_results(res_out, tiempos);
    }
}

void correr_test_nuestros()
{

    string output_path = "test_nuestros/med_tiempoEG.txt";
    // archivo de salida que me da el resultado de los test
    ofstream archivo_salida(output_path);
    archivo_salida << "Instancia \t \t \t \t  \t \t \t& n \t \t & m \t \t & Tiempo tardado  \t & test_result \t \t & Error Abs \t \t" << endl;
    double p = 0.75;
    int i = 0;

    for (string test : tests_implementaciones)
    {
        cout << "corriendo test " << test << " con p " << p << endl;

        // Medir tiempos
        auto inicio = chrono::high_resolution_clock::now();

        // Calculo de rankings
        Vector puntajes_finales;
        puntajes_finales = page_rank_EG("test_nuestros/" + test, p);

        // int iters;
        // tie(puntajes_finales,iters) = page_rank_Jacobi("test_nuestros/" + test, p);
        auto final = chrono::high_resolution_clock::now();

        chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;

        string res_path = test + ".out";

        // pair<bool, double> res_test = resultados_tests(res_path, puntajes_finales);

        string archivo = "test_nuestros/" + test;
        ifstream entrada(archivo);
        int n, m;
        entrada >> n >> m;

        archivo_salida << fixed << setprecision(8) << tiempo_ejecucion1.count() / 1000 << " seg \t   & \t" << endl;
        //<< res_output[res_test.first] << "& \t" << res_test.second <<endl;
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
        vector<double> tiemposEG(reps, 0);
        vector<double> tiemposJac(reps, 0);
        vector<double> tiemposGS(reps, 0);

        for (int j = 0; j < reps; j++)
        {
            // Medir tiempos
            auto inicio = chrono::high_resolution_clock::now();
            // Calculo de rankings
            Vector puntajes_finalesEG = page_rank_EG("test_nuestros/" + test, p);

            auto final = chrono::high_resolution_clock::now();
            chrono::duration<double, std::milli> tiempo_ejecucion1 = final - inicio;

            tiemposEG[j] = tiempo_ejecucion1.count();

            inicio = chrono::high_resolution_clock::now();
            Vector puntajes_finalesJac;
            int itersJ;
            tie(puntajes_finalesJac, itersJ) = page_rank_Jacobi("test_nuestros/" + test, p);
            final = chrono::high_resolution_clock::now();
            tiempo_ejecucion1 = final - inicio;

            tiemposJac[j] = tiempo_ejecucion1.count();

            inicio = chrono::high_resolution_clock::now();
            Vector puntajes_finalesGS;
            int itersGS;
            tie(puntajes_finalesGS, itersGS) = page_rank_GS("test_nuestros/" + test, p);
            final = chrono::high_resolution_clock::now();
            tiempo_ejecucion1 = final - inicio;

            tiemposGS[j] = tiempo_ejecucion1.count();
        }
        write_time_results((string) "tiempos_exp/" + test + "EG" + ".exp", tiemposEG);
        write_time_results((string) "tiempos_exp/" + test + "Jac" + ".exp", tiemposJac);
        write_time_results((string) "tiempos_exp/" + test + "GS" + ".exp", tiemposGS);
    }
}

void correr_test_nuestros_iteraciones(int reps)
{
    double p = 0.75;
    int i = 0;
    for (string test : tests_implementaciones)
    {
        cout << "corriendo test " << test << " con p " << p << endl;

        vector<double> epsilons(reps, 0);

        vector<int> iteracionesJ(reps, 0);
        vector<int> iteracionesGS(reps, 0);

        epsilon = 1e-15;
        for (int j = 0; j < reps; j++)
        {
            epsilon *= 10;
            epsilons[j] = epsilon;
            // Calculo de rankings

            Vector puntajes_finalesJac;
            int iterJ;
            tie(puntajes_finalesJac, iterJ) = page_rank_Jacobi("test_nuestros/" + test, p);
            iteracionesJ[j] = iterJ;

            Vector puntajes_finalesGS;
            int iterGS;
            tie(puntajes_finalesGS, iterGS) = page_rank_GS("test_nuestros/" + test, p);
            iteracionesGS[j] = iterGS;
        }
        write_iters_results((string) "iters_exp/" + test + "Jac" + ".exp", iteracionesJ, epsilons);
        write_iters_results((string) "iters_exp/" + test + "GS" + ".exp", iteracionesGS, epsilons);
        epsilon = 1e-15;
    }
    // reseteo al epsilon original
    epsilon = 1e-10;
}

int main()
{
    correr_test_catedra();
    // correr_test_catedra_experimentacion(1000);
    //correr_test_nuestros_experimentacion(1000);
    // correr_test_nuestros();
    //correr_test_nuestros_iteraciones(14);
    return 0;
}

/* //main entrega
int main(int argc, const char * argv[]){
    string archivo; 
    cin >> archivo; 
    double p;
    cin >> p;
    vector<double> puntajes_finales = calculo_rankings_CSR(archivo, p); 
    string res_path = archivo + ".out"; 
    resultados_tests(res_path, puntajes_finales);
    return 0; 

} */