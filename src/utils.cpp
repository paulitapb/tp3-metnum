#include "utils.h"

SparseMatrix read_test(string test_name)
{
    string archivo = test_name;

    ifstream entrada(archivo);

    double epsilon = pow(10, -5);

    int n, elem_no_nulos;
    int a, b;

    entrada >> n >> elem_no_nulos;

    std::vector<T> tripletList;
    tripletList.reserve(elem_no_nulos);

    for (int i = 0; i < elem_no_nulos; i++)
    {
        entrada >> a >> b;
        tripletList.push_back(T(b - 1, a - 1, 1));
    }

    SparseMatrix res(n, n);
    res.setFromTriplets(tripletList.begin(), tripletList.end());
    return res;
}

void print_sparce_matrix(SparseMatrix &m)
{
    cout << MatrixXd(m) << endl;
}

void print_vector(vector<double> const &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        cout << v[i] << " ";
    }
    cout << endl;
}

void out_eigvalues(Vector &eigvals, string path)
{
    ofstream output(path);
    for (Vector::iterator it = eigvals.begin(); it != eigvals.end(); it++)
    {
        output << *it << "\n";
    }
}

void out_eigvectors(Matrix eigvect, string path)
{
    ofstream output(path);
    for (int k = 0; k < eigvect.outerSize(); ++k)
    {
        for (Matrix::InnerIterator it(eigvect, k); it; ++it)
        {
            output << it.value() << "\t";
        }
        output << endl;
    }
}

void normalizar_vector(Vector &v)
{

    double suma = 0;
    for (int i = 0; i < v.size(); i++)
    {
        suma += fabs(v(i));
    }
    for (int i = 0; i < v.size(); i++)
    {
        v(i) = v(i) / suma;
    }
}

pair<bool, double> resultados_tests(string res_path, Vector const &v)
{
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
    for (int i = 0; i < v.size(); i++)
    {
        entrada >> res[i];
        salida << fixed << setprecision(6) << v(i) << endl;

        suma_error += fabs(v(i) - res[i]);
        if (fabs(v(i) - res[i]) > 1e-3)
        {
            paso_test &= false;
        }
    }
    return make_pair(paso_test, suma_error / v.size());
}

void write_time_results(string res_path, vector<double> tiempos)
{
    string archivo_out = res_path;
    ofstream salida(archivo_out);

    for (int i = 0; i < tiempos.size(); i++)
    {
        salida << fixed << setprecision(8) << tiempos[i] << endl;
    }
    salida.close();
}

void write_iters_results(string res_path, vector<int> iters, vector<double> es)
{
    string archivo_out = res_path;
    ofstream salida(archivo_out);

    for (int i = 0; i < iters.size(); i++)
    {
        salida << fixed << setprecision(15) << iters[i] << " " << es[i] << endl;
    }
    salida.close();
}