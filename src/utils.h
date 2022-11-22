#pragma once

#include <fstream> //Para leer el archivo
#include <iomanip> //Para setear la precision de la salida
#include<chrono> //Por el tiempo
#include <cassert> //Para los assert
#include <filesystem> // Para obtener los nombres de los archivos
#include <math.h>
#include <iostream>

#include "types.h"


void print_sparce_matrix(SparseMatrix &m);  
void print_vector(vector<double> const &v);

/*
Genera un archivo .out con los autovalores 
*/
void out_eigvalues(Vector &eigvals, string path); 
void out_eigvectors(Matrix eigvect, string path);

// Lee los test de la catedra  
SparseMatrix read_test(string test_name);
//Compara los test de la catedra con los .out
pair<bool, double> resultados_tests(string res_path, Vector const &v);

void write_time_results(string res_path, vector<double> tiempos);

//normaliza el vector
void normalizar_vector(Vector &v); 

