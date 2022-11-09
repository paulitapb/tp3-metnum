#pragma once

#include <fstream> //Para leer el archivo
#include <iomanip> //Para setear la precision de la salida
#include<chrono> //Por el tiempo
#include <cassert> //Para los assert
#include <filesystem> // Para obtener los nombres de los archivos
#include <math.h>
#include <iostream>

#include "types.h"


SparseMatrix read_matrix_karate(); 
SparseMatrix read_ego_facebook();

//SparseMatrix read_test();

void print_sparce_matrix(SparseMatrix &m);  
 
/*
Genera un archivo .out con los autovalores 
*/
void out_eigvalues(Vector &eigvals, string path); 
void out_eigvectors(Matrix eigvect, string path); 