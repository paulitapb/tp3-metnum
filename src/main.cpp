#include <iostream>

#include "types.h"
#include "utils.h"
#include "eigen.h"

using namespace std;

int main()
{
   SparseMatrix m(5, 5); 
   for(int i = 0; i< 5 ; i++){
        for(int j = 0; j  < 5; j++){
            if(i == j){
                m.coeffRef(i, j) = 2;
            }
             
        }
   }
   m.coeffRef(2,0) = 1;
   m.coeffRef(0,3) = 5;

   SparseMatrix m2(5, 5);
   m2.coeffRef(3, 4) = 17;
   m2.coeffRef(3, 2) = 91;
   m2.coeffRef(0, 1) = 1;
   m2.coeffRef(0, 4) = 32;
   m2.coeffRef(1, 4) = 32; 
   
   print_sparce_matrix(m2);
   elim_gauss(m2, 1e-10); 
   print_sparce_matrix(m2);  

    return 0;
}
