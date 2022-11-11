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
                m.coeffRef(i, j) = 1;
            }
             
        }
   }
   elim_gauss(m, 1e-10); 
   print_sparce_matrix(m);  

    return 0;
}
