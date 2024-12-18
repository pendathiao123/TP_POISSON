/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
if ((*kv)!=0){
  for (int i =0 ; i< (*la) ; i++ ){
    AB[i*(*lab)]  = 0;  

  }
  }

  for (int i =0; i< (*la) ; i++ ){
    AB[i*(*lab)+1] =-1; 

  }
  AB[(*kv)] =0 ; 
  // AB[(*la) ] =0 ; 

  for (int i =0 ; i< (*la) ; i++ ){
    AB[i*(*lab)+2] =2; 

  }

  for (int i =0 /*derniere ligne remplie de -1*/ ; i< (*la)  ; i++ ){
    AB[i*(*lab)+3] =-1; 
  }

  // AB[(*la)*(*kv)+3*(*la)] =0; 
  int last_element = (*lab ) * (*la) -1; 
  AB[last_element] =0; 
}



void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
    int diagonals = *la; // Nombre de diagonales
    int cols = *lab;     // Nombre de colonnes (taille de la matrice)

    for( int j=0; j<cols; j++){
        AB[j+cols]=1.0; //On met tous les élément de la diagonale principale à 1, le reste 0
    }
    
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
    //Cette fonction permet de définir le vecteur b, avec les conditions aux limites
    RHS[0] = *BC0;
    RHS[(*la)-1] = *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
    //EX_SOL = *BC0 + (*BC1-*BC0)*X
    
    int i;
    double T0 = *BC0;
    double T1 = *BC1;
    for (i = 0; i < *la; i++) 
    {
        EX_SOL[i] = X[i] * (T1 - T0) + T0;
    }
    
}  

void set_grid_points_1D(double* x, int* la){
    // Cette fonction permet de définir les point du maillage 
    // On définit dx comme étant la distance entre deux points consécutifs du maillage

    double dx = 1.0/((*la)-1);

    for(int i = 0; i<(*la); i++){
        x[i] = i*dx;
    }
}

double relative_forward_error(double* x, double* y, int* la){
    //Calcul de l'erreur relatif
    double norm_diff = 0.0;
    double norm_y = 0.0;

    for(int i=0; i< (*la); i++){
        norm_diff += (x[i] - y[i]) * (x[i] - y[i]); 
        norm_y +=y[i]*y[i];
    }

    return sqrt(norm_diff)/sqrt(norm_y);
}

int indexABCol(int i, int j, int *lab){
  return 0;
}


//Exercice6
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
    int i, k;
    double factor;
    
    // Vérification des paramètres
    if (*lab != 4 || *kl != 1 || *ku != 1 || *la < *n) {
        *info = -1;
        return *info;
    }

    *info = 0;

    // Initialisation du tableau des pivots
    for (i = 0; i < *n; i++) {
        ipiv[i] = i + 1; // Pas de permutation initialement
    }

    // Factorisation LU
    for (k = 0; k < *n - 1; k++) {
        // Vérification du pivot (diagonale principale)
        if (fabs(AB[2 + k * (*lab)]) < 1e-10) { // Pivot nul
            *info = k + 1;
            return *info;
        }

        // Calcul du facteur de pivot (sous-diagonale / diagonale principale)
        factor = AB[3 + k * (*lab)] / AB[2 + k * (*lab)];

        // Mise à jour de la sous-diagonale (elle contient désormais le facteur)
        AB[3 + k * (*lab)] = factor;

        // Mise à jour de la diagonale principale suivante
        AB[2 + (k + 1) * (*lab)] -= factor * AB[1 + k * (*lab)];
    }

    return *info;
    
}
