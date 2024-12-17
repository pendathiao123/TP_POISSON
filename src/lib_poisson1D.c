/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
    int diagonals = *la; // Nombre de diagonales
    int cols = *lab;     // Nombre de colonnes 

    //Ligne 0, que des zéros
    // Ligne 1 : Diagonale supérieure (-1)
    for (int j = 1; j < cols; j++) {
        AB[1 * cols + j] = -1.0;
    }

    // Ligne 2 : Diagonale principale (2)
    for (int j = 0; j < cols; j++) {
        AB[2 * cols + j] = 2.0;
    }

    // Ligne 3 : Diagonale inférieure (-1)
    for (int j = 0; j < cols - 1; j++) {
        AB[3 * cols + j] = -1.0;
    }
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
    //EX_SOL = *BC0 + (*BC1-*BC0)*X/L
    
    // Calcul de la longueur du domaine
    double L = X[*la - 1] - X[0];
    if (L <= 0.0) {
        fprintf(stderr, "Erreur : Longueur du domaine (L) invalide : %f\n", L);
        return;
    }

    // Calcul de la solution analytique
    for (int i = 0; i < *la; i++) {
        EX_SOL[i] = (*BC0) + ((*BC1) - (*BC0)) * (X[i] - X[0]) / L;
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
    int i, j;
    double factor;
    int N = *n;
    int KL = *kl;
    int KU = *ku;
    int LAB = *lab;

    for (i = 0; i < N; i++) {
        ipiv[i] = i + 1; // On suppose pas de permutation initialement
    }
    // Factorisation LU
    for (i = 0; i < N - 1; i++) {
        // La pivotisation se fait sur l'élément de la diagonale principale
        if (AB[i + LAB] == 0.0) {
            *info = i + 1;
            return *info; // Pivot nul
        }
        // Calcul du facteur de pivot
        factor = AB[i + LAB + 1] / AB[i + LAB];
        // Mettre à jour la sous-diagonale
        AB[i + LAB + 1] = factor;
        // Mettre à jour la diagonale supérieure
        AB[i + LAB + 2] -= factor * AB[i + LAB + 1];
        // La diagonale inférieure ne change pas
    }
    *info = 0; // Factorisation réussie
    return *info;
}
