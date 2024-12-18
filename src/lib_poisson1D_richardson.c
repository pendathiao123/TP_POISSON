/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
  double h=1.0/(*la+1.0);
  for(int k=0; k<*la; k++){
    eigval[k] = 4.0*pow(sin((k*M_PI*h)/2.0),2);
  }
}

double eigmax_poisson1D(int *la){
  double* eigval = malloc((*la) * sizeof(double)); // Allocation mémoire pour les valeurs propres
    if (eigval == NULL) {
        fprintf(stderr, "Erreur d'allocation mémoire pour eigval\n");
        exit(EXIT_FAILURE);
    }

    // Appel de la fonction pour remplir eigval
    eig_poisson1D(eigval, la);

    // Trouver la valeur maximale
    double max_val = eigval[0];
    for (int i = 1; i < *la; ++i) {
        max_val = fmax(max_val, eigval[i]);
    }

    free(eigval); // Libération de la mémoire
    return max_val;
}

double eigmin_poisson1D(int *la){
 double* eigval = malloc((*la) * sizeof(double)); // Allocation mémoire pour les valeurs propres
    if (eigval == NULL) {
        fprintf(stderr, "Erreur d'allocation mémoire pour eigval\n");
        exit(EXIT_FAILURE);
    }

    // Appel de la fonction pour remplir eigval
    eig_poisson1D(eigval, la);

    // Trouver la valeur maximale
    double max_val = eigval[0];
    for (int i = 1; i < *la; ++i) {
        max_val = fmin(max_val, eigval[i]);
    }

    free(eigval); // Libération de la mémoire
    return max_val;
}

double richardson_alpha_opt(int *la){
  double lamda_max = eigmax_poisson1D(la);
  double lamda_min = eigmin_poisson1D(la);
  
  double alpha = 2.0/(lamda_max+lamda_min);
  return alpha;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  int n = *la;
  double alpha = *alpha_rich;
  double res_norm;
  double *R = malloc(n * sizeof(double)); // Résidu

    for (int i = 0; i < n; i++) {
        R[i] = RHS[i];
        for (int j = 0; j < n; j++) {
            R[i] -= AB[i * n + j] * X[j];
        }
    }

    // Calcul de la norme du résidu (norme L2)
    res_norm = 0.0;
    for (int i = 0; i < n; i++) {
        res_norm += R[i] * R[i]; // Somme des carrés des éléments
    }
    res_norm = sqrt(res_norm); // Racine carrée pour la norme L2
    resvec[0] = res_norm;

    *nbite = 0;
    while (res_norm > *tol && *nbite < *maxit) {
        (*nbite)++;
        
        // Mise à jour de X : X = X + alpha * R
        for (int i = 0; i < n; i++) {
            X[i] += alpha * R[i];
        }

        // Recalcul du résidu : R = RHS - A*X
        for (int i = 0; i < n; i++) {
            R[i] = RHS[i]; // Commence avec RHS
            for (int j = 0; j < n; j++) {
                R[i] -= AB[i * n + j] * X[j]; // Soustraction de A*X
            }
        }

        // Calcul de la norme du résidu
        res_norm = 0.0;
        for (int i = 0; i < n; i++) {
            res_norm += R[i] * R[i]; // Somme des carrés
        }
        res_norm = sqrt(res_norm); // Norme L2
        resvec[*nbite] = res_norm;
    }

    free(R);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

