#include <stdio.h>

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv) {
    
    int n = *la;  
    int rows = *lab; 
    int diag = *kv;
   
    for (int i = 0; i < rows * n; i++) {
        AB[i] = 0.0;
    }

    for (int i = 0; i < n; i++) {
        AB[diag * n + i] = 2.0;  
    }
    for (int i = 0; i < n - 1; i++) {
        AB[(diag - 1) * n + i + 1] = -1.0;  
    }
    for (int i = 0; i < n - 1; i++) {
        AB[(diag + 1) * n + i] = -1.0;  
}}

int main() {
    int lab = 4;   
    int la = 15;  
    int kv = 1;   
    double AB[lab* la];  

    // Appel de la fonction pour remplir la matrice
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

    // Affichage des résultats pour vérifier le stockage
    printf("Les éléments de la matrice stockée en GB sont :\n");
    for (int i = 0; i < lab; i++) {
        for (int j = 0; j < la; j++) {
            printf("%6.2f ", AB[i * la + j]);
        }
        printf("\n");
    }

    return 0;
}
