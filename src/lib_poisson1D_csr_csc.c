#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Définir une structure pour une matrice CSR (Compressed Sparse Row)
typedef struct {
    int nb_lignes;
    int nb_colonnes;
    int nb_non_zeros;
    double* valeurs;
    int* indices_colonnes;
    int* pointeurs_lignes;
} CSRMatrice;

// Définir une structure pour une matrice CSC (Compressed Sparse Column)
typedef struct {
    int nb_lignes;
    int nb_colonnes;
    int nb_non_zeros;
    double* valeurs;
    int* indices_lignes;
    int* pointeurs_colonnes;
} CSCMatrice;

// Convertir une matrice CSR en matrice dense en ordre de colonne majeure
double* csr_vers_dense_col_major(CSRMatrice* csr) {
    double* dense = (double*)calloc(csr->nb_lignes * csr->nb_colonnes, sizeof(double));
    for (int i = 0; i < csr->nb_lignes; i++) {
        // Remplir une ligne
        for (int j = csr->pointeurs_lignes[i]; j < csr->pointeurs_lignes[i + 1]; j++) {
            dense[i + csr->indices_colonnes[j] * csr->nb_lignes] = csr->valeurs[j];
        }
    }
    return dense;
}

// Convertir une matrice dense en ordre de colonne majeure en matrice CSR
CSRMatrice dense_col_major_vers_csr(double* dense, int nb_lignes, int nb_colonnes) {
    CSRMatrice csr;
    csr.nb_lignes = nb_lignes;
    csr.nb_colonnes = nb_colonnes;
    csr.nb_non_zeros = 0;

    // Compter le nombre d'éléments non nuls
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_colonnes; j++) {
            if (dense[i + j * nb_lignes] != 0) {
                csr.nb_non_zeros++;
            }
        }
    }

    // Allouer la mémoire pour les tableaux CSR
    csr.valeurs = (double*)malloc(csr.nb_non_zeros * sizeof(double));
    csr.indices_colonnes = (int*)malloc(csr.nb_non_zeros * sizeof(int));
    csr.pointeurs_lignes = (int*)malloc((nb_lignes + 1) * sizeof(int));

    // Remplir les tableaux CSR
    int indice_non_zero = 0;
    for (int i = 0; i < nb_lignes; i++) {
        csr.pointeurs_lignes[i] = indice_non_zero;
        for (int j = 0; j < nb_colonnes; j++) {
            if (dense[i + j * nb_lignes] != 0) {
                csr.valeurs[indice_non_zero] = dense[i + j * nb_lignes];
                csr.indices_colonnes[indice_non_zero] = j;
                indice_non_zero++;
            }
        }
    }
    csr.pointeurs_lignes[nb_lignes] = indice_non_zero;

    return csr;
}

// Afficher une matrice CSR
void afficher_matrice_csr(CSRMatrice* csr) {
    for (int i = 0; i < csr->nb_lignes; i++) {
        for (int j = 0; j < csr->nb_colonnes; j++) {
            // Chercher si le nombre existe dans la matrice (0 sinon)
            int trouve = 0;
            for (int k = csr->pointeurs_lignes[i]; k < csr->pointeurs_lignes[i + 1]; k++) {
                if (csr->indices_colonnes[k] == j) {
                    printf("%.4lf\t", csr->valeurs[k]);
                    trouve = 1;
                    break;
                }
            }
            if (!trouve) {
                printf("%.4lf\t", 0.0);
            }
        }
        printf("\n");
    }
}

// Générer une matrice CSR à partir d'une matrice tridiagonale
CSRMatrice csr_de_tridiagonale(double* diag, double* inferieur, double* superieur, int nb_elements) {
    CSRMatrice csr;
    csr.nb_lignes = nb_elements;
    csr.nb_colonnes = nb_elements;
    csr.nb_non_zeros = 3 * nb_elements - 2; // Toutes les lignes ont 3 éléments sauf la première et la dernière avec 2 éléments

    // Allouer la mémoire pour les tableaux CSR
    csr.valeurs = (double*)malloc(csr.nb_non_zeros * sizeof(double));
    csr.indices_colonnes = (int*)malloc(csr.nb_non_zeros * sizeof(int));
    csr.pointeurs_lignes = (int*)malloc((nb_elements + 1) * sizeof(int));

    // Première ligne avec 2 et -1
    csr.valeurs[0] = diag[0];
    csr.indices_colonnes[0] = 0;
    csr.valeurs[1] = superieur[0];
    csr.indices_colonnes[1] = 1;
    csr.pointeurs_lignes[0] = 0;
    csr.pointeurs_lignes[1] = 2;

    // Lignes intermédiaires avec -1, 2, -1
    for (int i = 1; i < nb_elements - 1; i++) {
        csr.valeurs[3 * i - 1] = inferieur[i - 1];
        csr.indices_colonnes[3 * i - 1] = i - 1;
        csr.valeurs[3 * i] = diag[i];
        csr.indices_colonnes[3 * i] = i;
        csr.valeurs[3 * i + 1] = superieur[i];
        csr.indices_colonnes[3 * i + 1] = i + 1;
        csr.pointeurs_lignes[i + 1] = 3 * i + 2;
    }

    // Dernière ligne avec -1 et 2
    csr.valeurs[3 * nb_elements - 3] = inferieur[nb_elements - 2];
    csr.indices_colonnes[3 * nb_elements - 3] = nb_elements - 2;
    csr.valeurs[3 * nb_elements - 2] = diag[nb_elements - 1];
    csr.indices_colonnes[3 * nb_elements - 2] = nb_elements - 1;
    csr.pointeurs_lignes[nb_elements] = 3 * nb_elements - 2;

    return csr;
}

// Convertir une matrice CSC en matrice dense en ordre de colonne majeure
double* csc_vers_dense_col_major(CSCMatrice* csc) {
    double* dense = (double*)malloc(csc->nb_lignes * csc->nb_colonnes * sizeof(double));
    memset(dense, 0, csc->nb_lignes * csc->nb_colonnes * sizeof(double));
    for (int i = 0; i < csc->nb_colonnes; i++) {
        // Remplir une colonne
        for (int j = csc->pointeurs_colonnes[i]; j < csc->pointeurs_colonnes[i + 1]; j++) {
            dense[csc->indices_lignes[j] + i * csc->nb_lignes] = csc->valeurs[j];
        }
    }
    return dense;
}

// Convertir une matrice dense en ordre de colonne majeure en matrice CSC
CSCMatrice dense_col_major_vers_csc(double* dense, int nb_lignes, int nb_colonnes) {
    CSCMatrice csc;
    csc.nb_lignes = nb_lignes;
    csc.nb_colonnes = nb_colonnes;

    // Compter le nombre d'éléments non nuls
    csc.nb_non_zeros = 0;
    for (int i = 0; i < nb_colonnes; i++) {
        for (int j = 0; j < nb_lignes; j++) {
            if (dense[j + i * nb_lignes] != 0) {
                csc.nb_non_zeros++;
            }
        }
    }

    // Allouer la mémoire pour les tableaux CSC
    csc.valeurs = (double*)malloc(csc.nb_non_zeros * sizeof(double));
    csc.indices_lignes = (int*)malloc(csc.nb_non_zeros * sizeof(int));
    csc.pointeurs_colonnes = (int*)malloc((nb_colonnes + 1) * sizeof(int));

    // Remplir les tableaux CSC
    int indice_non_zero = 0;
    for (int i = 0; i < nb_colonnes; i++) {
        csc.pointeurs_colonnes[i] = indice_non_zero;
        for (int j = 0; j < nb_lignes; j++) {
            if (dense[j + i * nb_lignes] != 0) {
                csc.valeurs[indice_non_zero] = dense[j + i * nb_lignes];
                csc.indices_lignes[indice_non_zero] = j;
                indice_non_zero++;
            }
        }
    }
    csc.pointeurs_colonnes[nb_colonnes] = indice_non_zero;

    return csc;
}

// Afficher une matrice CSC
void afficher_matrice_csc(CSCMatrice* csc) {
    for (int i = 0; i < csc->nb_lignes; i++) {
        for (int j = 0; j < csc->nb_colonnes; j++) {
            // Chercher si le nombre existe dans la matrice (0 sinon)
            int trouve = 0;
            for (int k = csc->pointeurs_colonnes[j]; k < csc->pointeurs_colonnes[j + 1]; k++) {
                if (csc->indices_lignes[k] == i) {
                    printf("%lf\t", csc->valeurs[k]);
                    trouve = 1;
                    break;
                }
            }
            if (!trouve) {
                printf("%lf\t", 0.0);
            }
        }
        printf("\n");
    }
}

// Générer une matrice CSC à partir d'une matrice tridiagonale
CSCMatrice csc_de_tridiagonale(double* diag, double* inferieur, double* superieur, int nb_elements) {
    CSCMatrice csc;
    csc.nb_lignes = nb_elements;
    csc.nb_colonnes = nb_elements;
    csc.nb_non_zeros = 3 * nb_elements - 2; // Toutes les colonnes ont 3 éléments sauf la première et la dernière avec 2 éléments

    // Allouer la mémoire pour les tableaux CSC
    csc.valeurs = (double*)malloc(csc.nb_non_zeros * sizeof(double));
    csc.indices_lignes = (int*)malloc(csc.nb_non_zeros * sizeof(int));
    csc.pointeurs_colonnes = (int*)malloc((nb_elements + 1) * sizeof(int));

    // Première colonne avec 2 et -1
    csc.valeurs[0] = diag[0];
    csc.indices_lignes[0] = 0;
    csc.valeurs[1] = superieur[0];
    csc.indices_lignes[1] = 1;
    csc.pointeurs_colonnes[0] = 0;
    csc.pointeurs_colonnes[1] = 2;

    // Colonnes intermédiaires avec -1, 2, -1
    for (int i = 1; i < nb_elements - 1; i++) {
        csc.valeurs[3 * i - 1] = inferieur[i - 1];
        csc.indices_lignes[3 * i - 1] = i - 1;
        csc.valeurs[3 * i] = diag[i];
        csc.indices_lignes[3 * i] = i;
        csc.valeurs[3 * i + 1] = superieur[i];
        csc.indices_lignes[3 * i + 1] = i + 1;
        csc.pointeurs_colonnes[i + 1] = 3 * i + 2;
    }

    // Dernière colonne avec -1 et 2
    csc.valeurs[3 * nb_elements - 3] = inferieur[nb_elements - 2];
    csc.indices_lignes[3 * nb_elements - 3] = nb_elements - 2;
    csc.valeurs[3 * nb_elements - 2] = diag[nb_elements - 1];
    csc.indices_lignes[3 * nb_elements - 2] = nb_elements - 1;
    csc.pointeurs_colonnes[nb_elements] = 3 * nb_elements - 2;

    return csc;
}

// Cloner une matrice CSR
CSRMatrice cloner_csr(CSRMatrice* csr) {
    CSRMatrice clone;
    clone.nb_lignes = csr->nb_lignes;
    clone.nb_colonnes = csr->nb_colonnes;
    clone.nb_non_zeros = csr->nb_non_zeros;
    clone.valeurs = (double*)malloc(clone.nb_non_zeros * sizeof(double));
    clone.indices_colonnes = (int*)malloc(clone.nb_non_zeros * sizeof(int));
    clone.pointeurs_lignes = (int*)malloc((clone.nb_lignes + 1) * sizeof(int));
    memcpy(clone.valeurs, csr->valeurs, clone.nb_non_zeros * sizeof(double));
    memcpy(clone.indices_colonnes, csr->indices_colonnes, clone.nb_non_zeros * sizeof(int));
    memcpy(clone.pointeurs_lignes, csr->pointeurs_lignes, (clone.nb_lignes + 1) * sizeof(int));
    return clone;
}

// Cloner une matrice CSC
CSCMatrice cloner_csc(CSCMatrice* csc) {
    CSCMatrice clone;
    clone.nb_lignes = csc->nb_lignes;
    clone.nb_colonnes = csc->nb_colonnes;
    clone.nb_non_zeros = csc->nb_non_zeros;
    clone.valeurs = (double*)malloc(clone.nb_non_zeros * sizeof(double));
    clone.indices_lignes = (int*)malloc(clone.nb_non_zeros * sizeof(int));
    clone.pointeurs_colonnes = (int*)malloc((clone.nb_colonnes + 1) * sizeof(int));
    memcpy(clone.valeurs, csc->valeurs, clone.nb_non_zeros * sizeof(double));
    memcpy(clone.indices_lignes, csc->indices_lignes, clone.nb_non_zeros * sizeof(int));
    memcpy(clone.pointeurs_colonnes, csc->pointeurs_colonnes, (clone.nb_colonnes + 1) * sizeof(int));
    return clone;
}

// Libérer la mémoire d'une matrice CSR
void liberer_csr(CSRMatrice* csr) {
    free(csr->valeurs);
    free(csr->indices_colonnes);
    free(csr->pointeurs_lignes);
}

// Libérer la mémoire d'une matrice CSC
void liberer_csc(CSCMatrice* csc) {
    free(csc->valeurs);
    free(csc->indices_lignes);
    free(csc->pointeurs_colonnes);
}

// Obtenir un élément d'une matrice CSR
double element_csr_a(CSRMatrice* csr, int i, int j) {
    for (int indice_non_zero = csr->pointeurs_lignes[i]; indice_non_zero < csr->pointeurs_lignes[i + 1]; indice_non_zero++) {
        if (csr->indices_colonnes[indice_non_zero] == j) {
            return csr->valeurs[indice_non_zero];
        }
    }

    // Si l'élément n'est pas trouvé, c'est un 0
    return 0;
}

// Écrire un élément dans une matrice CSR
void ecrire_csr_a(CSRMatrice* csr, int i, int j, double valeur) {
    for (int indice_non_zero = csr->pointeurs_lignes[i]; indice_non_zero < csr->pointeurs_lignes[i + 1]; indice_non_zero++) {
        if (csr->indices_colonnes[indice_non_zero] == j) {
            csr->valeurs[indice_non_zero] = valeur;
            return;
        }
    }
}

// Obtenir un élément d'une matrice CSC
double element_csc_a(CSCMatrice* csc, int i, int j) {
    for (int indice_non_zero = csc->pointeurs_colonnes[j]; indice_non_zero < csc->pointeurs_colonnes[j + 1]; indice_non_zero++) {
        if (csc->indices_lignes[indice_non_zero] == i) {
            return csc->valeurs[indice_non_zero];
        }
    }

    // Si l'élément n'est pas trouvé, c'est un 0
    return 0;
}

// Écrire un élément dans une matrice CSC
void ecrire_csc_a(CSCMatrice* csc, int i, int j, double valeur) {
    for (int indice_non_zero = csc->pointeurs_colonnes[j]; indice_non_zero < csc->pointeurs_colonnes[j + 1]; indice_non_zero++) {
        if (csc->indices_lignes[indice_non_zero] == i) {
            csc->valeurs[indice_non_zero] = valeur;
            return;
        }
    }
}

// Générer une matrice CSR à partir d'une diagonale
CSRMatrice csr_de_diagonale(double* diag, int nb_elements) {
    CSRMatrice csr;
    csr.nb_lignes = nb_elements;
    csr.nb_colonnes = nb_elements;
    csr.nb_non_zeros = nb_elements;

    // Allouer la mémoire pour les tableaux CSR
    csr.valeurs = (double*)malloc(csr.nb_non_zeros * sizeof(double));
    csr.indices_colonnes = (int*)malloc(csr.nb_non_zeros * sizeof(int));
    csr.pointeurs_lignes = (int*)malloc((nb_elements + 1) * sizeof(int));

    // Tous les éléments sur la diagonale sont également dans le même ordre dans le tableau des valeurs CSR
    memcpy(csr.valeurs, diag, nb_elements * sizeof(double));

    // Remplir les tableaux CSR
    for (int i = 0; i < nb_elements; i++) {
        csr.indices_colonnes[i] = i;
        csr.pointeurs_lignes[i] = i;
    }
    csr.pointeurs_lignes[nb_elements] = nb_elements;

    return csr;
}

// Générer une matrice CSC à partir d'une diagonale
CSCMatrice csc_de_diagonale(double* diag, int nb_elements) {
    CSCMatrice csc;
    csc.nb_lignes = nb_elements;
    csc.nb_colonnes = nb_elements;
    csc.nb_non_zeros = nb_elements;

    // Allouer la mémoire pour les tableaux CSC
    csc.valeurs = (double*)malloc(csc.nb_non_zeros * sizeof(double));
    csc.indices_lignes = (int*)malloc(csc.nb_non_zeros * sizeof(int));
    csc.pointeurs_colonnes = (int*)malloc((nb_elements + 1) * sizeof(int));

    // Remplir les tableaux CSC
    for (int i = 0; i < nb_elements; i++) {
        csc.valeurs[i] = diag[i];
        csc.indices_lignes[i] = i;
        csc.pointeurs_colonnes[i] = i;
    }
    csc.pointeurs_colonnes[nb_elements] = nb_elements;

    return csc;
}

// Générer une matrice CSR à partir d'une matrice tridiagonale
CSRMatrice csr_de_tridiagonale(double* diag, double* inferieur, double* superieur, int nb_elements) {
    CSRMatrice csr;
    csr.nb_lignes = nb_elements;
    csr.nb_colonnes = nb_elements;
    csr.nb_non_zeros = 3 * nb_elements;  // Toutes les lignes ont 3 éléments
    csr.valeurs = (double*)malloc(csr.nb_non_zeros * sizeof(double));
    csr.indices_colonnes = (int*)malloc(csr.nb_non_zeros * sizeof(int));
    csr.pointeurs_lignes = (int*)malloc((nb_elements + 1) * sizeof(int));

    // Copier ligne par ligne
    for (int i = 0; i < nb_elements; i++) {
        csr.pointeurs_lignes[i] = 3 * i;
        csr.valeurs[3 * i] = diag[i];
        csr.indices_colonnes[3 * i] = i;
        if (i > 0) {
            csr.valeurs[3 * i - 1] = inferieur[i - 1];
            csr.indices_colonnes[3 * i - 1] = i - 1;
        }
        if (i < nb_elements - 1) {
            csr.valeurs[3 * i + 1] = superieur[i];
            csr.indices_colonnes[3 * i + 1] = i + 1;
        }
    }
    csr.pointeurs_lignes[nb_elements] = 3 * nb_elements;

    return csr;
}

// Générer une matrice CSC à partir d'une matrice tridiagonale
CSCMatrice csc_de_tridiagonale(double* diag, double* inferieur, double* superieur, int nb_elements) {
    CSCMatrice csc;
    csc.nb_lignes = nb_elements;
    csc.nb_colonnes = nb_elements;
    csc.nb_non_zeros = 3 * nb_elements;  // Toutes les colonnes ont 3 éléments
    csc.valeurs = (double*)malloc(csc.nb_non_zeros * sizeof(double));
    csc.indices_lignes = (int*)malloc(csc.nb_non_zeros * sizeof(int));
    csc.pointeurs_colonnes = (int*)malloc((nb_elements + 1) * sizeof(int));

    // Copier colonne par colonne
    for (int i = 0; i < nb_elements; i++) {
        csc.pointeurs_colonnes[i] = 3 * i;
        csc.valeurs[3 * i] = diag[i];
        csc.indices_lignes[3 * i] = i;
        if (i > 0) {
            csc.valeurs[3 * i - 1] = inferieur[i - 1];
            csc.indices_lignes[3 * i - 1] = i - 1;
        }
        if (i < nb_elements - 1) {
            csc.valeurs[3 * i + 1] = superieur[i];
            csc.indices_lignes[3 * i + 1] = i + 1;
        }
    }
    csc.pointeurs_colonnes[nb_elements] = 3 * nb_elements;

    return csc;
}

// Générer une matrice CSR à partir d'une matrice triangulaire inférieure
CSRMatrice csr_de_triangulaire_inferieure(double* mat, int nb_lignes) {
    CSRMatrice csr;
    csr.nb_lignes = nb_lignes;
    csr.nb_colonnes = nb_lignes;
    csr.nb_non_zeros = nb_lignes * (nb_lignes + 1) / 2;

    // Allouer la mémoire pour les tableaux CSR
    csr.valeurs = (double*)malloc(csr.nb_non_zeros * sizeof(double));
    csr.indices_colonnes = (int*)malloc(csr.nb_non_zeros * sizeof(int));
    csr.pointeurs_lignes = (int*)malloc((nb_lignes + 1) * sizeof(int));

    // mat est une matrice triangulaire inférieure stockée en 1D, tous les éléments sont dans le bon ordre
    memcpy(csr.valeurs, mat, csr.nb_non_zeros * sizeof(double));

    // Calculer les indices des colonnes
    int non_zero_idx = 0;
    int current_idx = 0;
    for (int i = 0; i < nb_lignes; i++) {
        csr.pointeurs_lignes[i] = non_zero_idx;
        for (int j = 0; j <= i; j++) {
            csr.indices_colonnes[non_zero_idx] = j;
            non_zero_idx++;
            current_idx++;
        }
    }
    csr.pointeurs_lignes[nb_lignes] = non_zero_idx;

    return csr;
}

// Générer une matrice CSC à partir d'une matrice triangulaire inférieure
CSCMatrice csc_de_triangulaire_inferieure(double* mat, int nb_lignes) {
    CSCMatrice csc;
    csc.nb_lignes = nb_lignes;
    csc.nb_colonnes = nb_lignes;
    csc.nb_non_zeros = nb_lignes * (nb_lignes + 1) / 2;

    // Allouer la mémoire pour les tableaux CSC
    csc.valeurs = (double*)malloc(csc.nb_non_zeros * sizeof(double));
    csc.indices_lignes = (int*)malloc(csc.nb_non_zeros * sizeof(int));
    csc.pointeurs_colonnes = (int*)malloc((nb_lignes + 1) * sizeof(int));

    // Calculer les pointeurs de colonnes et les valeurs
    // mat est une matrice triangulaire inférieure, mais nous ne pouvons pas copier car elle est dans le mauvais ordre
    int non_zero_idx = 0;
    int current_idx = 0;
    for (int i = 0; i < nb_lignes; i++) {
        csc.pointeurs_colonnes[i] = non_zero_idx;
        for (int j = 0; j <= i; j++) {
            csc.valeurs[non_zero_idx] = mat[current_idx];
            csc.indices_lignes[non_zero_idx] = j;
            non_zero_idx++;
            current_idx++;
        }
    }
    csc.pointeurs_colonnes[nb_lignes] = non_zero_idx;

    return csc;
}