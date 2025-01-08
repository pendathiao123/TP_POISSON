/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void eig_poisson1D(double *eigval, int *la)
{
  // Calcul de la taille de l'intervalle
  double h = 1.0 / (*la + 1);

  // Boucle pour calculer les valeurs propres
  for (int k = 1; k <= *la; k++)
  {
    // Calcul de sin(kπh / 2)
    double sin_theta = sin(k * M_PI * h / 2.0);
    // Calcul de la valeur propre λ^k = 4 sin²(kπh / 2)
    eigval[k - 1] = 4.0 * sin_theta * sin_theta;
  }
}

// Fonction pour calculer les valeurs propres d'une matrice
// La fonction prend en entrée un pointeur vers un entier représentant la taille de la matrice (la)
// et retourne un tableau de doubles contenant les valeurs propres.

double *compute_eigen_values(int *la)
{
  // Allouer de la mémoire pour les valeurs propres (eigval)
  // La taille du tableau est donnée par la valeur pointée par 'la'
  double *eigval = (double *)malloc(*la * sizeof(double));

  // Calculer les valeurs propres en appelant la fonction eig_poisson1D
  // Cette fonction remplit le tableau eigval avec les valeurs propres de la matrice
  eig_poisson1D(eigval, la);

  // Retourner le tableau de valeurs propres
  return eigval;
}

double eigmax_poisson1D(int *la)
{
  // Calculer les valeurs propres
  double *eigval = compute_eigen_values(la);

  // Trouver la valeur propre maximale
  double eigmax = eigval[0]; // Initialiser avec la première valeur propre
  for (int i = 1; i < *la; i++)
  {
    if (eigval[i] > eigmax)
    {
      eigmax = eigval[i]; // Mettre à jour eigmax si une valeur plus grande est trouvée
    }
  }

  // Libérer la mémoire allouée pour les valeurs propres
  free(eigval);

  // Retourner la valeur propre maximale
  return eigmax;
}

double eigmin_poisson1D(int *la)
{
  // Calculer les valeurs propres
  double *eigval = compute_eigen_values(la);

  // Trouver la valeur propre minimale
  double eigmin = eigval[0]; // Initialiser avec la première valeur propre
  for (int i = 1; i < *la; i++)
  {
    if (eigval[i] < eigmin)
    {
      eigmin = eigval[i]; // Mettre à jour eigmin si une valeur plus petite est trouvée
    }
  }

  // Libérer la mémoire allouée pour les valeurs propres
  free(eigval);

  // Retourner la valeur propre minimale
  return eigmin;
}

double richardson_alpha_opt(int *la)
{
  // Calculer l'alpha optimal pour la méthode de Richardson
  // en utilisant les valeurs propres minimale et maximale
  return 2.0 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  double alpha = *alpha_rich;      // Valeur de alpha pour la méthode de Richardson
  int nb_rows = *lab;              // Nombre de lignes de la matrice AB
  int nb_cols = *la;               // Nombre de colonnes de la matrice AB
  int upper_diags = *ku;           // Nombre de diagonales supérieures
  int lower_diags = *kl;           // Nombre de diagonales inférieures
  int max_iters = *maxit;          // Nombre maximal d'itérations
  double tolerance = *tol;         // Tolérance pour la convergence
  int *nb_iters_final = nbite;     // Nombre final d'itérations
  double *residual_norms = resvec; // Vecteur pour stocker les normes résiduelles

  // Allocation mémoire pour un vecteur temporaire
  double *tmp_vec = malloc(sizeof(double) * nb_cols);
  if (tmp_vec == NULL)
  {
    fprintf(stderr, "Memory allocation failed for temp\n");
    exit(EXIT_FAILURE);
  }

  // tmp_vec = RHS
  memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

  // tmp_vec = tmp_vec - AB * X
  cblas_dgbmv(CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_rows, X, 1, 1.0, tmp_vec, 1);

  // Calcul de la norme du résidu initial
  double norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
  residual_norms[0] = norm_res;

  // Itérer jusqu'à atteindre le nombre maximal d'itérations ou jusqu'à ce que la tolérance soit atteinte
  int iter;
  for (iter = 1; iter < max_iters; ++iter)
  {
    // X = X + alpha * tmp_vec
    cblas_daxpy(nb_cols, alpha, tmp_vec, 1, X, 1);

    // tmp_vec = RHS
    memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

    // tmp_vec = tmp_vec - AB * X
    cblas_dgbmv(CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_rows, X, 1, 1.0, tmp_vec, 1);

    // Calcul de la norme du résidu et stockage
    norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
    residual_norms[iter] = norm_res;

    // Vérification de la convergence
    if (norm_res < tolerance)
    {
      break;
    }
  }

  // Mise à jour du nombre final d'itérations
  *nb_iters_final = iter;

  // Libération de la mémoire allouée pour le vecteur temporaire
  free(tmp_vec);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv)
{
  int nb_rows = *lab;    // Nombre de lignes de la matrice AB
  int nb_cols = *la;     // Nombre de colonnes de la matrice AB
  int upper_diags = *ku; // Nombre de diagonales supérieures
  int lower_diags = *kl; // Nombre de diagonales inférieures
  int diagonals = *kv;   // Nombre total de diagonales

  // Tous les éléments de MB sont initialisés à 0 sauf la diagonale
  memset(MB, 0, nb_rows * diagonals * sizeof(double));

  // Puisque D est une matrice diagonale, nous pouvons simplement inverser les éléments diagonaux
  for (int i = 0; i < nb_cols; i++)
  {
    MB[1 + i * (diagonals + upper_diags + lower_diags)] = 1.0 / AB[1 + i * (diagonals + upper_diags + lower_diags)];
  }
}

int calculate_lower_triangular_index(int i, int j, int n)
{
  // Cette fonction calcule l'index linéaire correspondant aux coordonnées (i, j)
  // dans une matrice triangulaire inférieure stockée sous forme linéaire.
  return (i + 1) * i / 2 + j;
}

void calculate_inverse_lower_triangular(double *lower_triangular, double *inv_lower_triangular, int n)
{
  // Mettre les éléments diagonaux à leur réciproque
  for (int i = 0; i < n; i++)
  {
    inv_lower_triangular[calculate_lower_triangular_index(i, i, n)] = 1.0 / lower_triangular[calculate_lower_triangular_index(i, i, n)];
  }

  // Maintenant, nous calculons itérativement les autres éléments de haut en bas de sorte que A * inv_A = I
  for (int i = 1; i < n; i++)
  {
    // Calculer les éléments de la i-ème ligne
    for (int j = 0; j < i; j++)
    {
      double sum = 0.0;
      for (int k = j; k < i; k++)
      {
        sum += lower_triangular[calculate_lower_triangular_index(i, k, n)] *
               inv_lower_triangular[calculate_lower_triangular_index(k, j, n)];
      }
      inv_lower_triangular[calculate_lower_triangular_index(i, j, n)] =
          -sum / lower_triangular[calculate_lower_triangular_index(i, i, n)];
    }
  }
}
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv)
{
  int nb_rows = *lab;    // Nombre de lignes de la matrice AB
  int nb_cols = *la;     // Nombre de colonnes de la matrice AB
  int upper_diags = *ku; // Nombre de diagonales supérieures
  int lower_diags = *kl; // Nombre de diagonales inférieures
  int diagonals = *kv;   // Nombre total de diagonales

  // Passer de tridiagonale à triangulaire inférieure (D - E)
  double *lower_triangular = (double *)malloc(((nb_cols + 1) * nb_cols / 2) * sizeof(double));
  memset(lower_triangular, 0, ((nb_cols + 1) * nb_cols / 2) * sizeof(double));

  // Copier la diagonale
  int current_idx = 0;
  for (int i = 0; i < nb_cols; i++)
  {
    lower_triangular[current_idx] = AB[1 + i * (diagonals + upper_diags + lower_diags)];
    current_idx += i + 2;
  }

  // Copier la diagonale inférieure
  for (int i = 1; i < nb_cols; i++)
  {
    lower_triangular[calculate_lower_triangular_index(i, i - 1, nb_cols)] = AB[0 + i * (diagonals + upper_diags + lower_diags)];
  }

  // Inverser la matrice triangulaire inférieure
  double *inv_lower_triangular = (double *)malloc(((nb_cols + 1) * nb_cols / 2) * sizeof(double));
  calculate_inverse_lower_triangular(lower_triangular, inv_lower_triangular, nb_cols);

  // Tronquer inv_lower_triangular pour revenir à une matrice tridiagonale
  memset(MB, 0, nb_rows * (diagonals + upper_diags + lower_diags) * sizeof(double));
  current_idx = 0;

  // Copier la diagonale
  for (int i = 0; i < nb_cols; i++)
  {
    MB[1 + i * (diagonals + upper_diags + lower_diags)] = inv_lower_triangular[current_idx];
    current_idx += i + 2;
  }

  // Copier la diagonale inférieure
  for (int i = 1; i < nb_cols; i++)
  {
    MB[0 + i * (diagonals + upper_diags + lower_diags)] = inv_lower_triangular[calculate_lower_triangular_index(i, i - 1, nb_cols)];
  }

  // Libérer la mémoire allouée pour lower_triangular
  free(lower_triangular);
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  int nb_rows = *lab;              // Nombre de lignes de la matrice AB
  int nb_cols = *la;               // Nombre de colonnes de la matrice AB
  int upper_diags = *ku;           // Nombre de diagonales supérieures
  int lower_diags = *kl;           // Nombre de diagonales inférieures
  int max_iters = *maxit;          // Nombre maximal d'itérations
  double tolerance = *tol;         // Tolérance pour la convergence
  int *nb_iters_final = nbite;     // Nombre final d'itérations
  double *residual_norms = resvec; // Vecteur pour stocker les normes résiduelles

  // Initialiser X à 0
  memset(X, 0, nb_cols * sizeof(double));

  // Allocation mémoire pour un vecteur temporaire
  double *tmp_vec = malloc(sizeof(double) * nb_cols);
  if (tmp_vec == NULL)
  {
    fprintf(stderr, "Memory allocation failed for temp\n");
    exit(EXIT_FAILURE);
  }

  // tmp_vec = RHS
  memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

  // tmp_vec = tmp_vec - AB * X
  cblas_dgbmv(CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_rows, X, 1, 1.0, tmp_vec, 1);

  // Calculer la norme du résidu initial
  double norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
  residual_norms[0] = norm_res;

  // Itérer jusqu'à atteindre le nombre maximal d'itérations et arrêter si la tolérance est atteinte
  int iter;
  for (iter = 1; iter < max_iters; ++iter)
  {
    // X = X + MB * tmp_vec
    cblas_dgbmv(CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, 1.0, MB, nb_rows, tmp_vec, 1, 1.0, X, 1);

    // tmp_vec = RHS
    memcpy(tmp_vec, RHS, nb_cols * sizeof(double));

    // tmp_vec = tmp_vec - AB * X
    cblas_dgbmv(CblasColMajor, CblasNoTrans, nb_cols, nb_cols, lower_diags, upper_diags, -1.0, AB, nb_rows, X, 1, 1.0, tmp_vec, 1);

    // Calculer la norme du résidu et la stocker
    norm_res = cblas_dnrm2(nb_cols, tmp_vec, 1) / cblas_dnrm2(nb_cols, RHS, 1);
    residual_norms[iter] = norm_res;

    // Vérifier la convergence
    if (norm_res < tolerance)
    {
      break;
    }
  }

  // Mettre à jour le nombre final d'itérations
  *nb_iters_final = iter;

  // Libérer la mémoire allouée pour le vecteur temporaire
  free(tmp_vec);
}
