/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *kv)
{
  int ii, jj, kk;
  // Boucle sur chaque colonne de la matrice AB
  for (jj = 0; jj < (*la); jj++)
  {
    kk = jj * (*lab); // Calcul de l'indice de début de la colonne jj dans AB
    // Si kv est supérieur ou égal à 0, initialiser les premiers éléments de la colonne à 0.0
    if (*kv >= 0)
    {
      for (ii = 0; ii < *kv; ii++)
      {
        AB[kk + ii] = 0.0;
      }
    }
    // Définir les valeurs pour la tridiagonale de la matrice du système de Poisson 1D
    AB[kk + *kv] = -1.0;     // Élément de la sous-diagonale
    AB[kk + *kv + 1] = 2.0;  // Élément de la diagonale principale
    AB[kk + *kv + 2] = -1.0; // Élément de la sur-diagonale
  }
  // Gérer les conditions aux frontières pour le premier et le dernier élément
  AB[0] = 0.0; // Premier élément de la première colonne
  // Si kv est égal à 1, définir le deuxième élément de la première colonne à 0
  if (*kv == 1)
  {
    AB[1] = 0.0;
  }
  // Dernier élément de la dernière colonne
  AB[(*lab) * (*la) - 1] = 0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double *AB, int *lab, int *la, int *kv)
{
  int ii, jj, kk;

  // Boucle sur chaque colonne de la matrice AB
  for (jj = 0; jj < (*la); jj++)
  {
    kk = jj * (*lab); // Calcul de l'indice de début de la colonne jj dans AB

    // Si kv est supérieur ou égal à 0, initialiser les premiers éléments de la colonne à 0.0
    if (*kv >= 0)
    {
      for (ii = 0; ii < *kv; ii++)
      {
        AB[kk + ii] = 0.0;
      }
    }

    // Définir les valeurs pour la tridiagonale de l'opérateur identitaire
    AB[kk + *kv] = 0.0;     // Élément de la sous-diagonale
    AB[kk + *kv + 1] = 1.0; // Élément de la diagonale principale
    AB[kk + *kv + 2] = 0.0; // Élément de la sur-diagonale
  }

  // Gérer les conditions aux frontières pour le premier et le dernier élément
  AB[1] = 0.0;                  // Deuxième élément de la première colonne
  AB[(*lab) * (*la) - 1] = 0.0; // Dernier élément de la dernière colonne
}

void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1)
{
  int jj;

  // Définir les valeurs de la condition aux limites pour le premier élément
  RHS[0] = *BC0;

  // Définir les valeurs de la condition aux limites pour le dernier élément
  RHS[(*la) - 1] = *BC1;

  // Initialiser les autres éléments de RHS à 0.0
  for (jj = 1; jj < (*la) - 1; jj++)
  {
    RHS[jj] = 0.0;
  }
}

void set_analytical_solution_DBC_1D(double *EX_SOL, double *X, int *la, double *BC0, double *BC1)
{
  int jj;
  double h, DELTA_T;

  // Calculer la différence entre BC1 et BC0
  DELTA_T = (*BC1) - (*BC0);

  // Boucle sur chaque point de la grille
  for (jj = 0; jj < (*la); jj++)
  {
    // Calculer la solution analytique en chaque point de la grille
    EX_SOL[jj] = (*BC0) + X[jj] * DELTA_T;
  }
}

void set_grid_points_1D(double *x, int *la)
{
  int jj;
  double h;

  // Calculer la distance entre les points de la grille
  h = 1.0 / (1.0 * ((*la) + 1));

  // Boucle sur chaque point de la grille
  for (jj = 0; jj < (*la); jj++)
  {
    // Calculer la position de chaque point de la grille
    x[jj] = (jj + 1) * h;
  }
}

double relative_forward_error(double *x, double *y, int *la)
{
  // Calcul de l'erreur relative
  double norm_diff = 0.0; // Norme de la différence entre x et y
  double norm_y = 0.0;    // Norme de y

  // Boucle pour calculer les normes
  for (int i = 0; i < (*la); i++)
  {
    norm_diff += (x[i] - y[i]) * (x[i] - y[i]); // Ajoute le carré de la différence à la norme de la différence
    norm_y += y[i] * y[i];                      // Ajoute le carré de y[i] à la norme de y
  }

  // Retourne l'erreur relative, qui est la racine carrée de la norme de la différence divisée par la racine carrée de la norme de y
  return sqrt(norm_diff) / sqrt(norm_y);
}

int indexABCol(int i, int j, int *lab)
{
  // Cette fonction calcule l'index linéaire correspondant aux coordonnées (i, j)
  // dans une matrice stockée en colonne majeure (column-major order).
  // i : l'indice de la ligne
  // j : l'indice de la colonne
  // lab : pointeur vers le nombre de lignes de la matrice

  return (j + 1) * (*lab - 1) + i - 1;
}

int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
  int i, j, k;
  double temp;

  // Vérifier les conditions d'entrée
  if (*lab != 4 || *kl != 1 || *ku != 1 || *la < *n)
  {
    *info = -1; // Code d'erreur pour les conditions d'entrée invalides
    return *info;
  }

  *info = 0; // Initialiser le code d'erreur à 0 (pas d'erreur)

  // Boucle principale pour la factorisation LU
  for (k = 0; k < *n - 1; k++)
  {
    int p_i = k;
    double p_v = AB[indexABCol(k, k, lab)];

    // Recherche du pivot dans la colonne k
    for (int i = k + 1; i < *n; i++)
    {
      if (fabs(AB[indexABCol(i, k, lab)]) > fabs(p_v))
      {
        p_i = i;
        p_v = AB[indexABCol(i, k, lab)];
      }
    }

    // Échange des lignes si nécessaire
    if (p_i != k)
    {
      for (j = 0; j < *lab; j++)
      {
        temp = AB[indexABCol(p_i, j, lab)];
        AB[indexABCol(p_i, j, lab)] = AB[indexABCol(k, j, lab)];
        AB[indexABCol(k, j, lab)] = temp;
      }
      ipiv[k] = p_i;
    }
    else
    {
      ipiv[k] = k + 1;
    }

    // Mise à jour de la matrice
    for (i = k + 1; i < *n; i++)
    {
      AB[indexABCol(i, k, lab)] /= AB[indexABCol(k, k, lab)];
      for (j = k + 1; j < *n; j++)
      {
        AB[indexABCol(i, j, lab)] -= AB[indexABCol(i, k, lab)] * AB[indexABCol(k, j, lab)];
      }
    }
  }
  ipiv[*n - 1] = *n; // Dernière permutation

  return *info; // Retourner le code d'erreur (0 si pas d'erreur)
}

double calculate_error(double *SOL, double *EX_SOL, int size)
{
  double erreur = 0.0; // Initialiser l'erreur à 0.0

  // Boucle pour parcourir chaque élément des tableaux SOL et EX_SOL
  for (int i = 0; i < size; ++i)
  {
    // Calculer la différence absolue entre les éléments correspondants de SOL et EX_SOL
    double diff = fabs(SOL[i] - EX_SOL[i]);

    // Ajouter le carré de la différence à l'erreur
    erreur += diff * diff;
  }

  // Retourner la racine carrée de la somme des carrés des différences
  return sqrt(erreur);
}
