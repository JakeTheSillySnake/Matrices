#include "matrix_extra.h"

int is_valid(matrix_t *A) {
  int result = 0;
  if (A->rows && A->columns) {
    result = 1;
  }
  return result;
}

int is_square(matrix_t *A) {
  int result = 0;
  if (A->rows == A->columns) {
    result = 1;
  }
  return result;
}

void s21_fill_matrix(double number, matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; number += 1.0, j++)
      A->matrix[i][j] = number;
  }
}

void s21_minor_matrix(int row, int column, matrix_t *A, matrix_t *result) {
  if (A->rows == 1) {
    s21_create_matrix(1, 1, result);
    result->matrix[0][0] = 1;
    return;
  }
  int minor_i = 0;
  int minor_j = 0;
  s21_create_matrix(A->rows - 1, A->columns - 1, result);
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if (i != row && j != column) {
        result->matrix[minor_i][minor_j] = A->matrix[i][j];
        minor_j++;
        if (minor_j == A->columns - 1) {
          minor_j = 0;
          minor_i++;
        }
      }
    }
  }
}