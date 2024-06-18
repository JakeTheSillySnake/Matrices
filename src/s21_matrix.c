#include "s21_matrix.h"

#include "matrix_extra.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  result->matrix = NULL;
  int output = OK;

  if (rows && columns) {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));
  }
  if (result->matrix) {
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = calloc(columns, sizeof(double));
      if (!result->matrix[i]) {
        output = INCORRECT_MATRIX;
      }
    }
  } else {
    output = INCORRECT_MATRIX;
  }
  return output;
}

void s21_remove_matrix(matrix_t *A) {
  for (int i = 0; i < A->rows; i++) {
    free(A->matrix[i]);
  }
  free(A->matrix);
  A->rows = 0;
  A->columns = 0;
  A->matrix = NULL;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int output = SUCCESS;
  if (A->rows != B->rows || A->columns != B->columns || !is_valid(A) ||
      !is_valid(B)) {
    output = FAILURE;
  } else {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) >= 1e-7) {
          output = FAILURE;
        }
      }
    }
  }
  return output;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int output = OK;
  if (A->rows == B->rows && A->columns == B->columns) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else {
      output = INCORRECT_MATRIX;
    }
  } else {
    output = CALCULATION_ERROR;
  }
  return output;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int output = OK;
  if (A->rows == B->rows && A->columns == B->columns) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else {
      output = INCORRECT_MATRIX;
    }
  } else {
    output = CALCULATION_ERROR;
  }
  return output;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int output = OK;
  if (s21_create_matrix(A->rows, A->columns, result) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  } else {
    output = INCORRECT_MATRIX;
  }
  return output;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int output = OK;
  if (A->columns != B->rows) {
    return CALCULATION_ERROR;
  }
  if (is_valid(A) && is_valid(B) &&
      s21_create_matrix(A->rows, B->columns, result) == OK) {
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        for (int k = 0; k < A->columns; k++)
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
      }
    }
  } else {
    output = INCORRECT_MATRIX;
  }
  return output;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int output = OK;
  if (s21_create_matrix(A->columns, A->rows, result) == OK) {
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  } else {
    output = INCORRECT_MATRIX;
  }
  return output;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int output = OK;
  if (is_square(A)) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          matrix_t minor_matrix = {0};
          double det = 0;
          s21_minor_matrix(i, j, A, &minor_matrix);
          s21_determinant(&minor_matrix, &det);
          s21_remove_matrix(&minor_matrix);
          result->matrix[i][j] = det * ((i + j) % 2 == 0 ? 1 : -1);
        }
      }
    } else {
      output = INCORRECT_MATRIX;
    }
  } else {
    output = CALCULATION_ERROR;
  }
  return output;
}

int s21_determinant(matrix_t *A, double *result) {
  int output = OK;
  if (is_square(A) && is_valid(A)) {
    *result = 0;
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    } else {
      int sign = 1;
      for (int i = 0; i < A->rows; i++) {
        matrix_t minor_matrix = {0};
        double det = 0;
        s21_minor_matrix(0, i, A, &minor_matrix);
        s21_determinant(&minor_matrix, &det);
        s21_remove_matrix(&minor_matrix);
        *result += A->matrix[0][i] * det * sign;
        sign *= -1;
      }
    }
  } else if (!is_valid(A)) {
    output = INCORRECT_MATRIX;
  } else {
    output = CALCULATION_ERROR;
  }
  return output;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int output = OK;
  if (is_square(A)) {
    if (A->rows == 1 && A->columns == 1) {
      s21_create_matrix(A->rows, A->columns, result);
      result->matrix[0][0] = 1 / A->matrix[0][0];
    } else if (is_valid(A)) {
      double det;
      s21_determinant(A, &det);
      if (fabs(det) >= 1e-7) {
        matrix_t comps, comps_transposed;
        if (s21_calc_complements(A, &comps) == OK) {
          s21_transpose(&comps, &comps_transposed);
          s21_mult_number(&comps_transposed, 1 / det, result);
          s21_remove_matrix(&comps_transposed);
        }
        s21_remove_matrix(&comps);
      } else {
        output = CALCULATION_ERROR;
      }
    } else {
      output = INCORRECT_MATRIX;
    }
  } else {
    output = CALCULATION_ERROR;
  }
  return output;
}