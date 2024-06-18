#ifndef EXTRA_H
#define EXTRA_H

#include "s21_matrix.h"

int is_valid(matrix_t *A);
int is_square(matrix_t *A);
void s21_fill_matrix(double number, matrix_t *A);
void s21_minor_matrix(int row, int column, matrix_t *A, matrix_t *result);

#endif