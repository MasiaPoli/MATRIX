#ifndef MATR_H
#define MATR_H
#include <stdio.h>
#include <stdlib.h>
typedef struct matrix matrix;
typedef double element_t;
struct matrix;
/*��������� ����� �����*/
size_t matrix_m(matrix* a);
/*��������� ����� ��������*/
size_t matrix_n(matrix* a);
/*�������� �������*/
matrix* alloc(size_t M, size_t N);
/*�������� ������� �������*/
matrix* null_alloc(size_t M, size_t N);
/*�������� ��������� �������*/
matrix* e_alloc(size_t M, size_t N);
/*��������� �������� �� �������*/
void change(matrix* a, size_t i, size_t j, element_t x);
/*��������� �������� �� �������*/
element_t poind(matrix* a, size_t i, size_t j);
/*��������� ������� �������*/
matrix* re_alloc(matrix* a, size_t m, size_t n);
/*������������ ������ � �������� �������*/
void matrix_free(matrix* a);
#endif // MATR_H
