#ifndef MATR_H
#define MATR_H
#include <stdio.h>
#include <stdlib.h>
typedef struct matrix matrix;
typedef double element_t;
struct matrix;
/*Получение числа строк*/
size_t matrix_m(matrix* a);
/*Получение числа столбцов*/
size_t matrix_n(matrix* a);
/*Создание матрицы*/
matrix* alloc(size_t M, size_t N);
/*Создание нулевой матрицы*/
matrix* null_alloc(size_t M, size_t N);
/*Создание единичной матрицы*/
matrix* e_alloc(size_t M, size_t N);
/*Изменение элемента по индексу*/
void change(matrix* a, size_t i, size_t j, element_t x);
/*Получение элемента по индексу*/
element_t poind(matrix* a, size_t i, size_t j);
/*Изменение размера матрицы*/
matrix* re_alloc(matrix* a, size_t m, size_t n);
/*Освобождение памяти и удаление матрицы*/
void matrix_free(matrix* a);
#endif // MATR_H
