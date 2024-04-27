#ifndef MATRIX_MY_H
#define MATRIX_MY_H
typedef struct matrix matrix;
typedef double element_t;
struct matrix;
matrix* sume(matrix* a, matrix* b);
matrix* sum(matrix* a, matrix* b, matrix* c);
matrix* matrix_copy(matrix* a, matrix* b);
void e_swap(element_t* a, element_t* b);
matrix* trans(matrix* a, matrix* b);
matrix* self_t(matrix* a);
matrix* matrix_umn(matrix* a, matrix* b, matrix* c);
element_t m_norm(matrix* a);
matrix* onnumber(matrix* a, element_t b);
matrix* m_exp(matrix* a, element_t eps);
void menmest(matrix* a, size_t x, size_t y);
void mins(matrix* a, size_t x, size_t y, element_t c);
void strumn(matrix* a, size_t x, element_t c);
matrix* gauss (matrix* a);
#endif // MATRIX_MY_H
