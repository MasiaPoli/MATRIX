#ifndef MATRIX_MY_H
#define MATRIX_MY_H
typedef struct matrix matrix;
typedef double element_t;
struct matrix;
matrix* alloc(size_t m, size_t n);
matrix* null_alloc(size_t m, size_t n);
matrix* e_alloc(size_t m, size_t n);
void change(matrix* a, size_t i, size_t j, element_t x);
element_t poind(matrix* a, size_t i, size_t j);
matrix* sume(matrix* a, matrix* b);
matrix* re_alloc(matrix* a, size_t m, size_t n);
void matrix_free(matrix* a);
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
