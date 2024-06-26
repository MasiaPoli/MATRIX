#include "matr.h"
#include <malloc.h>
struct matrix
{
    size_t MAKS;
    size_t m;
    size_t n;
    element_t* mat;
};
size_t matrix_m(matrix* a)
{
    return a->m;
}
size_t matrix_n(matrix* a)
{
    return a->n;
}
matrix* alloc(size_t M, size_t N)
{
    matrix* a= malloc(sizeof(matrix));
    size_t x=M*N;
    if(x<10)
    {
        x=10;
    }
    a->MAKS=x;
    a->mat=malloc(x*sizeof(element_t));
    a->m=M;
    a->n=N;
    return a;
}
matrix* null_alloc(size_t M, size_t N)
{
    matrix* a= alloc(M, N);
    for(size_t i=0; i<M*N; i++)
    {
        a->mat[i]=0;
    }
    return a;
}
matrix* e_alloc(size_t M, size_t N)
{
    matrix* a= alloc(M, N);
    for(size_t i=0; i<M*N; i++)
    {
        size_t x=i/N;
        size_t y=i-x*N;
        if(x!=y)
        {
            a->mat[i]=0;
        }
        else
        {
            a->mat[i]=1;
        }
    }
    return a;
}
void change(matrix* a, size_t i, size_t j, element_t x)
{
    size_t N=a->n;
    a->mat[i*N+j]=x;
}
element_t poind(matrix* a, size_t i, size_t j)
{
    size_t N=a->n;
    return a->mat[i*N+j];
}
matrix* re_alloc(matrix* a, size_t M, size_t N)
{
    if(M*N>a->MAKS)
    {
        size_t x=M*N;
        if(x<2*a->MAKS)
        {
            x=2*a->MAKS;
        }
        a->mat=realloc(a->mat, x*sizeof(element_t));
         a->MAKS=x;
    }
    a->m=M;
    a->n=N;
    return a;
}
void matrix_free(matrix* a)
{
    free(a->mat);
    free(a);
}
