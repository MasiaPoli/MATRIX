#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_my.h"
#include <malloc.h>
#include "matr.h"
matrix* sume(matrix* a, matrix* b)
{
    size_t ma=matrix_m(a);
    size_t mb=matrix_m(b);
    size_t na=matrix_n(a);
    size_t nb=matrix_n(b);
    if(!a)
    {
        return NULL;
    }
    if(!b)
    {
        return NULL;
    }
    if(ma!=mb)
    {
        return NULL;
    }
    if(na!=nb)
    {
        return NULL;
    }
    for(size_t i=0; i<ma; i++)
    {
        for(size_t j=0; j<na; j++)
        {
            element_t x=poind(a, i, j);
            element_t y=poind(b, i, j);
            change(a, i, j, x+y);
        }
    }
    return a;
}
matrix* matrix_copy(matrix* a, matrix* b)
{
    if(!a)
    {
        return NULL;
    }
    size_t ma=matrix_m(a);
    size_t na=matrix_n(a);
    if(!b)
    {
        b=alloc(ma, na);

    }
    else
    {
        b=re_alloc(b, ma, na);
    }
    if(!b)
    {
        return NULL;
    }
    for(size_t i=0; i<ma; i++)
    {
        for(size_t j=0; j<na; j++)
        {
            element_t x=poind(a, i, j);
            change(b, i, j, x);
        }
    }
    return b;
}
matrix* sum(matrix* a, matrix* b, matrix* c)
{
    c=matrix_copy(a, c);
    c=sume(c, b);
    if(!c)
    {
        return NULL;
    }
    return c;
}
void e_swap(element_t* a, element_t* b)
{
    element_t c=*a;
    *a=*b;
    *b=c;
}
matrix* trans(matrix* a, matrix* b)
{
    size_t ma=matrix_m(a);
    size_t na=matrix_n(a);
    if(!b)
    {
        b=alloc(na, ma);
    }
    else
    {
        if((matrix_m(b)!=na) ||(matrix_n(b)!=ma))
        {
            b=re_alloc(b, na, ma);
        }
    }
    if(!b)
    {
        return NULL;
    }
    for(size_t i=0; i<ma; i++)
    {
        for(size_t j=0; j<na; j++)
        {
            element_t x=poind(a, i, j);
            change(b, j, i, x);
        }
    }
    return b;
}
matrix* self_t(matrix* a)
{
    if(!a)
    {
        return NULL;
    }
    size_t ma=matrix_m(a);
    size_t na=matrix_n(a);
    matrix* b=null_alloc(na, ma);
    b=trans(a, b);
    if(!b)
    {
        return NULL;
    }
    matrix* c=a;
    a=b;
    matrix_free(c);
    return a;
}
matrix* matrix_umn(matrix* a, matrix* b, matrix* c)
{
    size_t ma=matrix_m(a);
    size_t mb=matrix_m(b);
    size_t na=matrix_n(a);
    size_t nb=matrix_n(b);
    if(na!=mb)
    {
        return NULL;
    }
    if(!c)
    {
        c=alloc(ma, nb);
    }
    else
    {
        c=re_alloc(c, ma, nb);
    }
    if(!c)
    {
        return NULL;
    }
    for(size_t i=0; i<ma; i++)
    {
        for(size_t j=0; j<nb; j++)
        {
            change(c, i, j,0);
            for(size_t k=0; k<na; k++)
            {
                element_t xc=poind(c, i, j);
                element_t xa=poind(a, i, k);
                element_t xb=poind(b, k, j);
                change(c, i, j, xa*xb+xc);
            }
        }
    }
    return c;
}
element_t m_norm(matrix* a)
{
    size_t ma=matrix_m(a);
    size_t na=matrix_n(a);
    element_t norm=0;
    for(size_t i=0; i<na; i++)
    {
        element_t curnorm=0;
        for(size_t j=0; j<ma; j++)
        {
            curnorm+=fabs(poind(a, j, i));
        }
        if(curnorm>norm)
        {
            norm=curnorm;
        }
    }
    return norm;
}
matrix* onnumber(matrix* a, element_t b)
{
    if(!a)
    {
        return NULL;
    }
    size_t ma=matrix_m(a);
    size_t na=matrix_n(a);
    for(size_t i=0; i<ma; i++)
    {
        for(size_t j=0; j<na; j++)
        {
            element_t x=poind(a, i, j);
            change(a, i, j, x*b);
        }
    }
    return a;
}
matrix* m_exp(matrix* a, element_t eps)
{
    if(eps!=eps)
    {
        return NULL;
    }
    size_t ma=matrix_m(a);
    size_t na=matrix_n(a);
    if(ma!=na)
    {
        return NULL;
    }
    element_t A=1.0/0.0;
    if(A==eps)
    {
        return NULL;
    }
    matrix* deg=e_alloc(ma, na);
    matrix* cur=null_alloc(ma, na);
    double N=1;
    while(m_norm(deg)>= eps)
    {
        matrix* x=null_alloc(ma, na);
        cur=sume(cur, deg);
        x=matrix_umn(deg, a, x);
        x=onnumber(x, 1/N);
        N++;
        matrix* y=deg;
        deg=x;
        matrix_free(y);
    }
    cur=sume(cur, deg);
    matrix_free(deg);
    return cur;
}
void menmest(matrix* a, size_t x, size_t y)
{
    size_t na=matrix_n(a);
    element_t* b=malloc(na*sizeof(element_t));
    for(size_t j=0; j<na; j++)
    {
        b[j]=poind(a, x, j);
    }
    for(size_t j=0; j<na; j++)
    {
        element_t w=poind(a, y, j);
        change(a, x, j, w);
    }
    for(size_t j=0; j<na; j++)
    {
        change(a, y, j, b[j]);
    }
    free(b);
}
void mins(matrix* a, size_t x, size_t y, element_t c)
{
    size_t na=matrix_n(a);
    for(size_t j=0; j<na; j++)
    {
        element_t k=poind(a, x, j);
        element_t l=poind(a, y, j);
        change(a, x, j, k-c*l);
    }
}
void strumn(matrix* a, size_t x, element_t c)
{
     size_t na=matrix_n(a);
    for(size_t j=0; j<na; j++)
    {
        element_t k=poind(a, x, j);
        change(a, x, j, k*c);
    }
}
matrix* gauss (matrix* a)
{
    size_t ma=matrix_m(a);
    size_t na=matrix_n(a);
    matrix* b=alloc(ma, 1);
    for(size_t i=0; i<ma; i++)
    {
        if(poind(a, i, i)==0.0)
        {
            for(size_t j=i; j<ma; j++)
            {
                if(poind(a, j, i)!=0.0)
                {
                    menmest(a, i, j);
                    break;
                }
            }
        }
        if(poind(a, i, i)==0.0)
        {
            change(b, 0, 0, 0.0/0.0);
            return b;
        }
        element_t k=poind(a, i, i);
        for(size_t j=i+1; j<ma; j++)
        {
            element_t l=poind(a, j, i);
            mins(a, j, i, l/k);
        }
    }
    if(poind(a, ma-1, ma-1)==0.0)
    {

        change(b, 0, 0, 0.0/0.0);
        return b;
    }
    for(size_t j=ma; j>0; j--)
    {
        size_t i=j-1;
        element_t x=poind(a, i, na-1);
        for(size_t j=ma-1; j>i; j--)
        {
            x-=poind(b, j, 0)*poind(a, i, j);
        }
        element_t u=poind(a, i, i);
        change(b, i, 0, x/u);
    }
    return b;
}
