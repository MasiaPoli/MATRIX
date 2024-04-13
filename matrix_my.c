#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_my.h"
#include <malloc.h>
struct matrix
{
    size_t m;
    size_t n;
    element_t mat[10][10];
};
matrix* alloc(size_t mm, size_t nn)
{
    if(mm>100)
    {
        return NULL;
    }
    if(nn>100)
    {
        return NULL;
    }
    matrix* ma;
    ma=malloc(sizeof(matrix));
    if(!ma)
        return NULL;
    ma->m=mm;
    ma->n=nn;
    return ma;
}
matrix* null_alloc(size_t m, size_t n)
{
    matrix* ma=alloc(m, n);
    if(!ma)
    {
        return NULL;
    }
    for(size_t i=0; i<m; i++)
    {
        for(size_t j=0; j<n; j++)
        {
            ma->mat[i][j]=0;
        }
    }
    return ma;
}
matrix* e_alloc(size_t m, size_t n)
{
    matrix* ma=alloc(m, n);
    if(!ma)
    {
        return NULL;
    }
    for(size_t i=0; i<m; i++)
    {
        for(size_t j=0; j<n; j++)
        {
            ma->mat[i][j]=0;
            if(i==j)
            {
                ma->mat[i][j]=1;
            }
        }
    }
    return ma;
}
void change(matrix* a, size_t i, size_t j, element_t x)
{
    a->mat[i][j]=x;
}
element_t poind(matrix* a, size_t i, size_t j)
{
    return a->mat[i][j];
}
matrix* sume(matrix* a, matrix* b)
{
    if(!a)
    {
        return NULL;
    }
    if(!b)
    {
        return NULL;
    }
    if(a->m!=b->m)
    {
        return NULL;
    }
    if(a->n!=b->n)
    {
        return NULL;
    }
    for(size_t i=0; i<a->m; i++)
    {
        for(size_t j=0; j<a->n; j++)
        {
            a->mat[i][j]+=b->mat[i][j];
        }
    }
    return a;
}
matrix* re_alloc(matrix* a, size_t M, size_t N)
{
    if((M>100) || (N>100))
    {
        return NULL;
    }
    a->m=M;
    a->n=N;
    return a;
}
void matrix_free(matrix* a)
{
    free(a);
}
matrix* matrix_copy(matrix* a, matrix* b)
{
    if(!a)
    {
        return NULL;
    }
    if(!b)
    {
        b=alloc(a->m, a->n);

    }
    else
    {
        b=re_alloc(b, a->m, a->n);
    }
    if(!b)
    {
        return NULL;
    }
    for(size_t i=0; i<a->m; i++)
    {
        for(size_t j=0; j<a->n; j++)
        {
            b->mat[i][j]=a->mat[i][j];
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
    if(!b)
    {
        b=alloc(a->n, a->m);
    }
    else
    {
        if((b->m!=a->n) ||(b->n!=a->m))
        {
            b=re_alloc(b, a->n, a->m);
        }
    }
    if(!b)
    {
        return NULL;
    }
    for(size_t i=0; i<a->m; i++)
    {
        for(size_t j=0; j<a->n; j++)
        {
            b->mat[j][i]=a->mat[i][j];
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
    matrix* b;
    b=null_alloc(a->n, a->m);
    b=trans(a, b);
    if(!b)
    {
        return 1;
    }
    matrix* c=a;
    a=b;
    matrix_free(c);
    return a;
}
matrix* matrix_umn(matrix* a, matrix* b, matrix* c)
{
    if(a->n!=b->m)
    {
        return NULL;
    }
    if(!c)
    {
        c=alloc(a->m, b->n);
    }
    else
    {
        c=re_alloc(c, a->m, b->n);
    }
    if(!c)
    {
        return NULL;
    }
    for(size_t i=0; i<c->m; i++)
    {
        for(size_t j=0; j<c->n; j++)
        {
            c->mat[i][j]=0;
            for(size_t k=0; k<a->n; k++)
            {
                c->mat[i][j]+=a->mat[i][k]*b->mat[k][j];
            }
        }
    }
    return c;
}
matrix* m_deg(matrix* a, int d)
{
    if(a->m!=a->n)
    {
        return NULL;
    }
    matrix* deg=matrix_copy(a, deg);
    matrix* cur=e_alloc(a->m, a->n);
    for(int i=0; i<31; i++)
    {
        int k=d>>i;
        k&=1;
        matrix* e;
        if(k==1)
        {
            e=matrix_umn(cur, deg, e);
            matrix* g=cur;
            cur=e;
            matrix_free(g);
        }
        e=matrix_umn(deg, deg, e);
        matrix* h=deg;
        deg=e;
        matrix_free(h);
    }
    matrix_free(deg);
    return cur;
}
element_t m_norm(matrix* a)
{
    element_t norm=0;
    for(size_t i=0; i<a->n; i++)
    {
        element_t curnorm=0;
        for(size_t j=0; j<a->m; j++)
        {
            curnorm+=fabs(a->mat[j][i]);
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
    for(size_t i=0; i<a->m; i++)
    {
        for(size_t j=0; j<a->n; j++)
        {
            a->mat[i][j]*=b;
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
    if(a->m!=a->n)
    {
        return NULL;
    }
    element_t A=1.0/0.0;
    if(A==eps)
    {
        return NULL;
    }
    matrix* deg=e_alloc(a->m, a->n);
    matrix* cur=null_alloc(a->m, a->n);
    double N=1;
    while(m_norm(deg)>= eps)
    {
        matrix* x=null_alloc(a->m, a->n);
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
    int b[100];
    for(size_t j=0; j<a->n; j++)
    {
        b[j]=a->mat[x][j];
    }
    for(size_t j=0; j<a->n; j++)
    {
        a->mat[x][j]=a->mat[y][j];
    }
    for(size_t j=0; j<a->n; j++)
    {
        a->mat[y][j]=b[j];
    }
}
void mins(matrix* a, size_t x, size_t y, element_t c)
{
    for(size_t j=0; j<a->n; j++)
    {
        a->mat[x][j]-=a->mat[y][j]*c;
    }
}
void strumn(matrix* a, size_t x, element_t c)
{
    for(size_t j=0; j<a->n; j++)
    {
        a->mat[x][j]*=c;
    }
}
matrix* gauss (matrix* a)
{
    matrix* b=alloc(a->m, 1);
    for(size_t i=0; i<a->m; i++)
    {
        if(a->mat[i][i]==0.0)
        {
            for(size_t j=i; j<a->m; j++)
            {
                if(a->mat[j][i]!=0.0)
                {
                    menmest(a, i, j);
                    break;
                }
            }
        }
        if(a->mat[i][i]==0.0)
        {
            b->mat[0][0]=0.0/0.0;
            return b;
        }
        for(size_t j=i+1; j<a->m; j++)
        {
            mins(a, j, i, a->mat[j][i]/a->mat[i][i]);
        }
    }
    if(a->mat[a->m-1][a->m-1]==0.0)
    {

        b->mat[0][0]=0.0/0.0;
        return b;
    }
    for(size_t i=a->m-1; i>=0; i--)
    {
        element_t x=a->mat[i][a->n-1];
        for(size_t j=a->m-1; j>i; j--)
        {
            x-=b->mat[j][0]*a->mat[i][j];
        }
        b->mat[i][0]=x/a->mat[i][i];
        if(i==0)
        {
            break;
        }
    }
    return b;
}
