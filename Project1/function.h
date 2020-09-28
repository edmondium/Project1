#pragma once
#ifndef FUNCTION_
#define FUNCTION_
#ifdef T_DEBUG
#define T_LOG(x) x
#else
#define T_LOG(x)
#endif // T_DEBUG
#define HPoffset 8
#include <memory.h>
#include <cassert>
using namespace std;
template<class T> class function;
template<class T> class function1D;
template<class T> class FunProxy;
template<class T> class function2D;

template<class T> class function
{
protected:
    T* f;
    int N0;
    int N;
    function() : f(NULL), N0(0), N(0) {};
    explicit function(int N_) : N0(N_), N(N_) {};
    ~function() {};
    function(const function&) {};
    template<class U> friend class function2D;

public:
    T& operator[](int i)
    {
        if (i<N)
        {
            assert("Out of range in function[]");
            return f[i]; //trying to put here first
        }
    }
    const T& operator[](int i) const
    {
        if (i<N)
        {
            assert("Out of range in function[]");
            return f[i]; //trying to put here first
        }
    }
    const T& last() const
    { 
        return f[N - 1];
    }
    int size() const
    { 
        return N;
    }
    int FullSize() const
    { 
        return N0;
    }
    T* MemPt()
    { 
        return f;
    }
    const T* MemPt() const
    { 
        return f;
    }
    function& operator+=(const function& m);
    function& operator*=(const T& m);
    function& operator=(const T& c);
};

template<class T>
inline function<T>& function<T>::operator+=(const function& m)
{
    T_LOG(if (N != m.size())
    {
        cerr << "Functions not of equal length! Can't sum!" << endl;
    })
    for (int i = 0; i < N; i++)
    {
        f[i] += m[i];
    }
    return *this;
}

template<class T>
inline function<T>& function<T>::operator*=(const T& m)
{
    for (int i = 0; i < N; i++)
    {
        f[i] *= m;
    }
    return *this;
}

template <class T>
inline function<T>& function<T>::operator=(const T& c)
{
    T_LOG(if (N <= 0)
    {
        cerr << "Size of function is non positive! " << N << endl;
    })
    for (int i = 0; i < N; i++)
    {
        f[i] = c;
    }
    return *this;
}
//**************************************************************************************************//
// One dimensional functions derived from function<T>. It has its own constructors and destructors. //
//**************************************************************************************************//
template <class T> class function1D : public function<T>
{
public:
    function1D() {};
    void resize(int N_);
    explicit function1D(int N_);
    ~function1D();
    function1D(const function1D& m);
    function1D& operator=(const function1D& m);
    function1D& operator=(const T& c)
    { 
        function<T>::operator=(c);
        return *this;
    }
};

template<class T>
inline void function1D<T>::resize(int n)
{
    if (n > this->N0)
    {
        if (this->f)
        {
            delete[] this->f;
        }
        this->f = new T[n];
        this->N0 = n;
    }
    this->N = n;
}

template<class T>
inline function1D<T>::function1D(int N_) : function<T>(N_)
{
    this->f = new T[N_];
}

template<class T>
inline function1D<T>::~function1D()
{
    delete[] this->f;
    this->f = NULL;
}

template<class T>
inline function1D<T>::function1D(const function1D& m)
{
    resize(m.N);
    copy(m.f, m.f + this->N, this->f);
}

template <class T>
inline function1D<T>& function1D<T>::operator=(const function1D<T>& m)
{
    resize(m.N);
    copy(m.f, m.f + this->N, this->f);
    return *this;
}
//*******************************//
// A proxy class for function2D. //
//*******************************//
template <class T> class FunProxy : public function<T> {
public:
    void initialize(int N_, T* f_);
    void reinitialize(int N_, T* f_);
    void resize(int N_);
    FunProxy& operator=(const function<T>& m);
    ~FunProxy() {};
};

template <class T>
inline void FunProxy<T>::initialize(int N_, T* f_)
{
    this->N = this->N0 = N_;
    this->f = f_;
}

template <class T>
inline void FunProxy<T>::reinitialize(int N_, T* f_)
{
    this->N = N_;
    this->f = f_;
}

template <class T>
inline void FunProxy<T>::resize(int N_)
{
    if (N_ > this->N0)
    {
        cerr << "Can't resize FunProxy, to small FunProxy!" << endl;
    }
    else
    {
        this->N = N_;
    }
}

template <class T>
inline FunProxy<T>& FunProxy<T>::operator=(const function<T>& m)
{
    resize(m.size());
    copy(m.MemPt(), m.MemPt() + this->N, this->f);
    return *this;
}
//*******************************************************//
// Two dimentional function<T> derived from function<T>. //
//*******************************************************//
template<class T> class function2D
{
protected:
    void* memory;
    T* data;
    FunProxy<T>* f;
    int N0;
    int Nd0;
    int N;
    int Nd;
public:
    function2D() : memory(NULL), N0(0), Nd0(0), N(0), Nd(0) {};
    function2D(int N_, int Nd_);
    ~function2D();
    void resize(int N_, int Nd_);
    function2D& operator=(const function2D& m);
    function2D& operator=(const T& u);
    function2D& operator*=(const T& x);
    function2D& operator+=(double x);
    function2D& operator+=(const function2D& m);
    function2D& operator-=(double x);
    function2D& operator-=(const function2D& m);
    FunProxy<T>& operator[](int i)
    {
        if (i<N)
        {
            assert("Out of range in function2D[]");
            return f[i]; //trying to put here first
        }
    }
    const FunProxy<T>& operator[](int i) const
    {
        if (i<N)
        {
            assert("Out of range in function2D[]");
            return f[i]; //trying to put here first
        }
    }
    const T& operator()(int i, int j) const
    {
        if (i<N && j<Nd)
        {
            assert("Out of range in function2D(i,j)");
            return f[i].f[j]; //trying to put here first
        }
    }
    T& operator()(int i, int j)
    {
        if (i < N && j < Nd)
        {
            assert("Out of range in function2D(i,j)");
            return f[i].f[j]; //trying to put here first
        }
    }
    T* MemPt()
    { 
        return data;
    }
    const T* MemPt() const
    { 
        return data;
    }
    const int SizeN() const
    { 
        return N;
    }
    const int SizeNd() const
    { 
        return Nd;
    }
    const int FullSizeN() const
    { 
        return N0;
    }
    const int FullSizeNd() const
    { 
        return Nd0;
    }
    const int gw() const
    { 
        return Nd0;
    }
    void product(const std::string& transa, const std::string& transb, const function2D& A, const function2D& B, const T& alpha, const T& beta);
};

template<class T>
function2D<T>::function2D(int N_, int Nd_) : N0(N_), Nd0(Nd_), N(N_), Nd(Nd_)
{
    memory = operator new (sizeof(FunProxy<T>) * N0 + sizeof(T) * Nd0 * N0 + HPoffset);
    if (memory != NULL)
    {
        assert("Out of memory");
    }
    f = new (memory) FunProxy<T>[N0];
    int offset = sizeof(FunProxy<T>) * N0 + HPoffset;
    data = reinterpret_cast<T*>(static_cast<char*>(memory) + offset);
    for (int i = 0; i < N0; i++)
    {
        f[i].initialize(Nd0, data + i * Nd0);
    }
}

template<class T>
function2D<T>::~function2D()
{
    for (int i = 0; i < N0; i++)
    {
        f[i].~FunProxy<T>();
    }
    operator delete(memory);
    memory = NULL;
}

template <class T>
inline void function2D<T>::resize(int N_, int Nd_)
{
    if (N_ > N0 || Nd_ > Nd0)
    {
        int MemorySize = sizeof(FunProxy<T>) * N_ + sizeof(T) * Nd_ * N_ + HPoffset;
        operator delete(memory);
        memory = operator new (MemorySize);
        if (memory != NULL)
        {
            assert("Out of memory");
        }
        N = N0 = N_;
        Nd = Nd0 = Nd_;
        f = new (memory) FunProxy<T>[N];
        int offset = sizeof(FunProxy<T>) * N + HPoffset;
        data = reinterpret_cast<T*>(static_cast<char*>(memory) + offset);
        for (int i = 0; i < N; i++)
        {
            f[i].initialize(Nd, data + i * Nd);
        }
    }
    else
    {
        N = N_;
        Nd = Nd_;
    }
}

template <class T>
inline function2D<T>& function2D<T>::operator=(const function2D& m)
{
    if (m.N <= N0 && m.Nd <= Nd0)
    {
        N = m.N;
        Nd = m.Nd;
        for (int i = 0; i < N; i++)
        {
            memcpy(f[i].f, m.f[i].f, sizeof(T) * Nd);
        }
    }
    else
    {
        int MemorySize = sizeof(FunProxy<T>) * m.N + sizeof(T) * m.Nd * m.N + HPoffset;
        operator delete(memory);
        memory = operator new (MemorySize);
        if (memory != NULL)
        {
            assert("Out of memory");
        }
        memcpy(memory, m.memory, MemorySize);
        N = N0 = m.N;
        Nd = Nd0 = m.Nd;
        f = new (memory) FunProxy<T>[N];
        int offset = sizeof(FunProxy<T>) * N + HPoffset;
        data = reinterpret_cast<T*>(static_cast<char*>(memory) + offset);
        for (int i = 0; i < N; i++)
        {
            f[i].initialize(Nd, data + i * Nd);
        }
    }
    return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator=(const T& u)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < Nd; j++)
        {
            f[i].f[j] = u;
        }
    }
    return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator*=(const T& x)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < Nd; j++)
        {
            f[i][j] *= x;
        }
    }
    return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator+=(double x)
{
    if (N != Nd || !Nd || !N)
    {
        cerr << "Can't add number to non-square matrix!" << endl;
        return *this;
    }
    for (int i = 0; i < Nd; i++)
    {
        f[i][i] += x;
    }
    return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator+=(const function2D& m)
{
    if (N != m.N || Nd != m.Nd)
    {
        cerr << "Can't sum different matrices!" << endl;
        return *this;
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < Nd; j++)
        {
            f[i][j] += m[i][j];
        }
    }
    return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator-=(double x)
{
    if (N != Nd || !N || !Nd)
    {
        cerr << "Can't add number to non-square matrix!" << endl;
        return *this;
    }
    for (int i = 0; i < Nd; i++)
    {
        f[i][i] -= x;
    }
    return *this;
}

template <class T>
inline function2D<T>& function2D<T>::operator-=(const function2D& m)
{
    if (N != m.N || Nd != m.Nd)
    {
        cerr << "Can't sum different matrices!" << endl;
        return *this;
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < Nd; j++)
        {
            f[i][j] -= m[i][j];
        }
    }
    return *this;
}
// Test for non-qadratic matrices
void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, const double* alpha, const double* A, const int* lda, const double* B, const int* ldb, const double* beta, double* C, const int* ldc)
{

}

inline void xgemm(const char* transa, const char* transb, const int m, const int n, const int k, const double alpha, const double* A, const int lda, const double* B, const int ldb, const double beta, double* C, const int ldc)
{
    dgemm_(transa, transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

template <class T>
inline void function2D<T>::product(const string& transa, const string& transb, const function2D& A, const function2D& B, const T& alpha = 1, const T& beta = 0)
{
    if (transa != "N" && transa != "T")
    { 
        cerr << "Don't recognize your task. Specify how to multiply matrices in dgemm!" << endl;
        return;
    }
    if (transa == "N" && transb == "N")
    {
        if (A.Nd != B.N || !B.Nd || !A.N || !A.Nd || !B.N || N0 < A.N || Nd0 < B.Nd)
        {
            cerr << " Matrix sizes not correct" << endl;
        }
        xgemm("N", "N", B.Nd, A.N, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
        N = A.N;
        Nd = B.Nd;
    }
    else if (transa == "T" && transb == "N")
    {
        if (A.N != B.N || !B.Nd || !A.Nd || !A.N || !B.N || N0 < A.Nd || Nd0 < B.Nd)
        {
            cerr << " Matrix sizes not correct" << endl;
        }
        xgemm("N", "T", B.Nd, A.Nd, B.N, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
        N = A.Nd;
        Nd = B.Nd;
    }
    else if (transa == "N" && transb == "T")
    {
        if (A.Nd != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0 < A.N || Nd0 < B.N)
        {
            cerr << " Matrix sizes not correct" << endl;
        }
        xgemm("T", "N", B.N, A.N, B.Nd, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
        N = A.N;
        Nd = B.N;
    }
    else if (transa == "T" && transb == "T")
    {
        if (A.N != B.Nd || !B.N || !A.N || !A.Nd || !B.Nd || N0 < A.Nd || Nd0 < B.N)
        {
            cerr << " Matrix sizes not correct" << endl;
        }
        xgemm("T", "T", B.N, A.Nd, B.Nd, alpha, B.MemPt(), B.Nd0, A.MemPt(), A.Nd0, beta, MemPt(), Nd0);
        N = A.Nd;
        Nd = B.N;
    }
}
#endif // FUNCTION_