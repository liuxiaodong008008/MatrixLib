//    MIT License
//
//    Copyright (c) 2019 liu xiaodong
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy
//            of this software and associated documentation files (the "Software"), to deal
//    in the Software without restriction, including without limitation the rights
//            to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//            copies of the Software, and to permit persons to whom the Software is
//    furnished to do so, subject to the following conditions:
//
//    The above copyright notice and this permission notice shall be included in all
//            copies or substantial portions of the Software.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//            AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//    SOFTWARE.


#ifndef __MATRIX_OPERATION_H__
#define __MATRIX_OPERATION_H__
#include "matrix.h"
#include <initializer_list>

double epsilon( double e=-1) {
    static double eps = 1e-12;
    if(!(e<0)) eps = e;
    return eps;
}

template <template <typename> typename M1, \
        typename T1, typename Func,\
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto reduce_vertical(const M1<T1>& mat, const Func& func) {
    using RT = decltype(func(T1(),T1()));
    Matrix<RT> c(1,mat.cols());
    int rn = mat.rows();
    int cn = mat.cols();

    for(int ci=0;ci<cn;ci++) {
        RT rt = mat(0,ci);
        for(int ri=1;ri<rn;ri++) {
            rt=func(rt,mat(ri,ci));
        }
        c(0,ci)=rt;
    }
    return c;
}

template <template <typename> typename M1, \
        typename T1, typename Func,\
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto reduce_horizontal(const M1<T1>& mat, const Func& func) {
    using RT = decltype(func(T1(),T1()));
    Matrix<RT> c(mat.rows(),1);
    int rn = mat.rows();
    int cn = mat.cols();

    for(int ri=0;ri<rn;ri++) {
        RT rt = mat(ri,0);
        for(int ci=1;ci<cn;ci++) {
            rt=func(rt,mat(ri,ci));
        }
        c(ri,0)=rt;
    }
    return c;
}

template <template <typename> typename M1, \
        typename T1, typename Func,\
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto reduce_all(const M1<T1>& mat, const Func& func) {
    using RT = decltype(func(T1(),T1()));
    int n = mat.count();
    RT rt = mat(0);
    for(int i=1;i<n;i++) rt = func(rt,mat(i));

    return Matrix<RT>(1,1,rt);
}

enum class REDUCE_DIRECTION {
    HORIZONTAL,
    VERTICAL,
    ALL
};

template <template <typename> typename M1, \
        typename T1, typename Func,\
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto reduce(const M1<T1>& mat, const Func& func, REDUCE_DIRECTION rd) {
    switch (rd) {
        case REDUCE_DIRECTION::HORIZONTAL:
            return reduce_horizontal(mat,func);
            break;
        case REDUCE_DIRECTION::VERTICAL:
            return reduce_vertical(mat,func);
            break;
        case REDUCE_DIRECTION::ALL:
            return reduce_all(mat,func);
            break;
        default: // note: default is the same as case REDUCE_DIRECTION::ALL
            return reduce_all(mat, func);
            break;
    }
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto sum(const M1<T1>& mat, REDUCE_DIRECTION rd) {
    return reduce(mat,[](const T1&t1,const T1&t2){return t1+t2;},rd);
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto magnitude(const M1<T1>& mat) {
    return sqrt(sum(mat*mat,REDUCE_DIRECTION::ALL)(0));
}


template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto all(const M1<T1>& mat, REDUCE_DIRECTION rd) {
    return reduce(mat,[](const T1&t1,const T1&t2){return t1 && t2;},rd);
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto any(const M1<T1>& mat, REDUCE_DIRECTION rd) {
    return reduce(mat,[](const T1&t1,const T1&t2){return t1 || t2;},rd);
}


template <template <typename> typename M1,template <typename> typename M2, \
        typename T1, typename T2, \
        typename = std::enable_if_t< \
                (!std::is_base_of_v<IsMatrix,   T1>) \
                && (!std::is_base_of_v<IsMatrix,   T2>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>) \
                && (std::is_base_of_v<IMatrix<T2>, M2<T2>>),int>> \
Matrix<decltype(T1()*T2())> matmul(const M1<T1>&m1,const M2<T2>&m2) {
    assert(m1.cols()==m2.rows());
    using RT = decltype(T1()*T2());
    Matrix<RT> ret(m1.rows(),m2.cols());
    int cn = ret.cols();
    int rn = ret.rows();
    int kn = m1.cols();
    for(int r=0;r<rn;r++) {
        for(int c=0;c<cn;c++) {
            RT rt(0);
            for(int k=0;k<kn;k++) rt += m1(r,k)*m2(k,c);
            ret(r,c)=rt;
        }
    }
    return ret;
}

template <typename T1,typename T2, typename  ...Args>
auto matmul(const T1&m1, const T2&m2, Args ...args) {
    return matmul(matmul(m1,m2), args...);
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
void swap_row(const M1<T1>& mat, int r1, int r2) {
    int cn = mat.cols();
    T1 t;
    for(int c=0;c<cn;c++) {
        t = mat(r1,c);
        mat(r1,c) = mat(r2,c);
        mat(r2,c) = t;
    }
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto inverse(const M1<T1>& mat) {
    assert(mat.rows()==mat.cols());
    int n = mat.rows();
    Matrix<T1> b(n,2*n);
    b(Range(Point{0,0},Size{n,n})) = mat;
    b(Range(Point{n,0},Size{n,n})) = Matrix<T1>::eye(n);

    for(int r=0;r<n;r++) {
        T1 max_ = b(r,r);
        int max_r = r;
        for(int i=r;i<n;i++)
            if(max_<b(i,r)) {
                max_=b(i,r);
                max_r = i;
            }
        if(max_r!=r) swap_row(b,max_r,r);
        b.row(r)/=max_;
        for(int i=r+1;i<n;i++) b.row(i)-=b.row(r)*b(i,r);
        for(int i=r-1;i>=0;i--) b.row(i)-=b.row(r)*b(i,r);
    }

    auto re = b(Range(Point{n,0},Size{n,n})).copy();
    return re;
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
bool solve(const M1<T1>& A,const M1<T1>& b, M1<T1>& x) {
    assert(A.rows()==b.rows());
    assert(A.rows()==x.rows());
    assert(b.cols()==1);
    assert(x.cols()==1);

    int n = A.rows();
    int m = A.cols();
    Matrix<T1> Ab(n,m+1);
    Ab(Range(Point{0,0},Size{m,n})) = A;
    Ab(Range(Point{m,0},Size{1,n})) = b;

    for(int r=0;r<n;r++) {
        if(Ab(r,r)<T1(epsilon())) return false;
        Ab.row(r)/=Ab(r,r);
        for(int i=r+1;i<n;i++) b.row(i)-=b.row(r)*Ab(i,r);
        for(int i=r-1;i>=0;i--) b.row(i)-=b.row(r)*Ab(i,r);
    }

    x = Ab(Ab.cols()-1).copy();
    return true;
}


template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto upper_triangular(const M1<T1>& mat) {
    int n = mat.rows();
    Matrix<T1> b(n,n);
    b=mat;

    for(int r=0;r<n;r++) {
        T1 max_ = b(r,r);
        int max_r = r;
        for(int i=r;i<n;i++)
            if(max_<b(i,r)) {
                max_=b(i,r);
                max_r = i;
            }
        if(max_r!=r) swap_row(b,max_r,r);
        b.row(r)/=max_;
        for(int i=r+1;i<n;i++) b.row(i)-=b.row(r)*b(i,r);
    }
    return b;
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto lower_triangular(const M1<T1>& mat) {
    int n = mat.rows();
    Matrix<T1> b(n,n);
    b=mat;

    for(int r=n-1;r>=0;r--) {
        T1 max_ = b(r,r);
        int max_r = r;
        for(int i=r;i>=0;i--)
            if(max_<b(i,r)) {
                max_=b(i,r);
                max_r = i;
            }
        if(max_r!=r) swap_row(b,max_r,r);
        b.row(r)/=max_;
        for(int i=r-1;i>=0;i--) b.row(i)-=b.row(r)*b(i,r);
    }
    return b;
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto lower_triangular_inverse(const M1<T1>& mat) {
    assert(mat.rows()==mat.cols());
    int n = mat.rows();
    Matrix<T1> b(n,2*n);
    b(Range(Point{0,0},Size{n,n})) = mat;
    b(Range(Point{n,0},Size{n,n})) = Matrix<T1>::eye(n);

    for(int r=0;r<n;r++) {
        b.row(r)/=b(r,r);
        for(int i=r+1;i<n;i++) b.row(i)-=b.row(r)*b(i,r);
    }

    auto re = b(Range(Point{n,0},Size{n,n})).copy();
    return re;
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
void LU(const M1<T1>& mat, M1<T1>&L, M1<T1>&U) {
    assert(mat.getSize()==L.getSize());
    assert(mat.getSize()==U.getSize());
    int n = mat.rows();
    Matrix<T1> b(n,2*n);
    b(Range(Point{0,0},Size{n,n})) = mat;
    b(Range(Point{n,0},Size{n,n})) = Matrix<T1>::eye(n);

    for(int r=0;r<n;r++) {
        auto p = b(r,r);
        b.row(r)/=p;
        for(int i=r+1;i<n;i++) b.row(i)-=b.row(r)*b(i,r);
    }

    U = b.col(0,mat.cols());
    L = lower_triangular_inverse(b.col(mat.cols(),mat.cols()));
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
void LDU(const M1<T1>& mat, M1<T1>&L, M1<T1>&D, M1<T1>&U) {
    assert(mat.getSize()==L.getSize());
    assert(mat.getSize()==D.getSize());
    assert(mat.getSize()==U.getSize());
    int n = mat.rows();
    Matrix<T1> b(n,2*n);
    b(Range(Point{0,0},Size{n,n})) = mat;
    b(Range(Point{n,0},Size{n,n})) = D = Matrix<T1>::eye(n);

    for(int r=0;r<n;r++) {
        auto p = b(r,r);
        b.row(r)/=p;
        for(int i=r+1;i<n;i++) b.row(i)-=b.row(r)*b(i,r);
    }

    U = b.col(0,mat.cols());
    L = lower_triangular_inverse(b.col(mat.cols(),mat.cols()));
    for(int r=0;r<n;r++) {
        auto p = L(r,r);
        L.col(r)/=p;
        D(r,r) = p;
    }
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
bool LL(const M1<T1>& mat, M1<T1>& L) {
    assert(mat.rows()==mat.cols());
    assert(L.getSize()==mat.getSize());
    int n = mat.rows();
    auto & b = L;
    b=mat;

    for(int r=0;r<n;r++) {
        auto p = b(r,r);
        for(int i=r+1;i<n;i++) b.row(i)-=b.row(r)*(b(i,r)/p);
        if(p>=0) b.row(r)/=sqrt(p);
        else return false;
    }

    L = b.t();
    return true;
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto GramSchmit(const M1<T1>& mat) {
    Matrix<T1> Q=mat;
    int n = mat.rows();
    int m = mat.cols();

    for(int c=0;c<m;c++) {
        Q.col(c)/=magnitude(Q.col(c));
        auto p = matmul(Q.col(c),Q.col(c).t());
        for(int j=c+1;j<m;j++)
            Q.col(j)-=matmul(p,Q.col(j));
    }
    return Q;
}


template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto GramSchmitQR(const M1<T1>& mat, M1<T1>& Q, M1<T1>& R) {
    assert(mat.rows()==mat.cols());
    assert(Q.getSize()==mat.getSize());
    assert(R.getSize()==mat.getSize());

    Q = GramSchmit(mat);
    R = matmul(Q.t(),mat);
}

//template <typename T>
//auto sign(const T&t) {
//    if(t>=T(0)) return T(1);
//    else return T(-1);
//}
//
//template <template <typename> typename M1, typename T1, \
//        typename = std::enable_if_t< \
//                   (!std::is_base_of_v<IsMatrix,   T1>) \
//                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
//auto HouseholderQR(const M1<T1>& mat, M1<T1>& Q, M1<T1>& R) { //todo: form Q
//    assert(mat.rows()==mat.cols());
//    assert(Q.getSize()==mat.getSize());
//    assert(R.getSize()==mat.getSize());
//
//    int n = mat.cols();
//    int m = mat.rows();
//    R=mat;
//    for(int k=0;k<n;k++) {
//        Matrix<T1> x = R.row(k,m-k).col(k).copy();
//        Matrix<T1> v = x;
//        v(0) += sign(x(0))*magnitude(x);
//        v /= magnitude(v);
//        R.row(k,m-k).col(k,n-k) -= 2*matmul(v,matmul(v.t(),R.row(k,m-k).col(k,n-k)));
//    }
//}



template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto det(const M1<T1>& mat) {
    assert(mat.rows()==mat.cols());
    int n = mat.rows();
    auto b = mat;

    T1 ret(1);

    for(int r=0;r<n;r++) {
        auto p = b(r,r);
        for(int i=r+1;i<n;i++) b.row(i)-=b.row(r)*(b(i,r)/p);
        ret*=p;
    }

    return ret;
}


template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto repeat(const M1<T1>& mat, int rows, int cols) {
    Matrix<T1> R(mat.rows()*rows,mat.cols()*cols,0);
    for(int r=0;r<rows;r++) {
        for(int c=0;c<cols;c++) {
            R.row(r*mat.rows(),mat.rows()).col(c*mat.cols(),mat.cols())
            =mat;
        }
    }
    return R;
}

template <template <typename> typename M1, typename T1, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
auto stack(const std::initializer_list<std::initializer_list<M1<T1>>>& init) {
    int rows = 0;
    int cols = 0;

    for(auto rit=init.begin();rit!=init.end();rit++) {
        rows+=rit->begin()->rows();
    }
    for(auto cit=init.begin()->begin();cit!=init.begin()->end();cit++) {
        cols+=cit->cols();
    }

    Matrix<T1> ret(rows,cols,0);
    int r=0;
    for(auto rit=init.begin();rit!=init.end();rit++) {
        int c=0;
        for(auto cit=rit->begin();cit!=rit->end();cit++) {
            ret(Range({c,r},cit->getSize())) = *cit;
            c+=cit->cols();
        }
        r+=rit->begin()->rows();
    }
    return ret;
}

#endif //__MATRIX_OPERATION_H__
