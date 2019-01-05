//
// Created by liuxi on 2018/12/30.
//

#ifndef __MATRIX_OPERATION_H__
#define __MATRIX_OPERATION_H__
#include "matrix.h"

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

#endif //__MATRIX_OPERATION_H__
