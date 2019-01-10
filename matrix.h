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


#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <algorithm>
#include <cassert>
#include <ostream>
#include <string>
#include <typeinfo>
#include <type_traits>

struct Size;
struct Range;

struct HasSize {
    virtual Size getSize() const = 0;
    virtual ~HasSize() = default;
};

struct Point {
    int x,y;

    Point(int x, int y) : x(x), y(y) {}
    Point swap() const { return {y,x}; }

    Range operator,(const Size& size) const;
    Range operator,(const Point& other) const;

    Size operator-(const Point& other) const;

    operator Size() const;
    Point operator+(const Size& size) const;
    Point operator-(const Size& size) const;
    Point operator+() const { return {x,y}; };
    Point operator-() const { return {-x,-y}; };

    bool operator==(const Point &rhs) const {
        return x == rhs.x &&
               y == rhs.y;
    }

    bool operator!=(const Point &rhs) const {
        return !(rhs == *this);
    }

    friend std::ostream &operator<<(std::ostream &os, const Point &point) {
        os << "x: " << point.x << " y: " << point.y;
        return os;
    }
};

struct Size {
    int w,h;

    Size(int w, int h) : w(w), h(h) {}

    operator Point() const {
        return {w,h};
    }
    int count() const { return w*h; }
    int index(const Point& pt) const { return w*pt.y+pt.x; }
    int index(int row, int col) const { return w*row+col; }

    Size swap() const { return {h,w}; }

    bool operator==(const Size &rhs) const {
        return w == rhs.w &&
               h == rhs.h;
    }

    bool operator!=(const Size &rhs) const {
        return !(rhs == *this);
    }

    friend std::ostream &operator<<(std::ostream &os, const Size &size) {
        os << "w: " << size.w << " h: " << size.h;
        return os;
    }
};

struct Range : HasSize {
    Point beg;
    Size size;

    Range(const Point &beg, const Size &size) : beg(beg), size(size) {}
    Range(const Size &size) : beg(0,0), size(size) {}
    Range(int x,int y,int w,int h) : beg(x,y), size(w,h) {}

    Point beg_point() const {return beg;}
    Point end_point() const {return beg+size;}
    int left()   const {return beg.x;}
    int right()  const {return beg.x+size.w-1;}
    int top()    const {return beg.y;}
    int bottom() const {return beg.y+size.h-1;}

    Range background(const Range& bk) const {
        return Range{{bk.beg.x+this->beg.x,bk.beg.y+this->beg.y},
                     this->size};
    }

    int count()  const {return size.count();}

    bool operator==(const Range &rhs) const {
        return beg == rhs.beg &&
               size == rhs.size;
    }

    bool operator!=(const Range &rhs) const {
        return !(rhs == *this);
    }

    friend std::ostream &operator<<(std::ostream &os, const Range &range) {
        os << "beg: " << range.beg << " size: " << range.size;
        return os;
    }

    virtual Size getSize() const {
        return size;
    }
};

Range Point::operator,(const Size &size) const {
    return {*this,size};
}

Point::operator Size() const {
    return {x,y};
}

Point Point::operator+(const Size &size) const {
    return {x+size.w,y+size.h};
}

Point Point::operator-(const Size &size) const {
    return {x-size.w,y-size.h};
}

Range Point::operator,(const Point &other) const {
    return {*this,Point(other)-(*this)};
}

Size Point::operator-(const Point &other) const {
    return {this->x-other.x,this->y-other.y};
}

template <typename T>
struct Matrix;

template <typename T>
struct MatrixView;

template <typename T>
struct MatrixInputer;

template <typename T>
struct MatrixIterator;

template <typename T>
struct MatrixVectorWiseView;


struct IsMatrix {
    virtual ~IsMatrix() = default;
};

template <typename T>
struct IMatrix : HasSize,IsMatrix  {

    virtual ~IMatrix() = default;

    virtual Matrix<T> copy() const;

    virtual int count() const {return getSize().count();}
    virtual std::string classname() const { return "IMatrix"; }

    virtual T& operator()(const Point& idx2d) const {
        return (*this)(this->getSize().index(idx2d));
    }
    virtual T& operator()(int row, int col) const {
        return (*this)(this->getSize().index(row,col));
    }
    virtual T& operator()(int idx) const = 0;
    virtual Size getSize() const = 0 ;

    virtual T& at(const Point& idx2d) const { return (*this)(idx2d); }
    virtual T& at(int row, int col) const { return (*this)(row,col); }
    virtual T& at(int idx) const { return (*this)(idx); }

    virtual int rows() const { return getSize().h; }
    virtual int cols() const { return getSize().w; }

    virtual MatrixView<T> operator()(const Range& rng) const = 0;
    virtual MatrixView<T> view(Range rng) const { return (*this)(rng); }
    virtual MatrixView<T> view() const { return (*this)(Range{Point{0,0},getSize()}); }

    virtual MatrixView<T> col(int c) const { return (*this)(Range(c,0,1,this->rows())); }
    virtual MatrixView<T> row(int r) const { return (*this)(Range(0,r,this->cols(),1)); }

    virtual MatrixView<T> col(int from, int cnt) const { return (*this)(Range(from,0,cnt,this->rows())); }
    virtual MatrixView<T> row(int from, int cnt) const { return (*this)(Range(0,from,this->cols(),cnt)); }

    virtual MatrixVectorWiseView<T> colwise() const;
    virtual MatrixVectorWiseView<T> rowwise() const;

    virtual MatrixIterator<T> coliter() const;
    virtual MatrixIterator<T> rowiter() const;

    template <typename A>
    IMatrix<T>& operator=(const IMatrix<A>& _mat) const {
        assert(this->getSize()==_mat.getSize());
        int cnt = this->count();
        for(int i=0;i<cnt;i++) {
            (*this)(i)=_mat(i);
        }
        return (*this);
    }

    virtual void setValue(const T& t) {
        int n = this->count();
        for(int i=0;i<n;i++) (*this)(i)=t;
    }

    template <typename A=T> A atc(int row, int col) const {
        return (A)((*this)(row,col));
    }
    template <typename A=T> A atc(Point idx2d) const {
        return (A)((*this)(idx2d));
    }
    template <typename A=T> A atc(int idx) const {
        return (A)((*this)(idx));
    }
    template <typename A> Matrix<A> cast() const;

    virtual Matrix<T> t() const {
        Matrix<T> r(this->getSize().swap());
        int rn = this->rows();
        int cn = this->cols();
        for(int y=0;y<rn;y++)
            for(int x=0;x<cn;x++)
                r(x,y) = (*this)(y,x);
        return r;
    }

    friend std::ostream &operator<<(std::ostream &os, const IMatrix<T> &mat) {
        Size sz = mat.getSize();
        os << mat.classname() << ": " << sz.h << "x" << sz.w << "\n";
        for(int row = 0; row<sz.h; row++) {
            os << "        ";
            for(int col = 0; col<sz.w; col++)
                os << mat(row,col) << (col!=sz.w-1 ? ",\t" : "\n");
        }
        return os;
    }

    virtual MatrixInputer<T> operator<<(T t) const;


    static Matrix<T> ones(int _rows, int _cols) {
        return Matrix<T>(_rows, _cols, T(1));
    }

    static Matrix<T> zeros(int _rows, int _cols) {
        return Matrix<T>(_rows, _cols);
    }

    static Matrix<T> eye(int _rows, int _cols) {
        auto mat = Matrix<T>::zeros(_rows, _cols);
        auto m = std::min(_rows,_cols);
        for(int i=0;i<m;i++) mat(i,i)=T(1);
        return mat;
    }

    static Matrix<T> eye(int _n) {
        return Matrix<T>::eye(_n, _n);
    }


    static Matrix<T> ones_like(const HasSize& a) {
        return Matrix<T>(a.getSize(), T(1));
    }

    static Matrix<T> zeros_like(const HasSize& a) {
        return Matrix<T>(a.getSize());
    }

    static Matrix<T> eye_like(const HasSize& a) {
        auto mat = Matrix<T>::zeros(a.getSize());
        auto m = std::min(mat.getSize().w,mat.getSize().h);
        for(int i=0;i<m;i++) mat(i,i)=T(1);
        return mat;
    }
};

template <typename T>
struct MatrixInputer{
    int input_idx;
    IMatrix<T>& mat;
    MatrixInputer<T>&operator,(T t) {
        assert(input_idx<mat.count());
        mat(input_idx)=t;
        input_idx++;
        return *this;
    };
    friend struct IMatrix<T>;
private:
    MatrixInputer(int input_idx, const IMatrix<T> &mat)
            : input_idx(input_idx), mat(const_cast<IMatrix<T> &>(mat)) {}
    MatrixInputer(const IMatrix<T> &mat, T first_value)
            : input_idx(1), mat(const_cast<IMatrix<T> &>(mat)) {
        mat(0)=first_value;
    }
};

template<typename T>
MatrixInputer<T> IMatrix<T>::operator<<(T t) const {
    return MatrixInputer<T>(const_cast<IMatrix<T> &>(*this),t);
}

template <typename T>
struct MatrixIterator {
    Matrix<T> & _mat;
    Range _cur_rng;
    Range _iter_rng;
    Size _offset;

    MatrixIterator(const IMatrix<T> &_mat, const Range &_cur_rng, const Size &_offset)
            : _mat(const_cast<IMatrix<T> &>(_mat).view()._mat),
              _cur_rng(_cur_rng.background(_mat.view()._rng)),
              _iter_rng(_mat.view()._rng), _offset(_offset) {}

    virtual ~MatrixIterator() = default;

    virtual MatrixView<T> get() const;

    virtual void forward() {
        this->_cur_rng.beg = this->_cur_rng.beg+this->offset();
    }
    virtual bool valid() const {
        return _cur_rng.beg.x>=0
               && _cur_rng.beg.y>=0
               && _cur_rng.beg.x+_cur_rng.size.w<=_iter_rng.size.w
               && _cur_rng.beg.y+_cur_rng.size.h<=_iter_rng.size.h;
    }
    virtual Size offset() const { return _offset; }

    friend std::ostream &operator<<(std::ostream &os, const MatrixIterator &iterator) {
        os << "_mat: " << iterator._mat << " _cur_rng: " << iterator._cur_rng << " _iter_rng: " << iterator._iter_rng
           << " _offset: " << iterator._offset;
        return os;
    }
};


template <typename T>
struct MatrixView: IMatrix<T> {
    Matrix<T> & _mat;
    Range _rng;

    using IMatrix<T>::operator<<;
    using IMatrix<T>::operator();
    //using IMatrix<T>::operator=;

    virtual ~MatrixView() = default;

    virtual std::string classname() const { return "MatrixView"; }
    virtual Matrix<T> & refmat() const { return this->_mat; };
    virtual T& operator()(int row, int col) const;
    virtual T& operator()(int idx) const;
    virtual MatrixView<T> operator()(const Range& rng) const {
        return  MatrixView<T>(this->_mat,Range{this->_rng.beg+rng.beg,rng.size});
    }

    virtual Size getSize() const { return this->_rng.size; }

    friend struct IMatrix<T>;
    friend struct Matrix<T>;


    template <typename A>
    MatrixView<T>& operator=(const IMatrix<A>& mv){
        assert(this->getSize()==mv.getSize());
        int cnt = this->count();
        for(int i=0;i<cnt;i++) {
            (*this)(i)=mv(i);
        }
        return (*this);
    }

    MatrixView<T>& operator=(const MatrixView<T>& mv){
        assert(this->getSize()==mv.getSize());
        int cnt = this->count();
        for(int i=0;i<cnt;i++) {
            (*this)(i)=mv(i);
        }
        return (*this);
    }

    MatrixView(const MatrixView<T>& mv):_mat(mv._mat),_rng(mv._rng) {
    }
    MatrixView(const Matrix<T> &_mat, const Range &_rng)
            : _mat(const_cast<Matrix<T>&>(_mat)), _rng(_rng) {}
};


template<typename T>
MatrixView<T> MatrixIterator<T>::get() const {
    return _mat.view(_cur_rng);
}

template <typename T>
struct MatrixVectorWiseView: HasSize {
    Matrix<T> & _mat;
    Range _rng;
    Size _offset;
    Range _init_iter_rng;

    virtual ~MatrixVectorWiseView() = default;

    virtual MatrixIterator<T> iterator() const {
        return MatrixIterator<T>(_mat,_init_iter_rng,_offset);
    }

    friend struct IMatrix<T>;
    friend struct Matrix<T>;

    virtual Size getSize() const {
        return _rng.getSize();
    }

private:
    MatrixVectorWiseView(Matrix<T> &_mat, const Range &_rng,
                         const Size &_offset, const Range &_init_iter_rng) :
            _mat(_mat), _rng(_rng), _offset(_offset), _init_iter_rng(_init_iter_rng) {}
};


template <typename T>
struct Matrix:IMatrix<T>{
    Size _size;
    T * _data;

    using IMatrix<T>::operator<<;
    //using IMatrix<T>::operator=;
    using IMatrix<T>::operator();

    T& operator()(int idx) const {
        return _data[idx];
    }

    Matrix()
            : _size(1,1),
              _data(new T[_size.count()]{T{0}}) {}

    Matrix(Size size)
            : _size(size),
              _data(new T[size.count()]{T{0}}) {}

    Matrix(int _rows, int _cols)
            : _size(_cols,_rows),
              _data(new T[_rows*_cols]{T{0}}) {}

    Matrix(int _rows, int _cols, const T& value)
            : _size(_cols,_rows),
              _data(new T[_rows*_cols]{T{0}}) {
        std::fill_n(_data,_rows*_cols,value);
    }

    Matrix(const Size& size, const T& value)
            : _size(size),
              _data(new T[size.count()]{T(0)}) {
        std::fill_n(_data,this->count(),value);
    }

    Matrix(Matrix<T>&& _mat)
            : _size(_mat.getSize()),
              _data(_mat._data) {
        _mat._data = nullptr;
    }

    Matrix(const Matrix<T>& _mat)
            : _size(_mat.getSize()),
              _data(new T[_mat.count()]{T{0}}) {
        std::copy(_mat._data,_mat._data+this->count(),_data);
    }

    Matrix(const IMatrix<T>& _mat)
            : _size(_mat.getSize()),
              _data(new T[_mat.count()]{T(0)}) {
        int cnt = _mat.count();
        for(int i=0; i<cnt; i++)
            (*this)(i)=_mat(i);
    }

    virtual std::string classname() const { return "Matrix"; }


    T * data() const {
        return _data;
    }

    Size getSize() const override {
        return _size;
    }

    MatrixView<T> operator()(const Range& rng) const override {
        return MatrixView<T>(*this, rng);
    }

    Matrix<T>& operator=(Matrix<T>&& mat) {
        if (_data!= nullptr) {
            delete[] _data;
            _data = nullptr;
        }
        this->_size = mat._size;
        this->_data = mat._data;
        mat._data = nullptr;
        return *this;
    }

    template <typename A>
    Matrix<T>& operator=(const IMatrix<A>& _mat){
        if (this->getSize() != _mat.getSize()) {
            this->_size = _mat.getSize();
            if (_data!= nullptr) {
                delete[] _data;
                _data = nullptr;
            }
            _data = new T[_mat.count()]{T(0)};
        }

        int cnt = _mat.count();
        for(int i=0; i<cnt; i++)
            (*this)(i)=_mat(i);

        return *this;
    }

    Matrix<T>& operator=(const Matrix<T>& _mat){
        if (this->getSize() != _mat.getSize()) {
            this->_size = _mat.getSize();
            if (_data!= nullptr) {
                delete[] _data;
                _data = nullptr;
            }
            _data = new T[_mat.count()]{T(0)};
        }

        int cnt = _mat.count();
        for(int i=0; i<cnt; i++)
            (*this)(i)=_mat(i);

        return *this;
    }


    Matrix<T>& operator=(const T& a){
        if(_data!= nullptr) {
            int cnt = this->count();
            for(int i=0;i<cnt;i++) (*this)(i)=(a);
        }
        return *this;
    }


    virtual ~Matrix() {
        if (_data!= nullptr) {
            delete[] _data;
            _data = nullptr;
        }
    }
};


template<typename T>
template<typename A>
Matrix<A> IMatrix<T>::cast() const {
    Matrix<A> mat(this->getSize());
    int cnt = this->count();
    for(int i=0;i<cnt;i++) {
        mat(i)=(*this)(i);
    }
    return mat;
}

template<typename T>
Matrix<T> IMatrix<T>::copy() const {
    return Matrix<T>(*this);
}

template<typename T>
MatrixVectorWiseView<T> IMatrix<T>::colwise() const {
    auto v  = this->view();
    return MatrixVectorWiseView<T>(v._mat,v._rng,{1,0},{{0,0},{1,v.rows()}});
}

template<typename T>
MatrixVectorWiseView<T> IMatrix<T>::rowwise() const {
    auto v  = this->view();
    return MatrixVectorWiseView<T>(v._mat,v._rng,{0,1},{{0,0},{v.cols(),1}});
}

template<typename T>
MatrixIterator<T> IMatrix<T>::coliter() const {
    return MatrixIterator<T>(this->view(), this->col(0)._rng, Size(1,0));
}

template<typename T>
MatrixIterator<T> IMatrix<T>::rowiter() const {
    return MatrixIterator<T>(this->view(), this->row(0)._rng, Size(0,1));
}

template<typename T>
T &MatrixView<T>::operator()(int idx) const {
    return (*this)(idx/this->_rng.size.w, idx%this->_rng.size.w);
}

template<typename T>
T &MatrixView<T>::operator()(int row, int col) const {
    return this->_mat(row+this->_rng.beg.y, col+this->_rng.beg.x);
}

template<typename T1, typename T2, typename T3>
inline void mat_apply1(const IMatrix<T1>& a, IMatrix<T2>& dst,
                       T3 func) {
    assert(a.getSize()==dst.getSize());

    int cnt = a.count();
    for(int i=0;i<cnt;i++) dst(i) = func(a(i));
}


template<typename T1, typename T2, typename T3, typename T4>
inline void mat_apply2(const IMatrix<T1>& a, const IMatrix<T2>& b, IMatrix<T3>& dst,
                       T4 func) {
    assert(a.getSize()==b.getSize());
    assert(a.getSize()==dst.getSize());

    int cnt = a.count();
    for(int i=0;i<cnt;i++) dst(i) = func(a(i),b(i));
}


template<typename T1, typename T2, typename T3>
inline void mat_apply_inplace(IMatrix<T1>& a, const IMatrix<T2>& b,
                              T3 func) {
    assert(a.getSize()==b.getSize());
    int cnt = a.count();
    for(int i=0;i<cnt;i++) func(a(i),b(i));
}

template<typename T1, typename T2, typename T3, typename T4>
inline void mat_apply_with_right_number(const IMatrix<T1>& a, T2 b,
                                        IMatrix<T3>& dst, T4 func) {
    assert(a.getSize()==dst.getSize());

    int cnt = a.count();
    for(int i=0;i<cnt;i++)  dst(i) = func(a(i),b);
}

template<typename T1, typename T2, typename T3>
inline void mat_apply_with_right_number_inplace(IMatrix<T1>& a, T2 b, T3 func) {
    int cnt = a.count();
    for(int i=0;i<cnt;i++)  func(a(i),b);
}

template<typename T1, typename T2, typename T3, typename T4>
inline void mat_apply_with_left_number(T1 a, const IMatrix<T2>& b,
                                       IMatrix<T3>& dst, T4 func) {
    assert(b.getSize()==dst.getSize());
    int cnt = b.count();
    for(int i=0;i<cnt;i++)  dst(i) = func(a, b(i));
}

template<typename T1, typename T2, typename T3, typename T4>
inline void mat_apply_left_vectorwise(const MatrixVectorWiseView<T1>& a, const IMatrix<T2>& b, Matrix<T3>& c, T4 func) {
    assert(a.getSize()==c.getSize());
    auto ita = a.iterator();
    auto itc = MatrixIterator<T3>(c,{{0,0},ita._cur_rng.size},a._offset);
    assert(ita.valid() && ita._cur_rng.size == b.getSize());
    for(;ita.valid()&&itc.valid();ita.forward(),itc.forward()) {
        itc.get() = func(ita.get(),b);
    }
}

template<typename T1, typename T2, typename T3, typename T4>
inline void mat_apply_right_vectorwise(const IMatrix<T2>& a, const MatrixVectorWiseView<T1>& b, Matrix<T3>& c, T4 func) {
    assert(b.getSize()==c.getSize());
    auto itb = b.iterator();
    auto itc = MatrixIterator<T3>(c,{{0,0},itb._cur_rng.size},a._offset);
    assert(itb.valid() && itb._cur_rng.size == a.getSize());
    for(;itb.valid()&&itc.valid();itb.forward(),itc.forward()) {
        itc.get() = func(a,itb.get());
    }
}

#define SINGLE_MATRIX_OPERATOR(OP) \
template <template <typename> typename M1, \
        typename T1, \
        typename = std::enable_if_t< \
                (!std::is_base_of_v<IsMatrix,   T1>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
Matrix<decltype(OP T1())> operator OP(const M1<T1>& a) { \
    Matrix<decltype(OP T1())> c(a.getSize()); \
    mat_apply1(a,c,[](const T1&t1){ return OP t1; }); \
    return c; \
}


#define MATRIX_MATRIX_OPERATOR(OP) \
template <template <typename> typename M1,template <typename> typename M2, \
        typename T1, typename T2, \
        typename = std::enable_if_t< \
                (!std::is_base_of_v<IsMatrix,   T1>) \
                && (!std::is_base_of_v<IsMatrix,   T2>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>) \
                && (std::is_base_of_v<IMatrix<T2>, M2<T2>>),int>> \
Matrix<decltype(T1() OP T2())> operator OP(const M1<T1>& a, const M2<T2>& b) { \
    Matrix<decltype(T1() OP T2())> c(a.getSize()); \
    mat_apply2(a,b,c,[](const T1&t1,const T2&t2){ return t1 OP t2; }); \
    return c; \
}

#define MATRIX_MATRIX_ASSIGN_OPERATOR(OP) \
template <template <typename> typename M1,template <typename> typename M2, \
        typename T1, typename T2, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (!std::is_base_of_v<IsMatrix,   T2>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>) \
                && (std::is_base_of_v<IMatrix<T2>, M2<T2>>),int>> \
M1<T1>& operator OP(M1<T1>& a, const M2<T2>& b) { \
    mat_apply_inplace(a,b,[](T1&t1,const T2&t2){t1 OP t2;}); \
    return a; \
}

#define MATRIX_ELEMENT_OPERATOR(OP) \
template <template <typename> typename M, typename T1, typename T2, \
          typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (!std::is_base_of_v<IsMatrix,   T2>) \
                && (std::is_base_of_v<IMatrix<T1>, M<T1>>),int>> \
Matrix<decltype(T1() OP T2())> operator OP(const M<T1>& a, const T2& b) { \
    Matrix<decltype(T1() OP T2())> c(a.getSize()); \
    mat_apply_with_right_number(a,b,c,[](const T1&t1, const T2&t2){return t1 OP t2;});\
    return c; \
}

#define ELEMENT_MATRIX_OPERATOR(OP) \
template <template <typename> typename M, typename T1, typename T2, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (!std::is_base_of_v<IsMatrix,   T2>) \
                && (std::is_base_of_v<IMatrix<T2>, M<T2>>),int>> \
Matrix<decltype(T1() OP T2())> operator OP(const T1& a, const M<T2>& b) { \
    Matrix<decltype(T1() OP T2())> c(b.getSize()); \
    mat_apply_with_left_number(a,b,c,[](const T1&t1, const T2&t2){return t1 OP t2;});\
    return c; \
}

#define MATRIX_ELEMENT_ASSIGN_OPERATOR(OP) \
template <template <typename> typename M, typename T1, typename T2, \
          typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (!std::is_base_of_v<IsMatrix,   T2>) \
                && (std::is_base_of_v<IMatrix<T1>, M<T1>>),int>> \
M<T1>& operator OP(M<T1>& a, const T2& b) { \
    mat_apply_with_right_number_inplace(a,b,[](T1&t1, const T2&t2){t1 OP t2;});\
    return a; \
}

#define VECTORWISE_MATRIX_OPERATOR(OP) \
template <template <typename> typename M2, \
        typename T1, typename T2, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (!std::is_base_of_v<IsMatrix,   T2>) \
                && (std::is_base_of_v<IMatrix<T2>, M2<T2>>),int>> \
Matrix<decltype(T1() OP T2())> operator OP(const MatrixVectorWiseView<T1>& a, const M2<T2>& b) { \
    Matrix<decltype(T1() OP T2())> c(a.getSize());\
    mat_apply_left_vectorwise(a,b,c,[](const auto& t1, const auto& t2){ return t1 OP t2;});\
    return c; \
}

#define MATRIX_VECTORWISE_OPERATOR(OP) \
template <template <typename> typename M1, \
        typename T1, typename T2, \
        typename = std::enable_if_t< \
                   (!std::is_base_of_v<IsMatrix,   T1>) \
                && (!std::is_base_of_v<IsMatrix,   T2>) \
                && (std::is_base_of_v<IMatrix<T1>, M1<T1>>),int>> \
Matrix<decltype(T1() OP T2())> operator OP(const M1<T1>& a, const MatrixVectorWiseView<T2>& b) { \
    Matrix<decltype(T1() OP T2())> c(b.getSize());\
    mat_apply_right_vectorwise(a,b,c,[](const auto& t1, const auto& t2){ return t1 OP t2;});\
    return c; \
}

SINGLE_MATRIX_OPERATOR(+)
SINGLE_MATRIX_OPERATOR(-)
SINGLE_MATRIX_OPERATOR(~)
SINGLE_MATRIX_OPERATOR(!)

MATRIX_MATRIX_OPERATOR(+ )
MATRIX_MATRIX_OPERATOR(- )
MATRIX_MATRIX_OPERATOR(* )
MATRIX_MATRIX_OPERATOR(/ )
MATRIX_MATRIX_OPERATOR(% )
MATRIX_MATRIX_OPERATOR(& )
MATRIX_MATRIX_OPERATOR(| )
MATRIX_MATRIX_OPERATOR(^ )
MATRIX_MATRIX_OPERATOR(> )
MATRIX_MATRIX_OPERATOR(>=)
MATRIX_MATRIX_OPERATOR(< )
MATRIX_MATRIX_OPERATOR(<=)
MATRIX_MATRIX_OPERATOR(==)
MATRIX_MATRIX_OPERATOR(!=)
MATRIX_MATRIX_OPERATOR(&&)
MATRIX_MATRIX_OPERATOR(||)
MATRIX_MATRIX_OPERATOR(<<)
MATRIX_MATRIX_OPERATOR(>>)

MATRIX_MATRIX_ASSIGN_OPERATOR(+=)
MATRIX_MATRIX_ASSIGN_OPERATOR(-=)
MATRIX_MATRIX_ASSIGN_OPERATOR(*=)
MATRIX_MATRIX_ASSIGN_OPERATOR(/=)
MATRIX_MATRIX_ASSIGN_OPERATOR(%=)
MATRIX_MATRIX_ASSIGN_OPERATOR(&=)
MATRIX_MATRIX_ASSIGN_OPERATOR(|=)
MATRIX_MATRIX_ASSIGN_OPERATOR(^=)
MATRIX_MATRIX_ASSIGN_OPERATOR(<<=)
MATRIX_MATRIX_ASSIGN_OPERATOR(>>=)

MATRIX_ELEMENT_OPERATOR(+ )
MATRIX_ELEMENT_OPERATOR(- )
MATRIX_ELEMENT_OPERATOR(* )
MATRIX_ELEMENT_OPERATOR(/ )
MATRIX_ELEMENT_OPERATOR(% )
MATRIX_ELEMENT_OPERATOR(& )
MATRIX_ELEMENT_OPERATOR(| )
MATRIX_ELEMENT_OPERATOR(^ )
MATRIX_ELEMENT_OPERATOR(> )
MATRIX_ELEMENT_OPERATOR(>=)
MATRIX_ELEMENT_OPERATOR(< )
MATRIX_ELEMENT_OPERATOR(<=)
MATRIX_ELEMENT_OPERATOR(==)
MATRIX_ELEMENT_OPERATOR(!=)
MATRIX_ELEMENT_OPERATOR(&&)
MATRIX_ELEMENT_OPERATOR(||)
//MATRIX_ELEMENT_OPERATOR(<<)
//MATRIX_ELEMENT_OPERATOR(>>)

ELEMENT_MATRIX_OPERATOR(+ )
ELEMENT_MATRIX_OPERATOR(- )
ELEMENT_MATRIX_OPERATOR(* )
ELEMENT_MATRIX_OPERATOR(/ )
ELEMENT_MATRIX_OPERATOR(% )
ELEMENT_MATRIX_OPERATOR(& )
ELEMENT_MATRIX_OPERATOR(| )
ELEMENT_MATRIX_OPERATOR(^ )
ELEMENT_MATRIX_OPERATOR(> )
ELEMENT_MATRIX_OPERATOR(>=)
ELEMENT_MATRIX_OPERATOR(< )
ELEMENT_MATRIX_OPERATOR(<=)
ELEMENT_MATRIX_OPERATOR(==)
ELEMENT_MATRIX_OPERATOR(!=)
ELEMENT_MATRIX_OPERATOR(&&)
ELEMENT_MATRIX_OPERATOR(||)
//ELEMENT_MATRIX_OPERATOR(<<)
//ELEMENT_MATRIX_OPERATOR(>>)

MATRIX_ELEMENT_ASSIGN_OPERATOR(+=)
MATRIX_ELEMENT_ASSIGN_OPERATOR(-=)
MATRIX_ELEMENT_ASSIGN_OPERATOR(*=)
MATRIX_ELEMENT_ASSIGN_OPERATOR(/=)
MATRIX_ELEMENT_ASSIGN_OPERATOR(%=)
MATRIX_ELEMENT_ASSIGN_OPERATOR(&=)
MATRIX_ELEMENT_ASSIGN_OPERATOR(|=)
MATRIX_ELEMENT_ASSIGN_OPERATOR(^=)
//MATRIX_ELEMENT_ASSIGN_OPERATOR(<<=)
//MATRIX_ELEMENT_ASSIGN_OPERATOR(>>=)

VECTORWISE_MATRIX_OPERATOR(+ )
VECTORWISE_MATRIX_OPERATOR(- )
VECTORWISE_MATRIX_OPERATOR(* )
VECTORWISE_MATRIX_OPERATOR(/ )
VECTORWISE_MATRIX_OPERATOR(% )
VECTORWISE_MATRIX_OPERATOR(& )
VECTORWISE_MATRIX_OPERATOR(| )
VECTORWISE_MATRIX_OPERATOR(^ )
VECTORWISE_MATRIX_OPERATOR(> )
VECTORWISE_MATRIX_OPERATOR(>=)
VECTORWISE_MATRIX_OPERATOR(< )
VECTORWISE_MATRIX_OPERATOR(<=)
VECTORWISE_MATRIX_OPERATOR(==)
VECTORWISE_MATRIX_OPERATOR(!=)
VECTORWISE_MATRIX_OPERATOR(&&)
VECTORWISE_MATRIX_OPERATOR(||)
VECTORWISE_MATRIX_OPERATOR(<<)
VECTORWISE_MATRIX_OPERATOR(>>)

MATRIX_VECTORWISE_OPERATOR(+ )
MATRIX_VECTORWISE_OPERATOR(- )
MATRIX_VECTORWISE_OPERATOR(* )
MATRIX_VECTORWISE_OPERATOR(/ )
MATRIX_VECTORWISE_OPERATOR(% )
MATRIX_VECTORWISE_OPERATOR(& )
MATRIX_VECTORWISE_OPERATOR(| )
MATRIX_VECTORWISE_OPERATOR(^ )
MATRIX_VECTORWISE_OPERATOR(> )
MATRIX_VECTORWISE_OPERATOR(>=)
MATRIX_VECTORWISE_OPERATOR(< )
MATRIX_VECTORWISE_OPERATOR(<=)
MATRIX_VECTORWISE_OPERATOR(==)
MATRIX_VECTORWISE_OPERATOR(!=)
MATRIX_VECTORWISE_OPERATOR(&&)
MATRIX_VECTORWISE_OPERATOR(||)
MATRIX_VECTORWISE_OPERATOR(<<)
MATRIX_VECTORWISE_OPERATOR(>>)


#endif //__MATRIX_H__
