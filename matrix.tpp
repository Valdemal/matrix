#include <numeric>
#include "matrix.h"

/* Конструкторы */

template<typename baseType>
inline void matrix::Matrix<baseType>::init(size_t m, size_t n) {
    _m = m;
    _n = n;

    _values = std::vector<std::vector<baseType>>(m);
    for (auto &v : _values)
        v = std::vector<baseType>(n);
}

template<typename baseType>
matrix::Matrix<baseType>::Matrix(const matrix::Matrix<baseType> &copy) {
    init(copy._m, copy._n);

    for (size_t i = 0; i < _m; ++i)
        for (size_t j = 0; j < _n; ++j)
            _values[i][j] = copy._values[i][j];

}


/* Арифметические операции */

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::operator-() const {
    Matrix<baseType> result(_m, _n);
    for (size_t i = 0; i < _m; ++i)
        for (size_t j = 0; j < _n; ++j)
            result[i][j] = -_values[i][j];

    return result;
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::operator+(const matrix::Matrix<baseType> &b) const {
    if (haveEqualSize(*this, b)) {

        Matrix<baseType> result(_m, _n);

        for (size_t i = 0; i < _m; ++i)
            for (size_t j = 0; j < _n; ++j)
                result[i][j] = this->at(i, j) + b.at(i, j);

        return result;
    } else
        throw std::runtime_error("Can't sum non-equal size matrices!");
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::operator-(const matrix::Matrix<baseType> &b) const {
    if (haveEqualSize(*this, b)) {
        Matrix<baseType> result(_m, _n);

        for (size_t i = 0; i < _m; ++i)
            for (size_t j = 0; j < _n; ++j)
                result[i][j] = this->at(i, j) - b.at(i, j);

        return result;
    } else
        throw std::runtime_error("Can't difference non-equal size matrices!");
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::operator*(const matrix::Matrix<baseType> &b) const {
    if (this->_n == b._m) {
        Matrix<baseType> result(this->_m, b._n);

        for (size_t i = 0; i < this->_m; i++)
            for (size_t j = 0; j < b._n; ++j)
                for (size_t k = 0; k < this->_n; ++k)
                    result[i][j] += this->at(i, k) * b.at(k, j);

        return result;
    } else
        throw std::runtime_error("Matrix sizes are not suitable for multiplication!");
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::operator+(const baseType &val) const {
    Matrix<baseType> result(*this);

    for (size_t i = 0; i < _m; ++i)
        for (size_t j = 0; j < _n; ++j)
            result[i][j] += val;

    return result;
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::operator-(const baseType &val) const {
    Matrix<baseType> result(*this);

    for (size_t i = 0; i < _m; ++i)
        for (size_t j = 0; j < _n; ++j)
            result[i][j] -= val;

    return result;
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::operator*(const baseType &val) const {
    Matrix<baseType> result(*this);

    for (size_t i = 0; i < _m; ++i)
        for (size_t j = 0; j < _n; ++j)
            result[i][j] *= val;

    return result;
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::operator/(const baseType &val) const {
    Matrix<baseType> result(*this);

    for (size_t i = 0; i < _m; ++i)
        for (size_t j = 0; j < _n; ++j)
            result[i][j] /= val;

    return result;
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::mulRowVector(const std::vector<baseType> &v) const {
    Matrix row(1, v.size());

    for (size_t j = 0; j < row._n; ++j)
        row[0][j] = v[j];

    return row * *this;
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::mulColVector(const std::vector<baseType> &v) const {
    Matrix column(v.size(), 1);

    for (size_t i = 0; i < column._m; ++i)
        column[i][0] = v[i];

    return *this * column;
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::pow(size_t n) const {
    if (n == 0)
        return E(_n);
    else {
        matrix::Matrix<baseType> res = *this;

        for (size_t i = 1; i < n; ++i)
            res = res * *this;

        return res;
    }
}


/* Операторы модифицирующие размеры матрицы*/

template<typename baseType>
void matrix::Matrix<baseType>::insertColumn(size_t j, const std::vector<baseType> &column) {
    for (size_t i = 0; i < _m; ++i) {
        auto posIter = _values[i].begin() + j;
        _values[i].insert(posIter, column[i]);
    }
    _n++;
}

template<typename baseType>
void matrix::Matrix<baseType>::insertRow(size_t i, const std::vector<baseType> &row) {
    _values.insert(_values.begin() + i, row);
    _m++;
}

template<typename baseType>
void matrix::Matrix<baseType>::eraseColumn(size_t j) {
    for (size_t i = 0; i < _m; ++i) {
        auto posIter = _values[i].begin() + j;
        _values[i].erase(posIter);
    }
    _n--;
}

template<typename baseType>
void matrix::Matrix<baseType>::eraseRow(size_t i) {
    _values.erase(_values.begin() + i);

    _m--;
}


/* Характерные опрераторы матриц
 * (определитель, транспонированная матрица, обратная матрица)
 * */

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::transpose() const {
    Matrix<baseType> result(_n, _m);

    for (size_t i = 0; i < _n; ++i)
        for (size_t j = 0; j < _m; ++j)
            result[i][j] = at(j, i);

    return result;
}

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::inverse() const {
    // добавить проверку на квадратность матрицы и ненулевость определителя
    if (isSquare()) {
        Matrix temp(*this);

        auto e = Matrix::E(temp._m);

        // Присоединение еденичной матрицы
        for (size_t i = 0; i < e._m; ++i)
            temp.insertColumn(temp._n, e[i]);

        temp.gaussForward();

        Matrix res(this->_m, this->_n);

        for (size_t j = 0; j < _n; j++) {
            std::vector<double> t = temp.gaussReverseForSystem();

            if (t.empty())
                throw std::runtime_error("Can't find inverse to degenerate matrix");

            for (size_t i = 0; i < res._m; ++i)
                res[t.size() - i - 1][t.size() - j - 1] = t[i];

            temp.eraseColumn(temp.cols() - 1);
        }

        return res;
    } else
        throw std::runtime_error("Can't find inverse matrix to non-square matrix!");
}

template<typename baseType>
baseType matrix::Matrix<baseType>::det() const {
    if (this->isSquare()) {
        Matrix temp(*this);

        // Получение вектора сумм элементов строк
        std::vector<baseType> sums;
        for (size_t i = 0; i < _m; ++i)
            sums.push_back(
                    std::accumulate(temp[i].begin(), temp[i].end(), 0));

        // Добавление вектора к матрице в качестве столбца
        temp.insertColumn(temp.cols(), sums);

        size_t k = temp.gaussForward(); // Выполнение прямого хода Гаусса

        baseType diagonalMultiply = 1;
        for (size_t i = 0; i < _m; ++i)
            diagonalMultiply *= temp[i][i];

        return diagonalMultiply * (k % 2 == 0 ? 1 : -1);
    } else
        throw std::runtime_error("Can't find determinant of non-square matrix!");
}


/* Метод Гаусса */

template<typename baseType>
size_t matrix::Matrix<baseType>::gaussForward() {
    size_t permutationsCnt = 0;

    for (size_t k = 0; k < _m; ++k) {
        /* поиск индекса максимального по модулю элемента столбца,
         * начиная с k-й строки */
        int maxIndex = k;
        for (size_t i = k + 1; i < _m; ++i)
            if (abs(_values[i][k]) > abs(_values[maxIndex][k]))
                maxIndex = i;

        /*Выполнение перестановки,
         * если k - не строка с максимальным по модулю элементом  */
        if (maxIndex != k) {
            _values[k].swap(_values[maxIndex]);
            permutationsCnt++;
        }

        for (size_t i = k + 1; i < _m; ++i) {
            baseType fi = -(_values[i][k] / _values[k][k]);
            modifyStr(i, k, fi);
        }
    }

    return permutationsCnt;
}

template<typename baseType>
void matrix::Matrix<baseType>::modifyStr(size_t modifyStrI, size_t referenceStrI, baseType fi) {
    for (size_t j = 0; j < _n; ++j)
        _values[modifyStrI][j] += _values[referenceStrI][j] * fi;
}

template<typename baseType>
std::vector<baseType>
matrix::Matrix<baseType>::getSolutionSystemLinearEquations() {
    Matrix temp(*this);

    temp.gaussForward();

    return temp.gaussReverseForSystem();
}

template<typename baseType>
std::vector<baseType> matrix::Matrix<baseType>::gaussReverseForSystem() {
    std::vector<baseType> xVector;

    // проверка на отсутствие нулей на главной диагонали
    for (int i = 0; i < _m; ++i)
        if ((*this)[i][i] == 0)
            return xVector;

    for (int k = this->_m - 1; k >= 0; k--) {
        baseType p = (*this)[k][this->_n - 1];

        baseType sum = 0;

        for (size_t j = 0; j < xVector.size(); j++)
            sum += (*this)[k][this->_m - j - 1] * xVector.at(j);

        baseType x = (p - sum) / (*this)[k][k];

        xVector.push_back(x);
    }

    return xVector;
}


/* Статические методы */

template<typename baseType>
matrix::Matrix<baseType> matrix::Matrix<baseType>::E(size_t m) {
    Matrix<baseType> res(m, m);

    for (size_t i = 0; i < m; ++i)
        res[i][i] = 1;

    return res;
}

template<typename baseType>
bool matrix::Matrix<baseType>::haveEqualSize(const matrix::Matrix<baseType> &a, const matrix::Matrix<baseType> &b) {
    return a._m == b._m && a._n == b._n;
}