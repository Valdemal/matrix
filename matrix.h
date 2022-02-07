#ifndef MATRIX_LIB
#define MATRIX_LIB

#include <iostream>
#include <vector>

namespace matrix {

    template<typename baseType>
    class Matrix {
    private:
        std::vector<std::vector<baseType>> _values;
        size_t _m{}, _n{};

        // Инициализация матрицы размера m на n нулями
        inline void init(size_t m, size_t n);

        // Возвращает значение на i-го j-го элемента матрицы
        inline baseType at(size_t i, size_t j) const
        { return _values.at(i).at(j); };

    public:

        // Конструктор по умолчанию
        Matrix() = default;

        // Конструктор матрицы размера m на n
        Matrix(size_t m, size_t n) { init(m, n); };

        // Копирующий конструктор
        Matrix(const Matrix<baseType> &copy);

        // Вовзращает значение ИСТИНА, если матрица квадратная, иначе ЛОЖЬ
        inline bool isSquare() const { return (_m == _n); };

        // Вовзращает количество столбцов матрицы
        inline size_t cols() const { return _n; };

        // Возвращает количество строк матрицы
        inline size_t rows() const { return _m; };

        // Вставляет в матрицу столбец column на позицию j
        void insertColumn(size_t j, const std::vector<baseType> &column);

        // Вставляет в матрицу строку row на позицию i
        void insertRow(size_t i, const std::vector<baseType> &row);

        // Удаляет j-й столбец матрицы
        void eraseColumn(size_t j);

        // Удаляет i-ю строку матрицы
        void eraseRow(size_t i);

        // Возвращает ссылку на i-ю строку матрицы
        std::vector<baseType> &operator[](size_t i) { return _values[i]; }

        // Выводит матрицу a в поток out
        friend std::ostream& operator <<
        (std::ostream &out, const matrix::Matrix<baseType> &a){
            for (size_t i = 0; i < a._m; ++i) {
                for (size_t j = 0; j < a._n; j++)
                    out << a.at(i, j) << " ";

                out << std::endl;
            }
            return out;
        }

        // Считывает матрицу a из потока in
        friend std::istream &operator>>
        (std::istream &in, matrix::Matrix<baseType> &a) {
            for (size_t i = 0; i < a._m; ++i) {
                for (size_t j = 0; j < a._n; j++)
                    in >> a[i][j];
            }

            return in;
        }

        // Возвращает матрицу умноженную на -1
        Matrix<baseType> operator-() const;

        // Возращает сумму матриц
        Matrix<baseType> operator+(const Matrix<baseType> &b) const;

        // Возращает разность матриц
        Matrix<baseType> operator-(const Matrix<baseType> &b) const;

        // Возращает произведение матриц
        Matrix<baseType> operator*(const Matrix<baseType> &b) const;

        // Возращает сумму матрицы и значения val
        Matrix<baseType> operator+(const baseType &val) const;

        // Возращает разность матрицы и значения val
        Matrix<baseType> operator-(const baseType &val) const;

        // Возращает произведение матрицы и значения val
        Matrix<baseType> operator*(const baseType &val) const;

        // Возращает частное матрицы и значения val
        Matrix<baseType> operator/(const baseType &val) const;

        // Возвращает произведение матрицы на вектор-строку v
        Matrix<baseType> mulRowVector(const std::vector<baseType> &v) const;

        // Возвращает произведение матрицы на вектор-столбец v
        Matrix<baseType> mulColVector(const std::vector<baseType> &v) const;

        //Возвращает матрицу возведенную в степень n
        matrix::Matrix<baseType> pow(size_t n) const;

        // Возвращает транспонированную матрицу
        Matrix<baseType> transpose() const;

        // Возвращает определитель матрицы
        baseType det() const;

        // Возвращает обратную матрицу
        Matrix<baseType> inverse() const;

        /* Возвращает решение системы уравнений, заданной матрицей,
         * в виде вектора X. Если решений нет вернет пустой вектор */
        std::vector<baseType> getSolutionSystemLinearEquations();

        /* Выполняет прямой ход метода Гаусса, возвращает количество реальных
         * перестановок строк, произведенных при прямом ходе */
        size_t gaussForward();

        // Возвращает еденичную матрицу размера m на m
        static Matrix<baseType> E(size_t m);

        /* Возвращает значение ИСТИНА, если размеры матриц a и b совпадают,
         * иначе ЛОЖЬ */
        static bool haveEqualSize
        (const Matrix<baseType> &a, const Matrix<baseType> &b);

    private:
        /* Прибавляет к строке modifyStrI строку referenceStrI
         * умноженную на fi */
        void modifyStr(size_t modifyStrI, size_t referenceStrI, baseType fi);

        /* Выполняет обратный ход метода Гаусса для получения решения системы
         * уравнений в матрице, в которой был пременен прямой ход.
         * Возвращает вектор решений, или пустой вектор, если решений нет*/
        std::vector<baseType> gaussReverseForSystem();
    };
}

#include "matrix.tpp"

#endif