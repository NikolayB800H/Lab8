#include "exceptions.h"
#include "matrix.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <istream>
#include <limits>

#include <iostream>

template <class T>
static void swap(T &a, T &b) {
    T tmp = a;
    a = b;
    b = tmp;
}

static void swap(double *&a, double *&b) {
    double *tmp = a;
    a = b;
    b = tmp;
}

static bool almostEqual(double x, double y) {
    return std::fabs(x - y) <=
    #ifdef EPSILON
        EPSILON;
    #else
        std::numeric_limits<double>::epsilon() * std::fabs(x+y)
        || std::fabs(x - y) < std::numeric_limits<double>::min();
    #endif
}

static int minusOneNPower(size_t n) {
    return 1 - 2 * ((int) n % 2);
}

static size_t mtrxSize(size_t rows, size_t cols) {
    return sizeof(Matrix) + sizeof(double*) * rows
           + sizeof(double) * rows * cols;
}

static size_t matricesSizesSum(size_t from, size_t to) {  // Finish rows/cols
    size_t result = 0;                                   // and start rows/cols
    for (; from <= to; ++from) {
        result += mtrxSize(from, from);
    }
    return result;
}

static size_t factPow(size_t cols, size_t rows) {
    size_t result = 1;
    for (; cols > 0; --cols) {
        result *= cols;
    }
    size_t cpy = result;
    for (; rows > 1; --rows) {
        result *= cpy;
    }
    if (rows == 0) {
        result = 1;
    }
    return result;
}

void Matrix::fillSpecial() {
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        for (size_t col_i = 0; col_i < cols; ++col_i) {
            if (row_i == col_i) {
                data[row_i][col_i] = 1;
            } else if (row_i < col_i) {
                data[row_i][col_i] = 1.0 / factPow(col_i, row_i);;
            } else {
                data[row_i][col_i] = (1.0 - 2 * (row_i % 2)) / factPow(col_i, row_i);;
            }
        }
    }
}

Matrix::Matrix(size_t rows, size_t cols) : rows(rows), cols(cols) {
    data = new double *[rows];
    raw_data = new double[rows * cols];
    memset(data, 0, rows);
    memset(raw_data, 0, rows * cols);
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        data[row_i] = raw_data + row_i * cols;
    }
}

Matrix::Matrix(std::istream &is) {
    if (!is.good()) {
        throw InvalidMatrixStream();
    }
    is >> rows >> cols;
    if (is.fail()) {
        throw InvalidMatrixStream();
    }
    if (rows == 0 || cols == 0) {
        throw DimensionMismatch(*this);
    }
    data = new double *[rows];
    raw_data = new double[rows * cols];
    memset(data, 0, rows);
    memset(raw_data, 0, rows * cols);
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        data[row_i] = raw_data + row_i * cols;
    }
    const size_t display_row_cnt = DISPLAY_WIDTH / DOUBLE_WIDTH;
    const size_t parts = (cols - 1) / display_row_cnt + 1;
    for (size_t part_i = 0; part_i < parts; ++part_i) {
        for (size_t row_i = 0; row_i < rows; ++row_i) {
            size_t end = (part_i != parts - 1) ? (part_i + 1) * display_row_cnt : cols;
            for (size_t col_i = part_i * display_row_cnt; col_i < end; ++col_i) {
                is >> data[row_i][col_i];
                if (is.fail()) {
                    throw InvalidMatrixStream();
                }
            }
        }
    }
}

Matrix::Matrix(size_t rows, size_t cols, byte *mem) : rows(rows), cols(cols) {
    if (rows == 0 || cols == 0) {
        throw DimensionMismatch(*this);
    }
    if (mem == nullptr) {
        throw std::bad_alloc();
    }
    data = static_cast<double **>(static_cast<void *>(mem + sizeof(Matrix)));
    raw_data = static_cast<double *>(static_cast<void *>(static_cast<byte *>
                    (static_cast<void *>(data)) + rows * sizeof(double *)));
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        data[row_i] = raw_data + row_i * cols;
    }
}

Matrix::~Matrix() {
    delete[] data;
    delete[] raw_data;
}

Matrix::Matrix(const Matrix &rhs) : Matrix(rhs.rows, rhs.cols) {
    memcpy(raw_data, rhs.raw_data, rows * cols * sizeof(double));
}

Matrix &Matrix::operator=(const Matrix &rhs) {
    if (this == &rhs) {
        return *this;
    }
    delete[] data;
    delete[] raw_data;
    rows = rhs.rows;
    cols = rhs.cols;
    data = new double *[rows];
    raw_data = new double[rows * cols];
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        data[row_i] = raw_data + row_i * cols;
    }
    memcpy(raw_data, rhs.raw_data, rows * cols * sizeof(double));
    return *this;
}

size_t Matrix::getRows() const {
    return rows;
}

size_t Matrix::getCols() const {
    return cols;
}

double Matrix::operator()(size_t i, size_t j) const {
    if (i >= rows || j >= cols) {
        throw OutOfRange(i, j, *this);
    }
    return data[i][j];
}

double &Matrix::operator()(size_t i, size_t j) {
    if (i >= rows || j >= cols) {
        throw OutOfRange(i, j, *this);
    }
    return data[i][j];
}

bool Matrix::operator==(const Matrix &rhs) const {
    if (rhs.rows != rows || rhs.cols != cols) {
        return false;
    }
    bool is_equ = true;
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        for (size_t col_i = 0; col_i < cols; ++col_i) {
            if (!almostEqual(data[row_i][col_i], rhs.data[row_i][col_i])) {
                is_equ = false;
            }
        }
    }
    return is_equ;
}

bool Matrix::operator!=(const Matrix &rhs) const {
    return !(*this == rhs);
}

Matrix Matrix::operator+(const Matrix &rhs) const {
    if (rhs.rows != rows || rhs.cols != cols) {
        throw DimensionMismatch(*this, rhs);
    }
    return sumWithAdditionalSign(*this, rhs, PLUS);
}

Matrix Matrix::operator-(const Matrix &rhs) const {
    if (rhs.rows != rows || rhs.cols != cols) {
        throw DimensionMismatch(*this, rhs);
    }
    return sumWithAdditionalSign(*this, rhs, MINUS);
}

Matrix Matrix::operator*(const Matrix &rhs) const {
    if (cols != rhs.rows) {
        throw DimensionMismatch(*this, rhs);
    }
    Matrix result(rows, rhs.cols);
    for (size_t row_i = 0; row_i < result.rows; ++row_i) {
        for (size_t col_i = 0; col_i < result.cols; ++col_i) {
            double tmp = 0;
            for (size_t vec_i = 0; vec_i < cols; ++vec_i) {
                tmp += data[row_i][vec_i] * rhs.data[vec_i][col_i];
            }
            result.data[row_i][col_i] = tmp;
        }
    }
    return result;
}

Matrix Matrix::operator*(double val) const {
    Matrix result(*this);
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        for (size_t col_i = 0; col_i < cols; ++col_i) {
            result.data[row_i][col_i] *= val;
        }
    }
    return result;
}

void Matrix::deleteZeroStreaks() {
    size_t *use_row = new size_t[rows];
    size_t *use_col = new size_t[cols];
    memset(use_row, 0, rows * sizeof(size_t));
    memset(use_col, 0, cols * sizeof(size_t));
    size_t last_row = 0;
    size_t last_col = 0;
    for (size_t row_i = 0; row_i < cols; ++row_i) {
        for (size_t col_i = 0; col_i < rows; ++col_i) {
            if (!almostEqual(0, data[row_i][col_i])) {
                ++use_row[row_i];
                ++use_col[col_i];
            }
        }
    }
    size_t row_i = 0;
    size_t col_i = 0;
    size_t repeat_cnt = 0;
    size_t prev = 0;
    for (; row_i < rows && repeat_cnt < rows - row_i; ++row_i) {
        if (prev == row_i) {
            ++repeat_cnt;
        } else {
            repeat_cnt = 0;
            prev = row_i;
        }
        if (use_row[row_i] == 0) {
            for (size_t row_j = row_i; row_j < rows - 1; ++row_j) {
                swap(data[row_j], data[row_j + 1]);
                swap(use_row[row_j], use_row[row_j + 1]);
            }
            --row_i;
        }
    }
    for (; col_i < cols && repeat_cnt < cols - col_i; ++col_i) {
        if (prev == col_i) {
            ++repeat_cnt;
        } else {
            repeat_cnt = 0;
            prev = col_i;
        }
        if (use_col[col_i] == 0) {
            for (size_t col_j = col_i; col_j < cols - 1; ++col_j) {
                for (size_t row_i = 0; row_i < rows; ++row_i) {
                    swap(data[row_i][col_j], data[row_i][col_j + 1]);
                }
                swap(use_col[col_j], use_col[col_j + 1]);
            }
            --col_i;
        }
    }
    delete[] use_row;
    delete[] use_col;
    trunc(row_i, col_i);
}

void Matrix::trunc(size_t new_rows, size_t new_cols) {
    double **new_data = new double *[new_rows];
    double *new_raw_data = new double[new_rows * new_cols];
    for (size_t row_i = 0; row_i < new_rows; ++row_i) {
        new_data[row_i] = new_raw_data + row_i * new_cols;
        memcpy(new_data[row_i], data[row_i], new_cols * sizeof(double));
    }
    delete[] data;
    delete[] raw_data;
    rows = new_rows;
    cols = new_cols;
    data = new_data;
    raw_data = new_raw_data;
}

Matrix Matrix::transp() const {
    Matrix result(cols, rows);
    for (size_t row_i = 0; row_i < cols; ++row_i) {
        for (size_t col_i = 0; col_i < rows; ++col_i) {
            result.data[row_i][col_i] = data[col_i][row_i];
        }
    }
    return result;
}

double Matrix::det() const {
    if (rows != cols) {
        throw DimensionMismatch(*this);
    }
    if (rows == 0) {
        return BAD_RESULT;
    }
    if (rows == ONE_ELEM_MTRX_ROWS_CNT) {
        return data[0][0];
    }
    size_t all_size = 0;
    if (rows == MIN_RECUR_SQR_MTRX) {
        all_size = 0;
    } else {
        all_size = matricesSizesSum(MIN_RECUR_SQR_MTRX, rows - 1);
    }
    byte *mem_tmps_raws = new byte[all_size];
    memset(mem_tmps_raws, 0, all_size);
    double result = detWithRaws(mem_tmps_raws);
    delete[] mem_tmps_raws;
    return result;
}

double Matrix::detWithRaws(byte* mem_tmps_raws) const {  // mem_tmps_raws is
    if (rows != cols) {                                 // not const because
        throw DimensionMismatch(*this);                // later there will be
    }                                                 // placement matrices (*)
    if (mem_tmps_raws == nullptr) {
        return BAD_RESULT;
    }
    Matrix *tmp = nullptr;
    double result = 0;
    if (cols > MIN_RECUR_SQR_MTRX) {  // Otherwise don't need less matrix
        tmp = new(mem_tmps_raws) Matrix(rows - 1, cols - 1, mem_tmps_raws);  //
    } else if (cols < MIN_RECUR_SQR_MTRX) {                          // Here(*)
        return BAD_RESULT;
    }
    for (size_t i_of_minor = 0; i_of_minor < cols; ++i_of_minor) {
        double minor_i = 0;
        if (cols == MIN_RECUR_SQR_MTRX) {
            minor_i = data[MIN_RECUR_SQR_MTRX - 1][1 - i_of_minor];  // No less
        } else {                                                     // matrix
            createMinorPart(0, i_of_minor, tmp);
            minor_i = tmp->detWithRaws(mem_tmps_raws             // And here(*)
                                       + mtrxSize(tmp->rows, tmp->rows));
        }
        result += minor_i * data[0][i_of_minor] * minusOneNPower(i_of_minor);
    }
    return result;
}


Matrix Matrix::adj() const {
    if (rows != cols) {
        DimensionMismatch(*this);
    }
    Matrix transposed = transp();
    Matrix result(transposed.cols, transposed.rows);
    Matrix minor(transposed.rows - 1, transposed.cols - 1);
    for (size_t row_i = 0; row_i < result.rows; ++row_i) {
        for (size_t col_i = 0; col_i < result.cols; ++col_i) {
            transposed.createMinorPart(row_i, col_i, &minor);
            result.data[row_i][col_i] = minor.det() * minusOneNPower(row_i + col_i);
        }
    }
    return result;
}

Matrix Matrix::inv() const {
    if (rows != cols) {
        DimensionMismatch(*this);
    }
    double determinant = det();
    if (almostEqual(0, determinant)) {
        throw SingularMatrix();
    }
    return adj() * (1 / det());
}

void Matrix::createMinorPart(size_t row_to_del, size_t col_to_del,
                             Matrix *result) const {
    if ((row_to_del >= rows || col_to_del >= cols) || (rows < 2 || cols < 2)) {
        // throw OutOfRange(row_to_del, col_to_del, *this);
        throw DimensionMismatch(*this);
    }
    if (result->rows != rows - 1 || result->cols != cols - 1) {
        throw DimensionMismatch(*this, *result);
    }
    for (size_t row_i = 0; row_i < result->rows; ++row_i) {
        for (size_t col_i = 0; col_i < result->cols; ++col_i) {
            result->data[row_i][col_i] =
                data[((row_i < row_to_del) ? row_i : (row_i + 1))]
                    [((col_i < col_to_del) ? col_i : (col_i + 1))];
        }
    }
}

Matrix sumWithAdditionalSign(const Matrix &lhs, const Matrix &rhs, Sign sign) {
    if (rhs.rows != lhs.rows || rhs.cols != lhs.cols) {
        throw DimensionMismatch(lhs, rhs);
    }
    Matrix result(lhs);
    for (size_t row_i = 0; row_i < lhs.rows; ++row_i) {
        for (size_t col_i = 0; col_i < lhs.cols; ++col_i) {
            result.data[row_i][col_i] += rhs.data[row_i][col_i] * sign;
        }
    }
    return result;
}

Matrix operator*(double val, const Matrix &matrix) {
    return matrix * val;
}

std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
    os << matrix.rows << ' ' << matrix.cols << '\n';
    const size_t display_row_cnt = DISPLAY_WIDTH / DOUBLE_WIDTH;
    const size_t parts = (matrix.cols - 1) / display_row_cnt + 1;
    for (size_t part_i = 0; part_i < parts; ++part_i) {
        for (size_t row_i = 0; row_i < matrix.rows; ++row_i) {
            size_t end = (part_i != parts - 1) ? (part_i + 1) * display_row_cnt : matrix.cols;
            for (size_t col_i = part_i * display_row_cnt; col_i < end; ++col_i) {
                os << std::setw(DOUBLE_WIDTH) << matrix.data[row_i][col_i];
            }
            os << '\n';
        }
        os << '\n';
    }
    return os;
}
