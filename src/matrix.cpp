#include "exceptions.h"
#include "matrix.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <istream>
#include <limits>

#define LAB10_MSGS

#ifdef LAB10_MSGS
#include <iostream>
#endif

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

static double factPow(size_t cols, size_t rows) {
    double result = 1;
    for (; cols > 0; --cols) {
        result *= cols;
    }
    if (result == 1 || rows == 0) {
        return 0;
    }
    double cpy = result;
    result = 1 / result;
    for (; rows > 1; --rows) {
        result /= cpy;
    }
    if (rows == 0) {
        result = 1;
    }
    return result;
}

void Matrix::setPrecision(size_t new_precision) {
    if (new_precision > DISPLAY_WIDTH - ADDITIONAL_WIDTH) {
        throw BadMatrixPrecision();
    }
    precision = new_precision;
}

void Matrix::setStyle(DoubleStyle new_style) {
    if (new_style != SCIENTIFIC && new_style != FIXED) {
        throw BadMatrixStyle();
    }
    style = new_style;
}

void Matrix::fillSpecial() {
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        for (size_t col_i = 0; col_i < cols; ++col_i) {
            if (row_i == col_i) {
                data[row_i][col_i] = 1;
            } else if (row_i < col_i) {
                data[row_i][col_i] = factPow(col_i, row_i);
            } else {
                data[row_i][col_i] = (1.0 - 2 * (row_i % 2)) * factPow(col_i, row_i);
            }
        }
    }
}

void Matrix::fillE() {
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        for (size_t col_i = 0; col_i < cols; ++col_i) {
            data[row_i][col_i] = static_cast<double>(static_cast<int>(row_i == col_i));
        }
    }
}

Matrix::Matrix(size_t rows, size_t cols) : style(FIXED), precision(ADDITIONAL_WIDTH), rows(rows), cols(cols) {
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
    int tmp = 0;
    is >> rows >> cols >> precision >> tmp;
    style = static_cast<DoubleStyle>(tmp);
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
    if (precision > DISPLAY_WIDTH - ADDITIONAL_WIDTH) {
        throw BadMatrixPrecision();
    }
    size_t width = precision + ADDITIONAL_WIDTH;
    switch (style) {
        case SCIENTIFIC: width += ADDITIONAL_WIDTH; break;
        case FIXED: break;
        default: throw BadMatrixStyle();
    }
    const size_t display_row_cnt = DISPLAY_WIDTH / width;
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

Matrix::Matrix(size_t rows, size_t cols, byte *mem) : style(FIXED), precision(ADDITIONAL_WIDTH), rows(rows), cols(cols) {
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
    style = rhs.style;
    precision = rhs.precision;
    memcpy(raw_data, rhs.raw_data, rows * cols * sizeof(double));
}

Matrix &Matrix::operator=(const Matrix &rhs) {
    if (this == &rhs) {
        return *this;
    }
    delete[] data;
    delete[] raw_data;
    style = rhs.style;
    precision = rhs.precision;
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
    bool zero_matrix = true;
    for (size_t row_i = 0; row_i < cols; ++row_i) {
        for (size_t col_i = 0; col_i < rows; ++col_i) {
            if (!almostEqual(0, data[row_i][col_i])) {
                ++use_row[row_i];
                ++use_col[col_i];
                zero_matrix = false;
            }
        }
    }
    if (zero_matrix) {
        return;
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
    repeat_cnt = 0;
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
        throw DimensionMismatch(*this);
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
        throw DimensionMismatch(*this);
    }
    double determinant = det();
    if (almostEqual(0, determinant)) {
        throw SingularMatrix();
    }
    return adj() * (1 / det());
}

Matrix Matrix::invGauJor() const {
    if (rows != cols) {
        throw DimensionMismatch(*this);
    }
    Matrix original(*this);
    Matrix inverted(rows, cols);
    inverted.fillE();
    double tmp = 0;
    original.removeZerosFromMain(inverted);
    #ifdef LAB10_MSGS
    std::cout << "After removing zeros from main diag:\n" << original << "(original)\n";
    #endif
    for (size_t row_i = 0; row_i < inverted.rows; ++row_i) {
        tmp = original.data[row_i][row_i];
        if (almostEqual(0, tmp)) {
            throw SingularMatrix(); // removeZerosFromMain() should prevent this, but what if something went wrong...
        }
        for (size_t col_i = inverted.cols - 1; col_i != static_cast<size_t>(-1); --col_i) {
            inverted.data[row_i][col_i] /= tmp;
            original.data[row_i][col_i] /= tmp;
        }
        for (size_t row_j = row_i + 1; row_j < inverted.rows; ++row_j) {
            inverted.rowMinusRowCoef(row_j, original.data[row_j][row_i], row_i);
            original.rowMinusRowCoef(row_j, original.data[row_j][row_i], row_i);
        }
    }
    #ifdef LAB10_MSGS
    std::cout << "After forward move:\n" << original << "(original)\n" << inverted << "(inverted)\n";
    #endif
    for (size_t row_i = inverted.rows - 1; row_i != static_cast<size_t>(-1); --row_i) {
        for (size_t row_j = 0; row_j < row_i; ++row_j) {
            inverted.rowMinusRowCoef(row_j, original.data[row_j][row_i], row_i);
            original.rowMinusRowCoef(row_j, original.data[row_j][row_i], row_i);
        }
    }
    #ifdef LAB10_MSGS
    std::cout << "After backward move:\n" << original << "(original)\n" << inverted << "(inverted)\n";
    #endif
    return inverted;
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

void Matrix::rowMinusRowCoef(size_t row_j, double coef, size_t row_i) {
    for (size_t col_i = 0; col_i < cols; ++col_i) {
        data[row_j][col_i] -= data[row_i][col_i] * coef;
    }
}

void Matrix::removeZerosFromMain(Matrix &additional) {
    bool *if_row_not_zero_in_col = new bool[rows * cols];
    size_t *cnts = new size_t[rows];
    double **new_data = new double *[rows];
    double **new_data_additional = new double *[rows];
    size_t min_pos = 0;
    size_t min_cnt = rows + 1;
    for (size_t row_i = 0; row_i < rows; ++row_i) {
        cnts[row_i] = 0;
        for (size_t col_i = 0; col_i < cols; ++col_i) {
            if (!almostEqual(0, data[row_i][col_i])) {
                if_row_not_zero_in_col[row_i * cols + col_i] = true;
                ++cnts[row_i];
            } else {
                if_row_not_zero_in_col[row_i * cols + col_i] = false;
            }
        }
    }
    /*for (size_t row_i = 0; row_i < rows; ++row_i) {
        std::cout << cnts[row_i] << " | ";
        for (size_t col_i = 0; col_i < cols; ++col_i) {
            std::cout << if_row_not_zero_in_col[row_i * cols + col_i];
        }
        std::cout << std::endl;
    }*/
    size_t col_i_added = 0;
    bool is_singular = false;
    for (size_t row_i_adding = 0; row_i_adding < rows; ++row_i_adding) {
        min_cnt = rows + 1;
        is_singular = true;
        for (size_t row_i = 0; row_i < rows; ++row_i) {
            if (cnts[row_i] < min_cnt && if_row_not_zero_in_col[row_i * cols + row_i_adding]) {
                min_pos = row_i;
                min_cnt = cnts[row_i];
                is_singular = false;
            }
        }
        if (is_singular) {
            throw SingularMatrix();
        }
        new_data[row_i_adding] = data[min_pos];
        new_data_additional[row_i_adding] = additional.data[min_pos];
        col_i_added = row_i_adding;
        cnts[min_pos] = static_cast<size_t>(-1);
        for (size_t row_i = 0; row_i < rows; ++row_i) {
            if (if_row_not_zero_in_col[row_i * cols + col_i_added]) {
                if_row_not_zero_in_col[row_i * cols + col_i_added] = false;
                --cnts[row_i];
            }
        }
        //std::cout << "AFTER move #" << row_i_adding << ":\n" << *this << "-----\n";
    }
    delete[] data;
    data = new_data;
    delete[] additional.data;
    additional.data = new_data_additional;
    delete[] if_row_not_zero_in_col;
    delete[] cnts;
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
    if (matrix.precision > DISPLAY_WIDTH - ADDITIONAL_WIDTH) {
        throw BadMatrixPrecision();
    }
    os << matrix.rows << ' ' << matrix.cols << ' '
       << matrix.precision << ' ' << matrix.style << '\n';
    auto old_flags = os.flags();
    auto old_precision = os.precision(matrix.precision);
    size_t width = matrix.precision + ADDITIONAL_WIDTH;
    switch (matrix.style) {
        case SCIENTIFIC: width += ADDITIONAL_WIDTH; os.setf(std::ios::scientific); break;
        case FIXED: os.setf(std::ios::fixed); break;
        default: throw BadMatrixStyle();
    }
    const size_t display_row_cnt = DISPLAY_WIDTH / width;
    const size_t parts = (matrix.cols - 1) / display_row_cnt + 1;
    for (size_t part_i = 0; part_i < parts; ++part_i) {
        for (size_t row_i = 0; row_i < matrix.rows; ++row_i) {
            size_t end = (part_i != parts - 1) ? (part_i + 1) * display_row_cnt : matrix.cols;
            for (size_t col_i = part_i * display_row_cnt; col_i < end; ++col_i) {
                os << std::setw(width) << matrix.data[row_i][col_i];
            }
            os << '\n';
        }
        if (part_i + 1 < parts) {
            os << '\n';
        }
    }
    os.precision(old_precision);
    os.flags(old_flags);
    return os;
}
