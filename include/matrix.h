#pragma once  // NOLINT

#include <istream>

#define EPSILON 10e-5

typedef char byte;  // so byte[] is not cstring, it is for big memory part

constexpr size_t MIN_RECUR_SQR_MTRX = 2;
constexpr size_t ONE_ELEM_MTRX_ROWS_CNT = 1;
constexpr double BAD_RESULT = 0;
constexpr size_t ADDITIONAL_WIDTH = 5;
constexpr size_t DISPLAY_WIDTH = 79;

enum Sign {
    MINUS = -1,
    PLUS = 1
};

enum DoubleStyle {
    SCIENTIFIC,
    FIXED
};

class Matrix {
    DoubleStyle style;
    size_t precision;
    size_t rows;
    size_t cols;
    double** data;     // so we can use obj_ptr->data[i][j]
    double* raw_data;  // without memory fragmentation

    explicit Matrix(size_t rows, size_t cols, byte *mem);  // Placement Matrix
                                                        // for det optimization

    double detWithRaws(char* mem_tmps_raws) const;  // Using in det

    void createMinorPart(size_t row_to_del, size_t col_to_del,  // Using in adj
                         Matrix *result) const;           // and in detWithRaws

    void rowMinusRowCoef(size_t row_j, double coef, size_t row_i);

    void removeZerosFromMain(Matrix &additional);

 public:
    void setPrecision(size_t new_precision);

    void setStyle(DoubleStyle new_style);

    void fillSpecial();

    void fillE();

    explicit Matrix(size_t rows = 0, size_t cols = 0);

    explicit Matrix(std::istream &is);

    Matrix(const Matrix &rhs);

    Matrix &operator=(const Matrix &rhs);

    ~Matrix();

    size_t getRows() const;

    size_t getCols() const;

    double operator()(size_t i, size_t j) const;

    double &operator()(size_t i, size_t j);

    bool operator==(const Matrix &rhs) const;

    bool operator!=(const Matrix &rhs) const;

    friend
    Matrix sumWithAdditionalSign(const Matrix &lhs, const Matrix &rhs, Sign sign);

    Matrix operator+(const Matrix &rhs) const;

    Matrix operator-(const Matrix &rhs) const;

    Matrix operator*(const Matrix &rhs) const;

    Matrix operator*(double val) const;

    friend
    Matrix operator*(double val, const Matrix &matrix);

    friend
    std::ostream &operator<<(std::ostream &os, const Matrix &matrix);

    void deleteZeroStreaks();

    void trunc(size_t rows, size_t cols);

    Matrix transp() const;

    double det() const;  // on my pc ExtraTest passed for ~18.5 sec after its
                         // optimization
    Matrix adj() const;

    Matrix inv() const;

    Matrix invGauJor() const;
};

Matrix sumWithAdditionalSign(const Matrix &lhs, const Matrix &rhs, Sign sign);

Matrix operator*(double val, const Matrix &matrix);

std::ostream &operator<<(std::ostream &os, const Matrix &matrix);
