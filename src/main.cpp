#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "exceptions.h"
#include "matrix.h"

constexpr size_t PRECISION = 5;
constexpr size_t NAME_SIZE = 50;

enum SolveStatus {
    OK,
    CANT
};

static SolveStatus printSolveDSLAE(const Matrix &A, const Matrix &B) {
    std::cout << "Solving equation A * X = B:\n" << A << "*\n";
    size_t rows = A.getRows();
    size_t cols = A.getCols();
    for (size_t i = 0; i < rows; ++i) {
        std::cout << "  x" << i + 1 << '\n';
    }
    std::cout << "=\n";
    std::cout << B
              << "Finding A^-1:\n";
    Matrix inverted_A(rows, cols);
    try {
        inverted_A = A.invGauJor();
    } catch (SingularMatrix err) {
        std::cout << err.what() << " - no inverted matrix." << std::endl;
        return CANT;
    }
    std::cout << "A^-1 is:\n" << inverted_A;
    Matrix to_test(inverted_A * A);
    std::cout << "Testing, if A^-1 * A = E:\n" << to_test;
    Matrix testE(rows, cols);
    testE.fillE();
    std::cout << "Is it equal to E? - " << (to_test == testE)
              << "\nAnswer of the equation:\n" << inverted_A * B;
    return OK;
}

int main() {
    std::cout << "Running lab 10\n"
                 "For solving equation A * X = B enter file where is A from:"
              << std::endl;
    char A_file_name[NAME_SIZE];
    std::cin >> A_file_name;
    std::ifstream matrix_A(A_file_name);
    Matrix A(matrix_A);
    matrix_A.close();
    std::cout << "Enter file where is B from:" << std::endl;
    char B_file_name[NAME_SIZE];
    std::cin >> B_file_name;
    std::ifstream matrix_B(B_file_name);
    Matrix B(matrix_B);
    matrix_B.close();
    A.setPrecision(PRECISION);
    B.setPrecision(PRECISION);
    switch (printSolveDSLAE(A, B)) {
        case OK: std::cout << "Equation solved." << std::endl; break;
        case CANT: std::cout << "Can't solve equation." << std::endl; break;
        default: std::cout << "Incorrect solve status!" << std::endl;
    }
    return 0;
}
