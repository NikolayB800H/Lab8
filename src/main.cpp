#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "exceptions.h"
#include "matrix.h"

constexpr size_t PRECISION = 5;
constexpr size_t NAME_SIZE = 50;

int main() {
    std::cout << "Running lab 10\n"
                 "For solving equation A * X = B enter file where is A from:"
              << std::endl;
    char A_file_name[NAME_SIZE];
    std::cin >> A_file_name;
    std::cout << "Enter file where is B from:" << std::endl;
    char B_file_name[NAME_SIZE];
    std::cin >> B_file_name;
    std::ifstream matrix_A(A_file_name);
    Matrix A(matrix_A);
    matrix_A.close();
    std::cout << "Equation A * X = B:\n" << A << "*\n";
    size_t rows = A.getRows();
    size_t cols = A.getCols();
    for (size_t i = 0; i < rows; ++i) {
        std::cout << "  x" << i + 1 << '\n';
    }
    std::cout << "=\n";
    std::ifstream matrix_B(B_file_name);
    Matrix B(matrix_B);
    matrix_B.close();
    std::cout << B
              << "Finding A^-1:\n";
    A.setPrecision(PRECISION);
    B.setPrecision(PRECISION);
    Matrix inverted_A(A.invGauJor());
    std::cout << "A^-1 is:\n" << inverted_A;
    Matrix to_test(inverted_A * A);
    std::cout << "Testing, if A^-1 * A = E:\n" << to_test;
    Matrix testE(rows, cols);
    testE.fillE();
    std::cout << "Is it equal to E? - " << (to_test == testE)
              << "\nAnswer of the equation:\n" << inverted_A * B;
    return 0;
}
