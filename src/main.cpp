#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "exceptions.h"
#include "matrix.h"

constexpr size_t SIZE = 10;  // How many rows (and collums) in matrix (for lab8)

int main() {
    char yes_or_no = 'n';
    std::cout << "Check lab №8? (y/n)\n";
    std::cin >> yes_or_no;
    if (yes_or_no == 'y') {  // lab8
        std::cout << "Runing lab №8\nWhat precision?\n";
        size_t precision = 0;
        std::cin >> precision;
        size_t size_m = 0;
        std::cout << "What size?\n";
        std::cin >> size_m;
        std::cout << std::setw(81) << std::setfill('=') << "[80 chars len]=\n" << std::setfill(' ');
        Matrix m3(size_m, size_m);
        m3.fillSpecial();

        m3.setPrecision(precision);  // !!!
        //std::cout << "Test double: " << 1.1 << std::endl;
        m3.setStyle(SCIENTIFIC);  // !!!
        //std::cout << "Test double: " << 1.1 << std::endl;
        std::cout << m3;
        //std::cout << "Test double: " << 1.1 << std::endl;
        m3.setStyle(FIXED);  // !!!
        //std::cout << "Test double: " << 1.1 << std::endl;
        //std::cout << "Test double: " << 1.1 << std::endl;
        std::cout << m3;
        //std::cout << "Test double: " << 1.1 << std::endl;

        struct Bmaker {
            double *ptrs[SIZE];
            double B[SIZE][SIZE];
        } B;
        for (auto &i : B.B) {
            for (auto &j : i) {
                j = 0;
            }
        }
        Matrix m4(SIZE, SIZE, static_cast<byte *>(static_cast<void *>(B.ptrs)) - sizeof(Matrix));
        for (size_t row_i = 0; row_i < SIZE; ++row_i) {
            for (size_t col_i = 0; col_i < SIZE; ++col_i) {
                m4(row_i, col_i) = row_i * 10 + col_i;
            }
        }
        m4.setPrecision(0);  // !!!
        std::cout << m4;
        std::cout << B.B << std::endl;
        #define A B.B
        std::cout << A << std::endl;
        std::cout << "  " << A[2]+2       << "  " << *(A+1)  << std::endl;
        std::cout << (*A+2)  << "  " << *A[1] -1 << "  " << A[3][23] << std::endl;
        std::cout << (A+1) << "  " << A[7]-5      << std::endl;
        std::cout << A[0]+1 << "  " << *A[4]+4    << std::endl;
        std::cout << A[0][20]  << "  " << (A[0]+20) << "  " << A[2] << std::endl;
        m4.data = nullptr;
        m4.raw_data = nullptr;
    } else {  // delete zero streaks
        std::cout << "Runing no zero rows and cols\n";
        std::ifstream is("in");
        Matrix m1(is);
        std::cout << "Got from file:\n" << m1;
        m1.deleteZeroStreaks();
        std::cout << "Without zero streaks:\n" << m1;
        is.close();
    }
    return 0;
}
