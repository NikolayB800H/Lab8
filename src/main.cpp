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
        std::cout << "Runing lab №8\n";
        std::cout << std::setw(81) << std::setfill('=') << "[80 chars len]=\n" << std::setfill(' ');
        Matrix m3(SIZE, SIZE);
        m3.fillSpecial();
        std::cout << std::scientific << m3 << std::fixed
                  << std::setprecision(DOUBLE_WIDTH - ADDITIONAL_WIDTH) << m3;
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
        std::cout << std::resetiosflags(std::ios_base::floatfield) << m4;
        std::cout << B.B         << "  " << B.B[0]       << "  " << B.B[2]  << std::endl;
        std::cout << B.B[0][0]   << "  " << **B.B        << "  " << *B.B[0] << std::endl;
        std::cout << *(*(B.B+1)) << "  " << *B.B[1]      << std::endl;
        std::cout << *(B.B[0]+1) << "  " << *(*B.B+1)    << std::endl;
        std::cout << B.B[0][20]  << "  " << *(B.B[0]+20) << "  " << *B.B[2] << std::endl;
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
