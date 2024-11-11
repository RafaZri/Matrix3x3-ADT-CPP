#include "Mat3x3.h"
#include <iostream>
#include <cassert>

int main() {
    Mat3x3 im{1, 0, 0, 0, 1, 0, 0, 0, 1};
    Mat3x3 A{1, 1, -2, -3, -2, 5, -6, 4, 4};
    Mat3x3 m2{-28, -12, 1, -18, -8, 1, -24, -10, 1};
    Mat3x3 aInvSolution; // default ctor
    aInvSolution += 0.5 * m2; // op+=, op* (left scalar mult.)
    std::cout << "im: \n" << im << "\n"; // op<<
    std::cout << "A: \n" << A << "\n";
    std::cout << "aInvSolution: \n" << aInvSolution << "\n";
    assert(A * aInvSolution == im); // op*, op==
    Mat3x3 AinvComputed = A.inverse();
    std::cout << "AinvComputed: \n" << AinvComputed << "\n";
    assert(aInvSolution == AinvComputed); // inverse ok

    Mat3x3 B{-2, 5, -2, -2, 2, 0, -3, -2, 2};
    Mat3x3 m4{4, -6, 4, 4, -10, 4, 10, -19, 6};
    Mat3x3 bInvSolution = m4 * (-0.125); // op=, op* (right scalar mult.)
    std::cout << "B: \n" << B << "\n";
    std::cout << "bInvSolution: \n" << bInvSolution << "\n";
    assert(B * bInvSolution == im);
    Mat3x3 bInvComputed = B.inverse();
    std::cout << "bInvComputed: \n" << bInvComputed << "\n";
    assert(bInvSolution == bInvComputed); // inverse ok
    assert((A * B).inverse() == B.inverse() * A.inverse());

    Mat3x3 C{2, 7, 3, 1, 5, 8, 0, 4, 1};
    std::cout << "C: \n" << C << "\n";
    Mat3x3 D{3, 0, 1, 2, 1, 0, 1, 2, 4};
    std::cout << "D: \n" << D << "\n";
    Mat3x3 CxD{23, 13, 14, 21, 21, 33, 9, 6, 4};
    std::cout << "CxD: \n" << CxD << "\n";
    assert(C * D == CxD);

    std::cout << "CxD / D: \n" << (CxD / D) << "\n";
    assert(C == CxD / D);

    std::cout << "CxD / C: \n" << (CxD / C) << "\n";
    assert(C.inverse() * CxD == D);
    assert(--A == ++A); // unary + and -
    assert(-6 == (Mat3x3{1, 3, -2, 4, 1, -1, 5, -3, 2}.determinant()));
    assert(A.transpose().transpose() == A);
    assert(1.0 / (A.determinant()) == (A.inverse()).determinant());
    assert(A.transpose().determinant() == A.determinant());
    assert((A * B).transpose() == (B.transpose()) * (A.transpose()));

    Mat3x3 oldC{C};
    std::cout << "oldC: \n" << oldC << "\n";
    std::cout << "C: \n" << C << "\n";
    Mat3x3 Cpp{C++};
    std::cout << "Cpp: \n" << Cpp << "\n";
    std::cout << "C: \n" << C << "\n";
    assert(Cpp == oldC);
    assert(--C == oldC);
    std::cout << "C: \n" << C << "\n";
    std::cout << "oldC: \n" << oldC << "\n";

    return 0;
}
