#include <iostream>
#include <cassert>
#include <matrix.hpp>

using namespace std;
using namespace sm::math;

int main() {
    std::cout << "============ Test 0 =================================" << endl;
    {
        Mat44d A({ {
            {1, 10, 8,  4},
            {8, 2,  5,  3},
            {1, 2,  3,  4},
            {1, 2,  3,  4}
        } }),
            B({ {
            {1, 10, 8,  4},
            {8, 2,  5,  3},
            {1, 2,  3,  4},
            {1, 2,  1,  4}
        } });

        try
        {
            std::cout << "A: " << A << endl
                << "B: " << B << endl
                << "A+B: " << A + B << endl
                << "A-B: " << A - B << endl
                << "A*B: " << A * B << endl;
        }
        catch (const std::exception& e)
        {
            cerr << e.what() << endl;
        }
    }
    std::cout << endl << endl << endl;

    std::cout << "============ Test 1 =================================" << endl;
    {
        Matrix<double, 2, 4> A({ {
            {1, 10, 8,  4},
            {1, 2,  3,  4}
        } });
        Matrix<double, 4, 6>
            B({ {
            {1, 10, 8,  4,  2,  5},
            {8, 2,  5,  3,  2,  1},
            {1, 2,  3,  4,  1,  5},
            {1, 2,  1,  4,  2,  3}
        } });

        try
        {
            std::cout << "A: " << A << endl
                << "B: " << B << endl
                << "A*B: " << A * B << endl
                << "A.transpose: " << A.transpose() << endl
                << "A.inversion: " << A.inversion() << endl;
        }
        catch (const std::exception& e)
        {
            cerr << e.what() << endl;
        }
    }
    std::cout << endl << endl << endl;

    std::cout << "============ Test 2 =================================" << endl;
    {
        Mat44d A({ {
            {1, 10, 8,  4},
            {8, 2,  5,  3},
            {1, 2,  3,  4},
            {1, 2,  3,  4}
        } }),
            B({ {
            {1, 10, 8,  4},
            {8, 2,  5,  3},
            {1, 2,  3,  4},
            {1, 2,  1,  4}
        } });

        Matrix<double, 4, 4> E(1);

        try
        {
            std::cout << "A: " << A << endl
                << "A.determinant: " << A.determinant() << endl << endl
                << "A.inversion: " << A.inversion() << endl;
        }
        catch (const std::exception& e)
        {
            cerr << e.what() << endl;
        }
        std::cout << "-------------------------------" << endl;
        try
        {
            std::cout << "B: " << B << endl
                << "B.determinant: " << B.determinant() << endl << endl
                << "E: " << E << endl
                << "B.attach(E): " << B.attach(E) << endl
                << "B.toDiagonalMatrix: " << B.attach(E).toCanonicalStepMatrix() << endl
                << "B.inversion: " << B.inversion() << endl
                << "B.inversion.inversion: " << B.inversion().inversion();
        }
        catch (const std::exception& e)
        {
            cerr << e.what() << endl;
        }
    }
    std::cout << endl << endl << endl;

    std::cout << "============ Test 3 =================================" << endl;
    {
        Mat44d A({ {
            {1, 10, 8,  4},
            {8, 2,  5,  3},
            {1, 2,  3,  4},
            {1, 2,  3,  4}
        } }),
            B({ {
            { 6,  3,  1,  3},
            { 2,  6,  9,  6},
            { 2,  1,  3,  0},
            { 1,  1,  3,  2}
        } });
        Matrix<double,4,1> Y({ {
            {1},
            {2},
            {3},
            {4}
        } });

        try
        {
            std::cout << "A: " << A << endl
                << "A.determinant: " << A.determinant() << endl << endl
                << "X = A.SLAU(Y): " << A.SLAU(Y) << endl;
        }
        catch (const std::exception& e)
        {
            cerr << e.what() << endl;
        }
        std::cout << "-------------------------------" << endl;
        try
        {
            std::cout << "B: " << B << endl
                << "B.determinant: " << B.determinant() << endl << endl
                << "Y: " << Y << endl
                << "B.attach(Y): " << B.attach(Y) << endl
                << "B.toStepMatrix: " << B.attach(Y).toStepMatrix() << endl
                << "B.toCanonicalStepMatrix: " << B.attach(Y).toCanonicalStepMatrix() << endl
                << "X = B.SLAU(Y): " << B.SLAU(Y) << endl;
        }
        catch (const std::exception& e)
        {
            cerr << e.what() << endl;
        }
    }

    return 0;
}