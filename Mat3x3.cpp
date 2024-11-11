#include "Mat3x3.h"
#include <cmath>
#include <iomanip>
#include <initializer_list>

// Constructors

Mat3x3::Mat3x3() : values{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}} {}

/**
 * @brief Custom constructor. Accepts initial values for all nine entries.
 * @param list Initializer list containing nine elements.
 */
Mat3x3::Mat3x3(std::initializer_list<double> list) {
    if (list.size() != 9) {
        throw std::invalid_argument("Initializer list must contain exactly 9 elements.");
    }

    auto it = list.begin();
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            values[i][j] = *it++;
        }
    }
}

// Overloaded operators

/**
 * @brief Adds the given matrix to this matrix.
 * @param rhs The matrix to add.
 * @return Reference to this matrix.
 */
Mat3x3& Mat3x3::operator+=(const Mat3x3& rhs) {
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            values[i][j] += rhs.get(i, j);
        }
    }
    return *this;
}

/**
 * @brief Subtracts the given matrix from this matrix.
 * @param rhs The matrix to subtract.
 * @return Reference to this matrix.
 */
Mat3x3& Mat3x3::operator-=(const Mat3x3& rhs) {
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            values[i][j] -= rhs.get(i, j);
        }
    }
    return *this;
}

/**
 * @brief Multiplies this matrix by the given matrix.
 * @param rhs The matrix to multiply by.
 * @return Reference to this matrix.
 */
Mat3x3& Mat3x3::operator*=(const Mat3x3& rhs) {
    Mat3x3 result;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            result.values[i][j] = 0;
            for (size_t k = 0; k < 3; ++k) {
                result.values[i][j] += values[i][k] * rhs.get(k, j);
            }
        }
    }
    *this = result;
    return *this;
}

/**
 * @brief Divides this matrix by the given matrix.
 * @param rhs The matrix to divide by.
 * @return Reference to this matrix.
 */
Mat3x3& Mat3x3::operator/=(const Mat3x3& rhs) {
    *this *= rhs.inverse();
    return *this;
}

/**
 * @brief Adds the given scalar to this matrix.
 * @param scalar The scalar to add.
 * @return Reference to this matrix.
 */
Mat3x3& Mat3x3::operator+=(double scalar) {
    for (auto& row : values) {
        for (auto& val : row) {
            val += scalar;
        }
    }
    return *this;
}

/**
 * @brief Subtracts the given scalar from this matrix.
 * @param scalar The scalar to subtract.
 * @return Reference to this matrix.
 */
Mat3x3& Mat3x3::operator-=(double scalar) {
    for (auto& row : values) {
        for (auto& val : row) {
            val -= scalar;
        }
    }
    return *this;
}

/**
 * @brief Multiplies this matrix by the given scalar.
 * @param scalar The scalar to multiply by.
 * @return Reference to this matrix.
 */
Mat3x3& Mat3x3::operator*=(double scalar) {
    for (auto& row : values) {
        for (auto& val : row) {
            val *= scalar;
        }
    }
    return *this;
}

/**
 * @brief Divides this matrix by the given scalar.
 * @param scalar The scalar to divide by.
 * @return Reference to this matrix.
 * @throws std::overflow_error if scalar is approximately zero.
 */
Mat3x3& Mat3x3::operator/=(double scalar) {
    if (std::abs(scalar) < 1e-10) {
        throw std::overflow_error("Division by zero");
    }
    for (auto& row : values) {
        for (auto& val : row) {
            val /= scalar;
        }
    }
    return *this;
}

/**
 * @brief Prefix increment. Increments all entries by 1.
 * @return Reference to this matrix.
 */
Mat3x3& Mat3x3::operator++() {
    *this += 1;
    return *this;
}

/**
 * @brief Postfix increment. Increments all entries by 1.
 * @return Copy of this matrix before increment.
 */
Mat3x3 Mat3x3::operator++(int) {
    Mat3x3 temp = *this;
    ++(*this);
    return temp;
}

/**
 * @brief Prefix decrement. Decrements all entries by 1.
 * @return Reference to this matrix.
 */
Mat3x3& Mat3x3::operator--() {
    *this -= 1;
    return *this;
}

/**
 * @brief Postfix decrement. Decrements all entries by 1.
 * @return Copy of this matrix before decrement.
 */
Mat3x3 Mat3x3::operator--(int) {
    Mat3x3 temp = *this;
    --(*this);
    return temp;
}

/**
 * @brief Unary plus operator.
 * @return Copy of this matrix.
 */
Mat3x3 Mat3x3::operator+() const {
    return *this;
}

/**
 * @brief Unary minus operator.
 * @return Copy of this matrix with all entries negated.
 */
Mat3x3 Mat3x3::operator-() const {
    Mat3x3 result;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            result.set(i, j, -values[i][j]);
        }
    }
    return result;
}

/**
 * @brief Logical negation operator.
 * @return True if the matrix is non-invertible, false otherwise.
 */
bool Mat3x3::operator!() const {
    return !isInvertible();
}

double& Mat3x3::operator()(size_t row, size_t col) {
    if (row >= 3 || col >= 3) {
        throw std::invalid_argument("index out of bounds");
    }
    return values[row][col];
}

const double& Mat3x3::operator()(size_t row, size_t col) const {
    if (row >= 3 || col >= 3) {
        throw std::invalid_argument("index out of bounds");
    }
    return values[row][col];
}

// Accessor functions

/**
 * @brief Gets the value at the specified row and column.
 * @param row The row index.
 * @param col The column index.
 * @return The value at the specified row and column.
 * @throws std::invalid_argument if row or col is out of bounds.
 */
double Mat3x3::get(size_t row, size_t col) const {
    if (row >= 3 || col >= 3) {
        throw std::invalid_argument("index out of bounds");
    }
    return values[row][col];
}

/**
 * @brief Sets the value at the specified row and column.
 * @param row The row index.
 * @param col The column index.
 * @param value The value to set.
 * @throws std::invalid_argument if row or col is out of bounds.
 */
void Mat3x3::set(size_t row, size_t col, double value) {
    if (row >= 3 || col >= 3) {
        throw std::invalid_argument("index out of bounds");
    }
    values[row][col] = value;
}

// Member functions

/**
 * @brief Calculates the determinant of the matrix.
 * @return The determinant of the matrix.
 */
double Mat3x3::determinant() const {
    return values[0][0] * (values[1][1] * values[2][2] - values[2][1] * values[1][2])
           - values[0][1] * (values[1][0] * values[2][2] - values[2][0] * values[1][2])
           + values[0][2] * (values[1][0] * values[2][1] - values[2][0] * values[1][1]);
}

/**
 * @brief Calculates the trace of the matrix.
 * @return The trace of the matrix.
 */
double Mat3x3::trace() const {
    return values[0][0] + values[1][1] + values[2][2];
}

/**
 * @brief Checks if the matrix is antisymmetric.
 * @return True if the matrix is antisymmetric, false otherwise.
 */
bool Mat3x3::isAntisymmetric() const {
    return values[0][1] == -values[1][0] &&
           values[0][2] == -values[2][0] &&
           values[1][2] == -values[2][1];
}

/**
 * @brief Checks if the matrix is orthogonal.
 * @return True if the matrix is orthogonal, false otherwise.
 */
bool Mat3x3::isOrthogonal() const {
    Mat3x3 t = transpose();
    Mat3x3 identity = (*this) * t;
    return identity == Mat3x3{1, 0, 0, 0, 1, 0, 0, 0, 1};
}

/**
 * @brief Checks if the matrix is invertible.
 * @return True if the matrix is invertible, false otherwise.
 */
bool Mat3x3::isInvertible() const {
    return std::abs(determinant()) > 1e-10;
}

/**
 * @brief Checks if the matrix is symmetric.
 * @return True if the matrix is symmetric, false otherwise.
 */
bool Mat3x3::isSymmetric() const {
    return values[0][1] == values[1][0] &&
           values[0][2] == values[2][0] &&
           values[1][2] == values[2][1];
}

/**
 * @brief Transposes the matrix.
 * @return The transposed matrix.
 */
Mat3x3 Mat3x3::transpose() const {
    Mat3x3 result;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            result.set(j, i, values[i][j]);
        }
    }
    return result;
}

/**
 * @brief Inverts the matrix.
 * @return The inverted matrix.
 * @throws std::overflow_error if the matrix is not invertible.
 */
Mat3x3 Mat3x3::inverse() const {
    double det = determinant();
    if (std::abs(det) < 1e-10) {
        throw std::overflow_error("Inverse undefined");
    }

    Mat3x3 result;
    result.set(0, 0, values[1][1] * values[2][2] - values[2][1] * values[1][2]);
    result.set(0, 1, values[0][2] * values[2][1] - values[0][1] * values[2][2]);
    result.set(0, 2, values[0][1] * values[1][2] - values[0][2] * values[1][1]);

    result.set(1, 0, values[1][2] * values[2][0] - values[1][0] * values[2][2]);
    result.set(1, 1, values[0][0] * values[2][2] - values[0][2] * values[2][0]);
    result.set(1, 2, values[0][2] * values[1][0] - values[0][0] * values[1][2]);

    result.set(2, 0, values[1][0] * values[2][1] - values[1][1] * values[2][0]);
    result.set(2, 1, values[0][1] * values[2][0] - values[0][0] * values[2][1]);
    result.set(2, 2, values[0][0] * values[1][1] - values[0][1] * values[1][0]);

    result /= det;
    return result;
}

// Overloaded conversion operator

/**
 * @brief Converts the matrix to bool.
 * @return True if the matrix is invertible, false otherwise.
 */
Mat3x3::operator bool() const {
    return isInvertible();
}

// Function object

/**
 * @brief Function call operator.
 * @return The determinant of the matrix.
 */
double Mat3x3::operator()() const {
    return determinant();
}

// Non-member operators

Mat3x3 operator+(Mat3x3 lhs, const Mat3x3& rhs) {
    lhs += rhs;
    return lhs;
}

Mat3x3 operator-(Mat3x3 lhs, const Mat3x3& rhs) {
    lhs -= rhs;
    return lhs;
}

Mat3x3 operator*(Mat3x3 lhs, const Mat3x3& rhs) {
    lhs *= rhs;
    return lhs;
}

Mat3x3 operator/(Mat3x3 lhs, const Mat3x3& rhs) {
    lhs /= rhs;
    return lhs;
}

Mat3x3 operator+(Mat3x3 lhs, double scalar) {
    lhs += scalar;
    return lhs;
}

Mat3x3 operator-(Mat3x3 lhs, double scalar) {
    lhs -= scalar;
    return lhs;
}

Mat3x3 operator*(Mat3x3 lhs, double scalar) {
    lhs *= scalar;
    return lhs;
}

Mat3x3 operator/(Mat3x3 lhs, double scalar) {
    lhs /= scalar;
    return lhs;
}

Mat3x3 operator+(double scalar, Mat3x3 rhs) {
    rhs += scalar;
    return rhs;
}

Mat3x3 operator-(double scalar, Mat3x3 rhs) {
    rhs -= scalar;
    return rhs;
}

Mat3x3 operator*(double scalar, Mat3x3 rhs) {
    rhs *= scalar;
    return rhs;
}

Mat3x3 operator/(double scalar, Mat3x3 rhs) {
    rhs /= scalar;
    return rhs;
}

/**
 * @brief Checks if two matrices are equal.
 * @param lhs The first matrix.
 * @param rhs The second matrix.
 * @return True if the matrices are equal, false otherwise.
 */
bool operator==(const Mat3x3& lhs, const Mat3x3& rhs) {
    const double epsilon = 1e-10;
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            if (std::abs(lhs.get(i, j) - rhs.get(i, j)) >= epsilon) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Checks if two matrices are not equal.
 * @param lhs The first matrix.
 * @param rhs The second matrix.
 * @return True if the matrices are not equal, false otherwise.
 */
bool operator!=(const Mat3x3& lhs, const Mat3x3& rhs) {
    return !(lhs == rhs);
}

/**
 * @brief Checks if one matrix is less than another.
 * @param lhs The first matrix.
 * @param rhs The second matrix.
 * @return True if the first matrix is less than the second, false otherwise.
 */
bool operator<(const Mat3x3& lhs, const Mat3x3& rhs) {
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            if (lhs.get(i, j) >= rhs.get(i, j)) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Checks if one matrix is less than or equal to another.
 * @param lhs The first matrix.
 * @param rhs The second matrix.
 * @return True if the first matrix is less than or equal to the second, false otherwise.
 */
bool operator<=(const Mat3x3& lhs, const Mat3x3& rhs) {
    return (lhs < rhs) || (lhs == rhs);
}

/**
 * @brief Checks if one matrix is greater than another.
 * @param lhs The first matrix.
 * @param rhs The second matrix.
 * @return True if the first matrix is greater than the second, false otherwise.
 */
bool operator>(const Mat3x3& lhs, const Mat3x3& rhs) {
    return !(lhs <= rhs);
}

/**
 * @brief Checks if one matrix is greater than or equal to another.
 * @param lhs The first matrix.
 * @param rhs The second matrix.
 * @return True if the first matrix is greater than or equal to the second, false otherwise.
 */
bool operator>=(const Mat3x3& lhs, const Mat3x3& rhs) {
    return !(lhs < rhs);
}

// Stream operators

std::ostream& operator<<(std::ostream& os, const Mat3x3& mat) {
    for (const auto& row : mat.values) {
        for (const auto& val : row) {
            os << std::setw(10) << val << " ";
        }
        os << "\n";
    }
    return os;
}

std::istream& operator>>(std::istream& is, Mat3x3& mat) {
    for (auto& row : mat.values) {
        for (auto& val : row) {
            is >> val;
        }
    }
    return is;
}
