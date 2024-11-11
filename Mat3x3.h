#ifndef MAT3X3_H
#define MAT3X3_H

#include <array>
#include <iostream>
#include <stdexcept>

/**
 * @class Mat3x3
 * @brief A class to represent a 3x3 matrix and perform various matrix operations.
 */
class Mat3x3 {
public:
    /**
     * @brief Default constructor. Initializes all nine entries to zero.
     */
    Mat3x3();

    /**
     * @brief Custom constructor. Accepts initial values for all nine entries.
     * @param list Initializer list containing nine elements.
     */
    Mat3x3(std::initializer_list<double> list);

    // Default copy, move constructors, and assignment operators
    Mat3x3(const Mat3x3&) = default;
    Mat3x3(Mat3x3&&) = default;
    Mat3x3& operator=(const Mat3x3&) = default;
    Mat3x3& operator=(Mat3x3&&) = default;

    // Overloaded operators
    Mat3x3& operator+=(const Mat3x3& rhs);
    Mat3x3& operator-=(const Mat3x3& rhs);
    Mat3x3& operator*=(const Mat3x3& rhs);
    Mat3x3& operator/=(const Mat3x3& rhs);
    Mat3x3& operator+=(double scalar);
    Mat3x3& operator-=(double scalar);
    Mat3x3& operator*=(double scalar);
    Mat3x3& operator/=(double scalar);

    Mat3x3& operator++(); // Prefix increment
    Mat3x3 operator++(int); // Postfix increment
    Mat3x3& operator--(); // Prefix decrement
    Mat3x3 operator--(int); // Postfix decrement
    Mat3x3 operator+() const;
    Mat3x3 operator-() const;
    bool operator!() const; // Logical negation

    double& operator()(size_t row, size_t col);
    const double& operator()(size_t row, size_t col) const;

    // Member functions
    double determinant() const;
    double trace() const;
    bool isAntisymmetric() const;
    bool isOrthogonal() const;
    bool isInvertible() const;
    bool isSymmetric() const;
    Mat3x3 transpose() const;
    Mat3x3 inverse() const;

    // Overloaded conversion operator
    explicit operator bool() const;

    // Function object
    double operator()() const;

    // Friends
    friend std::ostream& operator<<(std::ostream& os, const Mat3x3& mat);
    friend std::istream& operator>>(std::istream& is, Mat3x3& mat);

    // Accessor functions
    double get(size_t row, size_t col) const;
    void set(size_t row, size_t col, double value);

private:
    std::array<std::array<double, 3>, 3> values; ///< 2D array to store matrix values
};

// Non-member operators
Mat3x3 operator+(Mat3x3 lhs, const Mat3x3& rhs);
Mat3x3 operator-(Mat3x3 lhs, const Mat3x3& rhs);
Mat3x3 operator*(Mat3x3 lhs, const Mat3x3& rhs);
Mat3x3 operator/(Mat3x3 lhs, const Mat3x3& rhs);
Mat3x3 operator+(Mat3x3 lhs, double scalar);
Mat3x3 operator-(Mat3x3 lhs, double scalar);
Mat3x3 operator*(Mat3x3 lhs, double scalar);
Mat3x3 operator/(Mat3x3 lhs, double scalar);
Mat3x3 operator+(double scalar, Mat3x3 rhs);
Mat3x3 operator-(double scalar, Mat3x3 rhs);
Mat3x3 operator*(double scalar, Mat3x3 rhs);
Mat3x3 operator/(double scalar, Mat3x3 rhs);

bool operator==(const Mat3x3& lhs, const Mat3x3& rhs);
bool operator!=(const Mat3x3& lhs, const Mat3x3& rhs);
bool operator<(const Mat3x3& lhs, const Mat3x3& rhs);
bool operator<=(const Mat3x3& lhs, const Mat3x3& rhs);
bool operator>(const Mat3x3& lhs, const Mat3x3& rhs);
bool operator>=(const Mat3x3& lhs, const Mat3x3& rhs);

#endif // MAT3X3_H
