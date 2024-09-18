#ifndef INTERVAL_H
#define INTERVAL_H

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>

class Interval {
public:
    double lower;
    double upper;

    // Default constructor
    Interval() : lower(0.0), upper(0.0) {}

    // Constructor initializing the lower and upper bounds of the interval
    Interval(double l, double u) : lower(l), upper(u) {}

    // Negating the interval
    Interval operator-() const {
        return Interval(-upper, -lower);
    }

    // Adding intervals
    Interval operator+(const Interval& other) const {
        return Interval(lower + other.lower, upper + other.upper);
    }

    // Subtracting intervals
    Interval operator-(const Interval& other) const {
        return Interval(lower - other.lower, upper - other.upper);
    }

    // Multiplying intervals
    Interval operator*(const Interval& other) const {
        double ll = lower * other.lower;
        double lu = lower * other.upper;
        double ul = upper * other.lower;
        double uu = upper * other.upper;
        return Interval(std::min({ll, lu, ul, uu}), std::max({ll, lu, ul, uu}));
    }

    // Dividing intervals
    Interval operator/(const Interval& other) const {
        if (other.lower <= 0 && other.upper >= 0) {
            throw std::domain_error("Division by interval containing zero.");
        }
        double ll = lower / other.lower;
        double lu = lower / other.upper;
        double ul = upper / other.lower;
        double uu = upper / other.upper;
        return Interval(std::min({ll, lu, ul, uu}), std::max({ll, lu, ul, uu}));
    }

    // Function to check if the interval contains zero
    bool containsZero() const {
        return lower <= 0 && upper >= 0;
    }

    // Raising the interval to an integer power
    Interval power(int exp) const {
        if (exp == 0) {
            return Interval(1, 1);
        }
        Interval result = *this;
        for (int i = 1; i < exp; ++i) {
            result = result * *this;
        }
        return result;
    }

    // Outputting the interval
    friend std::ostream& operator<<(std::ostream& os, const Interval& interval) {
        os << "[" << interval.lower << ", " << interval.upper << "]";
        return os;
    }
};

#endif // INTERVAL_H
