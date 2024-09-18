#include <iostream>
#include <map>
#include <string>
#include <regex>
#include <cmath>
#include <functional>

class Polynomial {
public:
    // Constructor that takes an equation string
    Polynomial(const std::string& equation) {
        parseEquation(equation);
    }

    // Function to evaluate the polynomial for given values of x and y
    double evaluate(double x, double y) const {
        double result = 0.0;
        for (const auto& term : terms) {
            int powerX = term.first.first;
            int powerY = term.first.second;
            result += term.second * std::pow(x, powerX) * std::pow(y, powerY);
        }
        return result;
    }

    // Create a callable function for polynomial evaluation
    std::function<double(double, double)> getFunction() const {
        return [this](double x, double y) {
            return this->evaluate(x, y);
        };
    }

    // Function to print the polynomial equation
    void print() const {
        bool firstTerm = true;
        for (auto it = terms.rbegin(); it != terms.rend(); ++it) {
            double coeff = it->second;
            int powerX = it->first.first;
            int powerY = it->first.second;

            if (coeff == 0) continue;

            if (!firstTerm && coeff > 0) std::cout << " + ";
            if (coeff < 0) std::cout << " - ";

            if (std::abs(coeff) != 1 || (powerX == 0 && powerY == 0)) std::cout << std::abs(coeff);
            if (powerX > 0) std::cout << "x";
            if (powerX > 1) std::cout << "^" << powerX;
            if (powerY > 0) std::cout << "y";
            if (powerY > 1) std::cout << "^" << powerY;

            firstTerm = false;
        }
        std::cout << std::endl;
    }

private:
    // Map to store terms as (powerX, powerY) -> coefficient
    std::map<std::pair<int, int>, double> terms;

    // Function to parse the input polynomial equation
    void parseEquation(const std::string& equation) {
        // Regex to match polynomial terms involving x, y, or constants
        std::regex termPattern(R"(([+-]?\s*\d*\.?\d*)(x(\^(\d+))?)?(y(\^(\d+))?)?)");
        auto begin = std::sregex_iterator(equation.begin(), equation.end(), termPattern);
        auto end = std::sregex_iterator();

        // Clear terms to prevent unwanted accumulation
        terms.clear();

        for (auto it = begin; it != end; ++it) {
            std::smatch match = *it;
            std::string coeffStr = match[1].str();
            std::string xPart = match[2].str();
            std::string xPowerStr = match[4].str();
            std::string yPart = match[5].str();
            std::string yPowerStr = match[7].str();

            // Remove spaces from coefficient string for correct conversion
            coeffStr.erase(std::remove_if(coeffStr.begin(), coeffStr.end(), ::isspace), coeffStr.end());

            // Skip if no valid match (like empty strings)
            if (coeffStr.empty() && xPart.empty() && yPart.empty()) {
                continue;
            }

            double coefficient = 1.0;  // Default coefficient
            if (!coeffStr.empty() && coeffStr != "+" && coeffStr != "-") {
                coefficient = std::stod(coeffStr);  // Parse the coefficient
            }
            else if (coeffStr == "-") {
                coefficient = -1.0;  // Handle negative sign alone
            }
            else if (coeffStr == "+") {
                coefficient = 1.0;   // Handle positive sign alone
            }

            int powerX = 0;  // Default power for x is 0 for constants
            if (!xPart.empty()) {
                if (xPowerStr.empty()) {
                    powerX = 1;  // Implicit power of x (i.e., "x" means x^1)
                }
                else {
                    powerX = std::stoi(xPowerStr);  // Extract explicit power of x
                }
            }

            int powerY = 0;  // Default power for y is 0 for constants
            if (!yPart.empty()) {
                if (yPowerStr.empty()) {
                    powerY = 1;  // Implicit power of y (i.e., "y" means y^1)
                }
                else {
                    powerY = std::stoi(yPowerStr);  // Extract explicit power of y
                }
            }

            // Assign terms, avoiding accumulating multiple matches for the same power pair
            terms[{powerX, powerY}] = coefficient;
        }
    }
};
