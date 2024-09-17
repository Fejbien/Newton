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

    // Function to evaluate the polynomial for a given value of x
    double evaluate(double x) const {
        double result = 0.0;
        for (const auto& term : terms) {
            result += term.second * std::pow(x, term.first);
        }
        return result;
    }

    // Create a callable function for polynomial evaluation
    std::function<double(double)> getFunction() const {
        return [this](double x) {
            return this->evaluate(x);
        };
    }

    // Function to print the polynomial equation
    void print() const {
        bool firstTerm = true;
        for (auto it = terms.rbegin(); it != terms.rend(); ++it) {
            double coeff = it->second;
            int power = it->first;

            if (coeff == 0) continue;

            if (!firstTerm && coeff > 0) std::cout << " + ";
            if (coeff < 0) std::cout << " - ";

            if (std::abs(coeff) != 1 || power == 0) std::cout << std::abs(coeff);
            if (power > 0) std::cout << "x";
            if (power > 1) std::cout << "^" << power;

            firstTerm = false;
        }
        std::cout << std::endl;
    }

private:
    std::map<int, double> terms; // Stores terms as power -> coefficient

    // Function to parse the input polynomial equation
    void parseEquation(const std::string& equation) {
        // Regex that matches polynomial terms, including correct handling of negative coefficients
        std::regex termPattern(R"(([+-]?\s*\d*\.?\d*)(x(\^(\d+))?)?)");
        auto begin = std::sregex_iterator(equation.begin(), equation.end(), termPattern);
        auto end = std::sregex_iterator();

        // Clear terms to prevent unwanted accumulation
        terms.clear();

        for (auto it = begin; it != end; ++it) {
            std::smatch match = *it;
            std::string coeffStr = match[1].str();
            std::string varPart = match[2].str();
            std::string powerStr = match[4].str();

            // Remove spaces from coefficient string for correct conversion
            coeffStr.erase(std::remove_if(coeffStr.begin(), coeffStr.end(), ::isspace), coeffStr.end());

            // Skip if no valid match (like empty strings)
            if (coeffStr.empty() && varPart.empty()) {
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

            int power = 0;  // Default power is 0 for constants
            if (!varPart.empty()) {
                if (powerStr.empty()) {
                    power = 1;  // Implicit power of x (i.e., "x" means x^1)
                }
                else {
                    power = std::stoi(powerStr);  // Extract explicit power
                }
            }

            // Assign terms, avoiding accumulating multiple matches for the same power
            terms[power] = coefficient;
        }
    }
};
