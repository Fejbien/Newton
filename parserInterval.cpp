#include <iostream>
#include <map>
#include <string>
#include <regex>
#include <cmath>
#include <functional>
#include "interval.h"  // Dodajemy nagłówek z definicją klasy Interval

class PolynomialInterval {
public:
    // Konstruktor, który przyjmuje równanie jako string
    PolynomialInterval(const std::string& equation) {
        parseEquation(equation);
    }

    // Funkcja do oceny wartości wielomianu dla podanych x i y (arytmetyka zmiennopozycyjna)
    double evaluate(double x, double y) const {
        double result = 0.0;
        for (const auto& term : terms) {
            int powerX = term.first.first;
            int powerY = term.first.second;
            result += term.second * std::pow(x, powerX) * std::pow(y, powerY);
        }
        return result;
    }

    // Funkcja do oceny wartości wielomianu dla przedziałów (arytmetyka przedziałowa)
    Interval evaluateAtInterval(const Interval& x, const Interval& y) const {
        Interval result(0.0, 0.0);
        for (const auto& term : terms) {
            int powerX = term.first.first;
            int powerY = term.first.second;

            // Obliczamy potęgowanie przedziału x i y
            Interval xPower = power(x, powerX);
            Interval yPower = power(y, powerY);

            // Mnożymy współczynnik przez x^powerX i y^powerY
            result = result + (xPower * yPower * Interval(term.second, term.second));
        }
        return result;
    }

    // Funkcja, która zwraca callable działający na arytmetyce zmiennopozycyjnej
    std::function<double(double, double)> getFunction() const {
        return [this](double x, double y) {
            return this->evaluate(x, y);
        };
    }

    // Funkcja, która zwraca callable działający na arytmetyce przedziałowej
    std::function<Interval(Interval, Interval)> getFunctionInterval() const {
        return [this](Interval x, Interval y) {
            return this->evaluateAtInterval(x, y);
        };
    }

    // Funkcja do wypisywania równania wielomianu
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
    // Mapa przechowująca składniki wielomianu jako (powerX, powerY) -> współczynnik
    std::map<std::pair<int, int>, double> terms;

    // Funkcja do parsowania równania wielomianu
    void parseEquation(std::string equation) {
        size_t equalSignPos = equation.find('=');
        std::string lhs = equation;
        std::string rhs;

        if (equalSignPos != std::string::npos) {
            lhs = equation.substr(0, equalSignPos);
            rhs = equation.substr(equalSignPos + 1);
            parseSide(rhs, true); // Parsujemy prawą stronę równania jako ujemne wartości
        }

        // Parsujemy lewą stronę
        parseSide(lhs, false);
    }

    // Funkcja pomocnicza do parsowania jednej strony równania
    void parseSide(const std::string& side, bool negate) {
        std::regex termPattern(R"(([+-]?\s*\d*\.?\d*)(x(\^(\d+))?)?(y(\^(\d+))?)?)");
        auto begin = std::sregex_iterator(side.begin(), side.end(), termPattern);
        auto end = std::sregex_iterator();

        for (auto it = begin; it != end; ++it) {
            std::smatch match = *it;
            std::string coeffStr = match[1].str();
            std::string xPart = match[2].str();
            std::string xPowerStr = match[4].str();
            std::string yPart = match[5].str();
            std::string yPowerStr = match[7].str();

            coeffStr.erase(std::remove_if(coeffStr.begin(), coeffStr.end(), ::isspace), coeffStr.end());

            if (coeffStr.empty() && xPart.empty() && yPart.empty()) continue;

            double coefficient = 1.0;
            if (!coeffStr.empty() && coeffStr != "+" && coeffStr != "-") {
                coefficient = std::stod(coeffStr);
            } else if (coeffStr == "-") {
                coefficient = -1.0;
            }

            if (negate) coefficient = -coefficient;

            int powerX = 0;
            if (!xPart.empty()) {
                if (xPowerStr.empty()) {
                    powerX = 1;
                } else {
                    powerX = std::stoi(xPowerStr);
                }
            }

            int powerY = 0;
            if (!yPart.empty()) {
                if (yPowerStr.empty()) {
                    powerY = 1;
                } else {
                    powerY = std::stoi(yPowerStr);
                }
            }

            terms[{powerX, powerY}] += coefficient;
        }
    }

    // Funkcja do potęgowania przedziału
    Interval power(const Interval& base, int exp) const {
        if (exp == 0) return Interval(1, 1);
        Interval result = base;
        for (int i = 1; i < exp; ++i) {
            result = result * base;
        }
        return result;
    }
};
