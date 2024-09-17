#include "mainwindow.h"

#include <QApplication>

// int main(int argc, char *argv[])
// {
//     QApplication a(argc, argv);
//     MainWindow w;
//     w.show();
//     return a.exec();
// }


#include <iostream>
#include <cmath>
#include <functional>
#include "parser.cpp"  // Include the custom parser

// Function to calculate the numerical derivative of a function
std::function<double(double)> derivative(std::function<double(double)> f, double h = 1e-6) {
    return [f, h](double x) -> double {
        return (f(x + h) - f(x)) / h;
    };
}

// Function to perform the Newton-Raphson method
double newtonRaphson(std::function<double(double)> f, double x0, double epsilon = 1e-6, int maxIterations = 100) {
    double x = x0;
    for (int i = 0; i < maxIterations; i++) {
        double fx = f(x);       // Evaluate f(x)
        auto df = derivative(f);
        double dfx = df(x);     // Evaluate f'(x)

        if (std::abs(fx) < epsilon) {
            std::cout << "Converged to solution after " << i << " iterations.\n";
            return x;
        }

        if (dfx == 0) {
            std::cerr << "Derivative is zero, stopping.\n";
            return x;
        }

        double x1 = x - fx / dfx;  // Newton-Raphson iteration

        if (std::abs(x1 - x) < epsilon) {
            std::cout << "Converged to solution after " << i + 1 << " iterations.\n";
            return x1;
        }

        x = x1;  // Update x for the next iteration
    }
    std::cerr << "Maximum iterations reached.\n";
    return x;
}

// Example usage
int main(int argc, char *argv[]){


    // Define the function as a string (e.g., "x^3 - x^2 + 17")
    std::string fncstring = "x^3 - x^2 + 17";  // Example polynomial

    // Use the custom parser to parse the polynomial equation
    Polynomial poly(fncstring);

    // Get the function that evaluates the polynomial
    auto f = poly.getFunction();

    // Initial guess
    double initialGuess = 1.0;

    // Call the Newton-Raphson method to find the root
    double root = newtonRaphson(f, initialGuess);

    std::cout << "The root is: " << root << std::endl;

    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();

    return 0;
}
