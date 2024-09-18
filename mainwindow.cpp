#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "parser.cpp"

#include <iostream>
#include <cmath>
#include <functional>
#include <QDoubleValidator>
#include <QIntValidator>
#include <algorithm>

// Jakby, co zlego to nie ja :3
/*
testy dla:

x^2 + y - 1
-y^2-x

*/

double tol = 1e-6;            // Tolerancja dla zbieżności
int maxIter = 10000;          // Maksymalna liczba iteracji

using namespace std;
int currentMethod = 1;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QDoubleValidator *val = new QDoubleValidator();
    val->setNotation(QDoubleValidator::ScientificNotation);
    ui->epsilonEdit->setValidator(val);
    ui->maxIterEdit->setValidator(new QIntValidator());

    ui->epsilonEdit->setText(QString::number(tol));
    ui->maxIterEdit->setText(QString::number(maxIter));
}

MainWindow::~MainWindow()
{
    delete ui;
}

// Tworzenie funkcji pochodnych (Jacobian). Dla uproszczenia liczony numerycznie.
double partialDerivative(function<double(double, double)> func, double x1, double x2, bool respectToX1) {
    double h = 1e-6;  // Krok dla różniczkowania numerycznego
    if (respectToX1) {
        return (func(x1 + h, x2) - func(x1, x2)) / h;
    }
    else {
        return (func(x1, x2 + h) - func(x1, x2)) / h;
    }
}

// Definiowanie Jacobianu numerycznie
vector<vector<double>> jacobian(const vector<double>& x, const function<double(double, double)>& f1, const function<double(double, double)>& f2) {
    vector<vector<double>> J(2, vector<double>(2));
    J[0][0] = partialDerivative(f1, x[0], x[1], true);  // d(f1)/d(x1)
    J[0][1] = partialDerivative(f1, x[0], x[1], false); // d(f1)/d(x2)
    J[1][0] = partialDerivative(f2, x[0], x[1], true);  // d(f2)/d(x1)
    J[1][1] = partialDerivative(f2, x[0], x[1], false); // d(f2)/d(x2)
    return J;
}

// Rozwiązywanie układu równań liniowych metodą Gaussa (dla kroków metody Newtona)
vector<double> gaussElimination(const vector<vector<double>>& A, const vector<double>& b) {
    int n = A.size();
    vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));

    // Tworzenie macierzy rozszerzonej
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmentedMatrix[i][j] = A[i][j];
        }
        augmentedMatrix[i][n] = b[i];
    }

    // Eliminacja Gaussa
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(augmentedMatrix[k][i]) > fabs(augmentedMatrix[maxRow][i])) {
                maxRow = k;
            }
        }
        for (int k = i; k < n + 1; k++) {
            swap(augmentedMatrix[maxRow][k], augmentedMatrix[i][k]);
        }

        for (int k = i + 1; k < n; k++) {
            double c = -augmentedMatrix[k][i] / augmentedMatrix[i][i];
            for (int j = i; j < n + 1; j++) {
                if (i == j) {
                    augmentedMatrix[k][j] = 0;
                }
                else {
                    augmentedMatrix[k][j] += c * augmentedMatrix[i][j];
                }
            }
        }
    }

    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = augmentedMatrix[i][n] / augmentedMatrix[i][i];
        for (int k = i - 1; k >= 0; k--) {
            augmentedMatrix[k][n] -= augmentedMatrix[k][i] * x[i];
        }
    }
    return x;
}

vector<string> splitByNewline(const string& input) {
    vector<string> lines;
    stringstream ss(input);
    string line;
    while (getline(ss, line, '\n')) {
        lines.push_back(line);
    }
    return lines;
}

// Define a tolerance for comparing doubles
const double TOLERANCE = 0.5;

// Custom function to check if two vectors are "equal" considering floating-point precision
bool areVectorsEqual(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        return false;  // Vectors of different sizes cannot be equal
    }

    for (size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i] - b[i]) > TOLERANCE) {  // Compare with tolerance
            return false;  // If any element differs more than the tolerance, they are not equal
        }
    }
    return true;
}

void MainWindow::on_resultButton_clicked()
{
    switch(currentMethod){
    case 1: {
        // Case 1: Newton's method with doubles
        double initialGuess = ui->initialGuessInput->text().toDouble();
        auto array = splitByNewline(ui->textEdit->toPlainText().toStdString());

        Polynomial p1(array.at(0));
        Polynomial p2(array.at(1));

        function<double(double, double)> func1 = p1.getFunction();
        function<double(double, double)> func2 = p2.getFunction();

        vector<double> x = { initialGuess, initialGuess };  // Initial guesses

        for (int iter = 0; iter <= maxIter; iter++) {
            vector<double> f = { func1(x[0], x[1]), func2(x[0], x[1]) };
            vector<vector<double>> J = jacobian(x, func1, func2);

            double normF = 0;
            for (double val : f) {
                normF += val * val;
            }
            normF = sqrt(normF);
            if (normF < tol) {
                QString str = "x = " + QString::number(x[0]) + ", y = " + QString::number(x[1]) + "\n\nFound in " + QString::number(iter) + " iterations"
                              + "\nEpsilon: " + QString::number(tol) + ", Max iterations: " + QString::number(maxIter);
                ui->resultBox->setText(str);
                return;
            }

            vector<double> delta_x = gaussElimination(J, f);
            for (int i = 0; i < x.size(); i++) {
                x[i] -= delta_x[i];
            }
        }

        ui->resultBox->setText("Maximum number of iterations exceeded.");
    } break;

    case 2: {
        // Case 1: Newton's method with doubles
        string input = ui->initialGuessInput->text().toStdString();

        // Remove square brackets and semicolon
        input = input.substr(1, input.size() - 2);  // Removes the [ and ]

        // Replace the semicolon with a space for easier parsing
        size_t semicolonPos = input.find(';');
        input[semicolonPos] = ' ';

        double lower, upper;
        // Use stringstream to extract the two numbers
        std::stringstream ss(input);
        ss >> lower >> upper;

        cout << lower << upper << endl;

        auto array = splitByNewline(ui->textEdit->toPlainText().toStdString());

        Polynomial p1(array.at(0));
        Polynomial p2(array.at(1));

        function<double(double, double)> func1 = p1.getFunction();
        function<double(double, double)> func2 = p2.getFunction();

        vector<vector<double>> results = vector<vector<double>>();

        // Tak zgadza sie sprawdza sie jest to zjebane
        for(int i = lower * 10; i < upper * 10; i++){
            cout << i << endl;
            vector<double> x = { ((double)i)/10, ((double)i)/10 };  // Initial guesses

            for (int iter = 0; iter <= maxIter; iter++) {
                vector<double> f = { func1(x[0], x[1]), func2(x[0], x[1]) };
                vector<vector<double>> J = jacobian(x, func1, func2);

                double normF = 0;
                for (double val : f) {
                    normF += val * val;
                }
                normF = sqrt(normF);
                if (normF < tol) {
                    results.push_back({x[0], x[1]});  // Use an initializer list to push both x[0] and x[1]
                    break;
                }

                vector<double> delta_x = gaussElimination(J, f);
                for (int i = 0; i < x.size(); i++) {
                    x[i] -= delta_x[i];
                }
            }

            ui->resultBox->setText("Maximum number of iterations exceeded.");
        }

        // Create a new vector to store unique vectors
        std::vector<std::vector<double>> uniqueResults;

        // Iterate over the results and add unique vectors to uniqueResults
        for (const auto& vec : results) {
            bool isUnique = true;
            for (const auto& uniqueVec : uniqueResults) {
                if (areVectorsEqual(vec, uniqueVec)) {
                    isUnique = false;
                    break;
                }
            }
            if (isUnique) {
                uniqueResults.push_back(vec);
            }
        }

        QString str = "";

        for(int i = 0; i < uniqueResults.size(); i++){
            str += "\nx = " + QString::number(uniqueResults[i][0]) + ", y = " + QString::number(uniqueResults[i][1]);
        }

        ui->resultBox->setText(str);

    } break;
    }
}


QString initialGuessInputOld1 = "";
QString polynomialInputOld1 = "";
QString resultBoxOld1 = "";

QString initialGuessInputOld2 = "";
QString polynomialInputOld2 = "";
QString resultBoxOld2 = "";

void MainWindow::on_radioNetwon_clicked()
{
    initialGuessInputOld2 = ui->initialGuessInput->text();
    polynomialInputOld2 = ui->textEdit->toPlainText();
    resultBoxOld2 = ui->resultBox->toPlainText();

    ui->initialGuessInput->setText(initialGuessInputOld1);
    ui->textEdit->setText(polynomialInputOld1);
    ui->resultBox->setText(resultBoxOld1);

     ui->label_2->setText("Podaj zmienna");
    currentMethod = 1;
}

void MainWindow::on_radioInterval_clicked()
{
    initialGuessInputOld1 = ui->initialGuessInput->text();
    polynomialInputOld1 = ui->textEdit->toPlainText();
    resultBoxOld1 = ui->resultBox->toPlainText();

    ui->initialGuessInput->setText(initialGuessInputOld2);
    ui->textEdit->setText(polynomialInputOld2);
    ui->resultBox->setText(resultBoxOld2);

    ui->label_2->setText("Podaj przedzial np: [-2;3]");
    currentMethod = 2;
}

void MainWindow::on_epsilonEdit_textChanged(const QString &arg1)
{
    string str = arg1.toStdString();
    for (int i = 0; i < str.size(); i++)
        if(str[i] == ',') str[i] = '.';

    double dbl = 0.0;
    std::istringstream num(str);

    num >> dbl;

    if(!num.fail() && num.eof()) {
        tol = dbl;
    }
    else {
        cout << "Nie dobrze XD" << endl;
    }
}

void MainWindow::on_maxIterEdit_textChanged(const QString &arg1)
{
    maxIter = arg1.toInt();
}
