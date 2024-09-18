#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "parser.cpp"
#include "parserInterval.cpp"
#include "interval.h"

#include <iostream>
#include <cmath>
#include <functional>
#include <QDoubleValidator>
#include <QIntValidator>

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

// Tworzenie funkcji pochodnych dla przedziałów
Interval partialDerivative(function<Interval(Interval, Interval)> func, const Interval& x1, const Interval& x2, bool respectToX1) {
    Interval h(1e-6, 1e-6);  // Krok dla różniczkowania numerycznego
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

// Definiowanie Jacobianu dla przedziałów
vector<vector<Interval>> jacobian(const vector<Interval>& x, const function<Interval(Interval, Interval)>& f1, const function<Interval(Interval, Interval)>& f2) {
    vector<vector<Interval>> J(2, vector<Interval>(2));
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

// Eliminacja Gaussa dla przedziałów
vector<Interval> gaussElimination(const vector<vector<Interval>>& A, const vector<Interval>& b) {
    int n = A.size();
    vector<vector<Interval>> augmentedMatrix(n, vector<Interval>(n + 1));

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
            if (fabs(augmentedMatrix[k][i].lower) > fabs(augmentedMatrix[maxRow][i].lower)) {
                maxRow = k;
            }
        }
        for (int k = i; k < n + 1; k++) {
            swap(augmentedMatrix[maxRow][k], augmentedMatrix[i][k]);
        }

        for (int k = i + 1; k < n; k++) {
            Interval c = -augmentedMatrix[k][i] / augmentedMatrix[i][i];
            for (int j = i; j < n + 1; j++) {
                if (i == j) {
                    augmentedMatrix[k][j] = Interval(0, 0);
                }
                else {
                    augmentedMatrix[k][j] = augmentedMatrix[k][j] + c * augmentedMatrix[i][j];
                }
            }
        }
    }

    vector<Interval> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = augmentedMatrix[i][n] / augmentedMatrix[i][i];
        for (int k = i - 1; k >= 0; k--) {
            augmentedMatrix[k][n] = augmentedMatrix[k][n] - augmentedMatrix[k][i] * x[i];
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

void MainWindow::on_resultButton_clicked()
{
    switch(currentMethod){
    case 1:{
        double initialGuess = ui->initialGuessInput->text().toDouble();
        auto array = splitByNewline(ui->textEdit->toPlainText().toStdString());

        Polynomial p1(array.at(0));
        Polynomial p2(array.at(1));

        function<double(double, double)> func1 = p1.getFunction();
        function<double(double, double)> func2 = p2.getFunction();

        vector<double> x = { initialGuess, initialGuess };  // Wstępne przybliżenia

        for (int iter = 0; iter <= maxIter; iter++) {
            vector<double> f = { func1(x[0], x[1]), func2(x[0], x[1]) };
            vector<vector<double>> J = jacobian(x, func1, func2);

            double normF = 0;
            for (double val : f) {
                normF += val * val;
            }
            normF = sqrt(normF);
            if (normF < tol) {
                QString str = "x = " + QString::number(x[0]) + ", y = " + QString::number(x[1]) + "\n\nZnaleziono w " + QString::number(iter) + " iteracjach"
                              + "\nEpsilon: " + QString::number(tol) + ", Max iteracji: " + QString::number(maxIter);
                ui->resultBox->setText(str);
                return;
            }

            vector<double> delta_x = gaussElimination(J, f);
            for (int i = 0; i < x.size(); i++) {
                x[i] -= delta_x[i];
            }
        }

        ui->resultBox->setText("Maximum number of iterations exceeded.");
    }break;

    case 2: {
        double initialLower = ui->initialGuessLowerInput->text().toDouble();
        double initialUpper = ui->initialGuessUpperInput->text().toDouble();
        Interval initialGuess(initialLower, initialUpper);

        auto array = splitByNewline(ui->textEdit->toPlainText().toStdString());
        int strPlacement = array.at(0).rfind('=');
        if (strPlacement == -1) strPlacement = array.at(0).size();
        array.at(0).resize(strPlacement);

        strPlacement = array.at(1).rfind('=');
        if (strPlacement == -1) strPlacement = array.at(1).size();
        array.at(1).resize(strPlacement);

        PolynomialInterval p1(array.at(0));
        PolynomialInterval p2(array.at(1));

        function<Interval(Interval, Interval)> func1 = p1.getFunctionInterval();
        function<Interval(Interval, Interval)> func2 = p2.getFunctionInterval();

        vector<Interval> x = { initialGuess, initialGuess };  // Wstępne przybliżenia

        for (int iter = 0; iter <= maxIter; iter++) {
            vector<Interval> f = { func1(x[0], x[1]), func2(x[0], x[1]) };
            vector<vector<Interval>> J = jacobian(x, func1, func2);

            if (f[0].containsZero() && f[1].containsZero()) {

                QString str = "x = [" + QString::number(x[0].lower) + ", " + QString::number(x[0].upper) + "], y = [" + QString::number(x[1].lower) + ", " + QString::number(x[1].upper) + "]\n\nZnaleziono w " + QString::number(iter) + " iteracjach"
                              + "\nEpsilon: " + QString::number(tol) + ", Max iteracji: " + QString::number(maxIter);
                ui->resultBox->setText(str);

                return;
            }

            vector<Interval> delta_x = gaussElimination(J, f);
            for (int i = 0; i < x.size(); i++) {
                x[i] = x[i] - delta_x[i];
            }
        }

        cout << "Maximum number of iterations exceeded." << endl;
    }break;
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
