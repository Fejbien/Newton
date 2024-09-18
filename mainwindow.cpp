#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "parser.cpp"
#include <iostream>
#include <cmath>
#include <functional>

// Jakby, co zlego to nie ja :3
/*
testy dla:

x^2 + y - 1
-y^2-x

*/

const double tol = 1e-6;            // Tolerancja dla zbieżności
const int maxIter = 10000;            // Maksymalna liczba iteracji

using namespace std;
int currentMethod = 1;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
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
        // Szukanie maksymalnego elementu w kolumnie
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(augmentedMatrix[k][i]) > fabs(augmentedMatrix[maxRow][i])) {
                maxRow = k;
            }
        }

        // Zamiana wierszy
        for (int k = i; k < n + 1; k++) {
            swap(augmentedMatrix[maxRow][k], augmentedMatrix[i][k]);
        }

        // Zerowanie elementów poniżej aktualnego wiersza
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

    // Rozwiązywanie układu równań
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
    stringstream ss(input);  // Create a stringstream object from the input string
    string line;

    // Use getline to split the string by '\n'
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
        // qDebug() << ui->textEdit->toPlainText();
        auto array = splitByNewline(ui->textEdit->toPlainText().toStdString());
        // Pozwala aby rownania konczyly sie na " = 0 "
        int strPlacement = array.at(0).rfind('=');
        if(strPlacement == -1) strPlacement = array.at(0).size();
        array.at(0).resize(strPlacement);

        strPlacement = array.at(1).rfind('=');
        if(strPlacement == -1) strPlacement = array.at(1).size();
        array.at(1).resize(strPlacement);

        // Parse the polynomials
        Polynomial p1(array.at(0));
        Polynomial p2(array.at(1));

        // Create callable functions from the parsed polynomials
        function<double(double, double)> func1 = p1.getFunction();
        function<double(double, double)> func2 = p2.getFunction();

        // Print the parsed polynomials
        // cout << "Parsed first equation: ";
        // p1.print();
        // cout << "Parsed second equation: ";
        // p2.print();

        // Initial guesses for the solution
        vector<double> x = { initialGuess, initialGuess };  // Wstępne przybliżenia


        // Newton's method loop
        for (int iter = 0; iter < maxIter; iter++) {
            // Calculate the function values
            vector<double> f = { func1(x[0], x[1]), func2(x[0], x[1]) };

            // Calculate the Jacobian matrix
            vector<vector<double>> J = jacobian(x, func1, func2);

            // Check convergence
            double normF = 0;
            for (double val : f) {
                normF += val * val;
            }
            normF = sqrt(normF);
            if (normF < tol) {
                // cout << "Solution found after " << iter << " iterations:" << endl;
                QString str = "x1 = " + QString::number(x[0]) + ", x2 = " + QString::number(x[1]) + "\n\nZnaleziono w " + QString::number(iter) + " iteracjach";
                ui->resultBox->setText(str);
                return;
            }

            // Calculate delta_x (Newton's step)
            vector<double> delta_x = gaussElimination(J, f);
            for (int i = 0; i < x.size(); i++) {
                x[i] -= delta_x[i];  // Update guesses
            }
        }

        ui->resultBox->setText("Maximum number of iterations exceeded.");

    }break;
    case 2:{
        //std::string str = ui->polynomialInput->text().toStdString();
        // Polynomial poly(str);
        // auto f = poly.getFunction();

        // double initialGuess = ui->initialGuessInput->text().toDouble();

        // double root = 0; // tutaj dasz inna funcke liczenia

        // ui->resultBox->setText(QString::number(root));
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

