#include <iostream>
#include <cmath>
#include <vector>

using namespace std;
const double PI = 3.14159265358979323846;

// Функція
double equation(double x) {
    return x * x + 5 * sin(x) - 1;
}

// Генерування вузлів Чебишова
vector<double> chebyshevNodes(int n, double a, double b) {
    vector<double> nodes;
    for (size_t i = 0; i < n; ++i) {
        nodes.push_back(((a + b) / 2) + (((b - a) / 2) * cos(((2.0 * i + 1.0) * PI) / (2.0 * n))));
    }
    return nodes;
}

// Функція для обчислення полінома Ньютона
double newtonInterpolation(const vector<double>& x, const vector<double>& y, double xi) {
    double result = 0.0;
    vector<double> coefficients(y.size(), 0.0);

    for (size_t k = 0; k < x.size(); k++) {
        for (size_t j = 0; j <= k; j++) {
            double temp = 1.0;
            for (size_t i = 0; i <= k; i++) {
                if (i != j) {
                    temp *= (x[j] - x[i]);
                }
            }
            coefficients[k] += y[j] / temp;
        }
    }

    for (size_t i = 0; i < coefficients.size(); ++i) {
        double term = 1;
        for (size_t j = 0; j < i; ++j) {
            term = term * (xi - x[j]);
        }
        result += term * coefficients[i];
    }

    return result;
}

// Функція для обчислення полінома Лагранжа
double lagrangeInterpolation(const vector<double>& x, const vector<double>& y, double xi) {
    double result = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double term = 1;
        for (size_t j = 0; j < x.size(); ++j) {
            if (j != i) {
                term = term * (xi - x[j]) / (x[i] - x[j]);
            }
        }
        result += term * y[i];
    }
    return result;
}

int main() {
    // Визначення проміжку та кількості вузлів
    double a = -3.0;
    double b = -1.0;
    int numNodes = 10;

    vector<double> x_values = chebyshevNodes(numNodes, a, b);
    vector<double> y_values;

    // Обчислення значень функції для відповідних x
    for (double x : x_values) {
        y_values.push_back(equation(x));
    }

    cout << "xk: {";
    for (size_t i = 0; i < x_values.size() - 1; i++) {
        cout << x_values[i] << ", ";
    }
    cout << x_values[x_values.size() - 1] << "}\n";

    cout << "yk: {";
    for (size_t i = 0; i < y_values.size() - 1; i++) {
        cout << y_values[i] << ", ";
    }
    cout << y_values[y_values.size() - 1] << "}\n";

    // Задання точки, в якій шукаємо значення функції
    // Тут були використані методи дихотомії, які були видалені

    return 0;
}