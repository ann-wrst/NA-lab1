// labachm.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <iomanip>
using namespace std;

double f1(double x) {
    return x * x * x - 5 * x * x - 4 * x + 20;
}
double df1(double x) {
    return 3 * x * x - 10 * x - 4;
}
double ddf1(double x) {
    return 6 * x - 10;
}
double f2(double x) {
    return x * x * x - 7 * x * x + 7 * x + 15;
}
double df2(double x) {
    return 3 * x * x - 14 * x + 7;
}
double f3(double x) {
    return x * x * x - 8 * x * x + 9 * x + 18;
}
double df3(double x) {
    return 3 * x * x - 16 * x + 9;
}
double ddf3(double x) {
    return 6 * x - 16;
}
double min(double a, double b, double (*f)(double)) {
    a *= 100.0;
    b *= 100.0;
    double min = (*f)(a/100.0);
    for (; a <= b; a += 1) {
        if ((*f)(a/100.0) < min)min = (*f)(a/100.0);
    }
    return min;
}

double max(double a, double b, double (*f)(double)) {
    a *= 100.0;
    b *= 100.0;
    double max = (*f)(a/100.0);
    for (; a <= b; a += 1) {
        if ((*f)(a/100.0) > max)max = (*f)(a/100.0);
    }
    return max;
}
void modnewton(double a,double b,double x_r, double x0, double eps = 1e-3) {
    cout << "Modified Newton\n" << "n\t\t\t\t\tx\t\t\t\t\t\tf(x)\n";
    if (a >= b) swap(a, b);
    if ((x0 > b) || (x0 < a)) {
        cout << "Out of range" << endl;
        return;
    }
    if (f1(a) * f1(b) >= 0 || f1(x0)*ddf1(x0)<=0) {
        cout << "Incorrect interval" << endl;
        return;
    }
    double m1 = abs(min(a, b, df1));
    double M2 = abs(min(a, b, ddf1));
    double q = (M2 * abs(x0 - x_r)) / (2 * m1);
    if (q >= 1) {
        cout << "Incorrect interval and x0" << endl;
        return;
    }
    int n = log((log(abs(x_r - x0) / eps) / log(1 / q)) + 1) / log(2) + 1;
    cout << "Number of iterations: n>=" << n << endl;
    int i = 0;
    cout << i << "\t\t\t\t\t" << x0 << "\t\t\t\t\t" << f1(x0) << endl;
    double x = x0 - f1(x0) / df1(x0);
    while (i<(n-1)) {
        i++;
        cout << i << "\t\t\t\t\t" << x << "\t\t\t\t\t" << f1(x) << endl;
        x = x - f1(x) / df1(x0);
    }
    cout << i+1 << "\t\t\t\t\t" << x << "\t\t\t\t\t" << f1(x) << endl;
    cout << "--------------------------------------------------------------------------------------------------------------" << endl;
}
void relaxation(double a, double b, double x0, double x_r, double eps = 1e-3) {
    cout<<"Relaxation\n" << "n\t\t\t\t\tx\t\t\t\t\t\tf(x)\n";
    if (a >= b) swap(a, b);
    if ((x0 > b)||(x0 < a)) {
        cout << "Out of range" << endl;
        return;
    }
    if (f2(a) * f2(b) >= 0) {
        cout << "Incorrect interval"<< endl;
        return;
    }
    double m1 = abs(min(a, b, df2));
    double M1 = abs(max(a, b, df2));

    double q = (M1 - m1) / (M1 + m1);
    int n = (log(eps * (1 - q)) / ((b - a)* log(q)));
    cout << "Number of iterations: n>=" << (n+1) << endl;
    int i = 0;
    cout << i << "\t\t\t\t\t" << x0 << "\t\t\t\t\t" << f2(x0) << endl;
    double tau = (double)2/(m1+M1);
    double x = x0 - tau * f2(x0);
    while (i < n){
        i++;
        cout << i << "\t\t\t\t\t" << x << "\t\t\t\t\t" << f2(x) << endl;
        x = round((x - tau * f2(x)) * 1000000) / 1000000;
    }
    cout << i + 1 << "\t\t\t\t\t" << x << "\t\t\t\t\t" << f2(x) << endl;
    cout << "--------------------------------------------------------------------------------------------------------------" << endl;
}

void secant(double a, double b,double x0, double x1, double x_r, double eps = 1e-3) {
    cout << "Secant method\n" << "n\t\t\t\t\tx\t\t\t\t\t\tf(x)\n";
    if (a >= b) swap(a, b);
    if ((x0 > b) || (x0 < a)) {
        cout << "Out of range" << endl;
        return;
    }
    if (f3(a) * f3(b) >= 0 || f3(x0) * ddf3(x0) <= 0) {
        cout << "Incorrect interval" << endl;
        return;
    }

    double m1 = abs(min(a, b, df3));
    double M2 = abs(min(a, b, ddf3));

    double q = (M2 * abs(x0 - x_r)) / (2 * m1);
    if (q >= 1) {
        cout << "Incorrect interval and x0" << endl;
        return;
    }
    int n = log((log(abs(x_r - x0) / eps) / log(1 / q)) + 1) / log(2) + 1;
    cout << "Number of iterations: n>=" << n << endl;
    int i = 0;
    cout << i << "\t\t\t\t\t" << x0 << "\t\t\t\t\t" << f3(x0) << endl;
    double x = x1 - ((x1 - x0) * f3(x1)) / (f3(x1) - f3(x0));
    double x_prev = x1;
    while (i < (n-1)) {
        i++;
        cout << i << "\t\t\t\t\t" << x << "\t\t\t\t\t" << f3(x) << endl;
        double x_temp=x;
        x = x - ((x - x_prev) * f3(x)) / (f3(x) - f3(x_prev));
        x_prev = x_temp;
    }
    cout << i + 1 << "\t\t\t\t\t" << x << "\t\t\t\t\t" << f3(x) << endl;
    cout << "--------------------------------------------------------------------------------------------------------------" << endl;
}

int main()
{
    cout << setprecision(5) << fixed;
    modnewton(-2.5,-1,-2, -2.5);
    relaxation(-2, 0, -2, -1);
    secant(-2,-0.5,-2, -1.5, -1);
    return 0;
}


// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
