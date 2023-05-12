#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

double y_precise(double t)
{
    return cos(t*t-1);
}

double f(double t, double y)
{
    return (y - 2*t*sin(t*t-1) - cos(t*t-1));
}

void solve_DE_rk(std::vector<double>& y, double f(double, double), double h, double t0, double T)
{
    double k1, k2, k3;
    size_t steps = (T-t0)/h;

    for(int i = 0; i < steps; i++)
    {
        k1 = h*f(t0 + i*h, y[i]);
        k2 = h*f(t0 + h*(i+0.5), y[i] + 0.5*k1);
        k3 = h*f(t0 + h*(i + 1), y[i] - k1 + 2*k2);
        y[i+1] = y[i] + 1.0/6*(k1 + 4*k2 + k3);
    }
}

double find_error(std::vector<double>& y, double y_presice(double), double t0, double h, int* max_error_index = nullptr)
{
    double error = 0;
    for (int i = 0; i < y.size(); i++)
    {
        double current_error = abs(y[i] - y_presice(t0 + i*h));
        if(current_error > error)
        {
            error = current_error;
            *max_error_index = i;
        }
    }

    return error;
}



void print_results(std::vector<double> y, double t0, double h)
{
    std::cout << std::setw(10) << "t" << std::setw(25) << "value\n";
    for(int i = 0; i <= y.size()-1; i++)
    {
        std::cout << std::setw(10) << std::setprecision(2) << std::scientific << (t0 + i*h)
            << std::setw(25) << std::setprecision(18) << std::fixed << y[i] << std::endl;
    }
}

int main()
{
    const double t0 = 1;
    const double T = 2;
    const double h = 0.1;
    const int steps = (T-t0)/h+1;
    std::vector<double> y(steps);
    y[0] = 1;

    std::cout << "Solving from on [" << t0 << ", " << T << "] with h = " << h 
        << ", steps: " << steps << std::endl;

    solve_DE_rk(y, f, h, t0, T);

    int max_error_index;
    double error = find_error(y, y_precise, t0, h, &max_error_index);

    print_results(y, t0, h);
    std::cout << "\nMax error: " << std::setprecision(2) << std::scientific << error << " on t= " << (t0 + max_error_index * h) << '\n';

    return 0;
}