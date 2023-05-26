#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <algorithm>

#include "../lab03/int_functions.h"

double y_precise(double t)
{
    return cos(t*t-1);
}

double f(double t, double y)
{
    return (y - 2*t*sin(t*t-1) - cos(t*t-1));
}

void solve_DE_rk(std::vector<double>& y, double f(double, double), size_t steps, double t0, double T)
{
    double k1, k2, k3;
    double h = (T-t0)/(steps);

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

void print_to_file(std::vector<double> data, size_t precision, std::string file_name)
{
    std::ofstream output_file;
    output_file.open(file_name);
    output_file << std::setprecision(precision);
    for(auto val : data)
    {
        output_file << val << '\n';
    }
    output_file.close();
}

void print_to_file(std::vector<std::vector<double>> data, size_t precision, std::string file_name)
{
    std::ofstream output_file;
    output_file.open(file_name);
    output_file << std::setprecision(precision);
    for(int i = 0; i < data[0].size(); i++)
    {
        for(int j = 0; j < data.size(); j++)
        {
            output_file << data[j][i] << ' ';
        }
        output_file << '\n';
    }
    output_file.close();
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
    double t0 = 1;
    const double T = 2;
    const int steps = 10;
    const double h = (T-t0)/steps;
    std::vector<double> y(steps);
    y[0] = 1;

    std::cout << "Solving from on [" << t0 << ", " << T << "] with h = " << h 
        << ", steps: " << steps << std::endl;

    solve_DE_rk(y, f, steps, t0, T);

    int max_error_index;
    double error = find_error(y, y_precise, t0, h, &max_error_index);

    print_results(y, t0, h);
    
    std::vector<double> t(y.size());
    std::generate(t.begin(), t.end(), [h, a= (t0 - h)]() mutable {return a += h; });
    std::vector<std::vector<double>> data_to_print = {t, y};
    print_to_file(data_to_print, 18, "res.txt");
    std::cout << "\nMax error: " << std::setprecision(2) << std::scientific << error << " on t= " << (t0 + max_error_index * h) << '\n';

    return 0;
}