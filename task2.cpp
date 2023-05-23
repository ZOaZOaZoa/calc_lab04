#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <algorithm>

double f(double t, double y, double gamma)
{
    return (y*(1-y) - gamma*y/0.3);
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

void solve_EK(std::vector<double>& y, double f(double, double, double), size_t steps, double t0, double T, double gamma)
{
    double y_predicted;
    double h = (T-t0)/(steps);

    for(int i = 0; i < steps; i++)
    {
        y_predicted = y[i] + h*f(t0 + i*h, y[i], gamma);
        y[i+1] = y[i] + h * (f(t0 + i*h, y[i], gamma) + f(t0 + (i+1)*h, y_predicted, gamma))/2;
    }
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
    double t0 = 0;
    const double T = 50;
    const int steps = 501;
    const double h = (T-t0)/(steps-1);
    double gamma0 = 0;
    const int gamma_steps = 11;
    const double gamma_h = 1.0/(gamma_steps - 1);
    std::string output_file = "res.txt";
    
    std::vector<std::vector<double>> data(gamma_steps + 1);
    data[0] = std::vector<double>(steps);
    std::generate(data[0].begin(), data[0].end(), [h, a= (t0 - h)]() mutable {return a += h; });
    for(int i = 1; i < data.size(); i++)
    {
        data[i] = std::vector<double>(steps);
        data[i][0] = 0.1;
    }

    std::cout << "Solving from on [" << t0 << ", " << T << "] with h = " << h 
        << ", steps: " << steps << std::endl;

    for(int i = 0; i < gamma_steps; i++)
    {
        std::cout << (gamma0 + i*gamma_h) << '\n';
        solve_EK(data[i+1], f, steps-1, t0, T, gamma0 + i*gamma_h);
    }

    //print_results(y, t0, h);
    print_to_file(data, 18, output_file);
    std::cout << "Calculations finished. See " << output_file << " for results\n";

    return 0;
}