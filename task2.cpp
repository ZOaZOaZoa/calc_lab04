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

void print_to_file(std::vector<std::vector<double>>& data, size_t precision, std::string file_name)
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

void print_errors(std::vector<double>& errors, double gamma0, double gamma_h, std::vector<double>& V)
{
    std::cout << std::setw(10) << "gamma" << std::setw(20) << "error" << std::setw(15) << "V" << std::endl;
    for(int i = 0; i < errors.size(); i++)
    {
        std::cout << std::setw(10) << std::fixed << std::setprecision(6) << (gamma0 + i*gamma_h)
                << std::setw(20) << std::scientific << std::setprecision(3) << errors[i] 
                << std::setw(15) << std::fixed << std::setprecision(10) << V[i] << '\n';
    }
}

double integrate_simpson(std::vector<double>& y, double T)
{
    size_t N = y.size();
    double h = T / (N - 1);
    double steps = (N - 1)/2;

    double res = 0;
    for(int i = 1; i < steps; i++)
    {
        res += 4 * y[2*i - 1];
        res += 2 * y[2*i];
    }
    res -= 2 * y[N - 1];

    res += y[0];
    res += y[N-1];

    return h/6*res;
}

//3 знака
int main()
{
    double t0 = 0;
    const double T = 50;
    const int steps = 501;
    const double h = (T-t0)/(steps-1);
    double gamma0 = 0;
    //const double gamma_steps = 11;
    const int gamma_steps = 10001;
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

    //errors[i] <=> data[i+1]
    //V[i] <=> data[i+1]
    std::vector<double> errors(gamma_steps);
    size_t double_step = (steps-1)/2 + 1;
    std::vector<double> rung_values(double_step + 1);
    
    std::vector<double> V(gamma_steps);

    double gamma_max = 0;
    double V_max = 0;
    for(int i = 0; i < gamma_steps; i++)
    {
        double current_gamma = gamma0 + i*gamma_h;
        solve_EK(data[i+1], f, steps-1, t0, T, current_gamma);
        solve_EK(rung_values, f, double_step, t0, T, current_gamma);

        errors[i] = (data[i+1][0] - rung_values[0]) / 3;
        for(int j = 1; j < double_step; j++)
        {
            double cur_error = (data[i+1][2*j] - rung_values[j]) / 3;
            if (cur_error > errors[i])
            {
                errors[i] = cur_error;
            }
        }
        
        V[i] = current_gamma * integrate_simpson(data[i+1], T);
        if(V[i] > V_max)
        {
            V_max = V[i];
            gamma_max = current_gamma;
        }
    }
    std::cout << "----------------------------------------------\n";
    //print_results(y, t0, h);
    print_to_file(data, 18, output_file);
    std::cout << "Calculations finished. See " << output_file << " for results\n";
    print_errors(errors, gamma0, gamma_h, V);

    std::cout << "Max V is " << V_max << " reached by gamma = " << gamma_max << '\n';

    return 0;
}