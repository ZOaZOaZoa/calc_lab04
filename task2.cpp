#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
#include <algorithm>

#include "funcs.h"

double f(double t, double y, double gamma)
{
    return (y*(1-y) - gamma*y/0.3);
}

double t0 = 0;
const double T = 50;
double calc_V(double gamma);

int main()
{
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

    //errors[i] <=> data[i+1]
    //V[i] <=> data[i+1]
    std::vector<double> errors(gamma_steps);
    size_t double_step = (steps-1)/2 + 1;
    std::vector<double> rung_values(double_step + 1);
    
    std::vector<double> V(gamma_steps);

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
    }

    double gamma_max = opt_golden_ratio(calc_V, 0, 1, 5e-4, true);
    double V_max = calc_V(gamma_max);

    std::cout << "----------------------------------------------\n";
    print_to_file(data, 18, output_file);
    std::cout << "Calculations finished. See " << output_file << " for results\n";
    print_errors(errors, gamma0, gamma_h, V);

    std::cout << "Max V is " << V_max << " reached by gamma = " << gamma_max << '\n';

    return 0;
}

double calc_V(double gamma)
{
    const double steps = 501;
    std::vector<double> solution(steps);
    solution[0] = 0.1;

    solve_EK(solution, f, steps-1, t0, T, gamma);
    return gamma * integrate_simpson(solution, T);
}