#include "funcs.h"

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

double opt_golden_ratio(double f(double), double a, double b, double eps, bool to_max)
{
    int sign = pow(-1, to_max);

    double alpha = a + 2 * (b-a)/(3 + sqrt(5));
    double beta = a + 2 * (b-a)/(1 + sqrt(5));
    double f_alpha = sign*f(alpha);
    double f_beta = sign*f(beta);

    while(abs(b-a)/2 > eps)
    {
        if(f_alpha > f_beta)
        {
            a = alpha;
            alpha = beta;
            f_alpha = f_beta;

            beta = a + 2 * (b-a)/(1 + sqrt(5));
            f_beta = sign*f(beta);
        }
        else
        {
            b = beta;
            beta = alpha;
            f_beta = f_alpha;
            
            alpha = a + 2 * (b-a)/(3 + sqrt(5));
            f_alpha = sign*f(alpha);
        }
    }

    return (b+a)/2;
}
