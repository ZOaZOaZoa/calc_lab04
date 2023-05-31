#ifndef FUNCS_H
#define FUNCS_H

#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <string.h>

void print_to_file(std::vector<std::vector<double>>& data, size_t precision, std::string file_name);
void solve_EK(std::vector<double>& y, double f(double, double, double), size_t steps, double t0, double T, double gamma);
void print_errors(std::vector<double>& errors, double gamma0, double gamma_h, std::vector<double>& V);
double integrate_simpson(std::vector<double>& y, double T);
double opt_golden_ratio(double f(double), double a, double b, double eps, bool to_max = false);


#endif //FUNCS_H