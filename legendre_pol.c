#include <stdio.h>
#include <math.h>

// Legendre P_n(x)'s polynomial using Bonnet's recursion
void legendre_polynomials(int n, double x, double *P) {
    P[0] = 1.0; // P_0(x) = 1
    if (n > 0) {
        P[1] = x; // P_1(x) = x
        for (int k = 2; k <= n; k++) {
            P[k] = ((2.0 * k - 1.0) * x * P[k-1] - (k - 1.0) * P[k-2]) / k;
        }
    }
}

void legendre_expansion(double (*f)(double), int n, double *coefficients) {
    for (int i = 0; i <= n; i++) {
        coefficients[i] = 0.0;
    }
    
    int num_samples = 1000;
    for (int j = 0; j < num_samples; j++) {
        double x = -1.0 + 2.0 * j / (num_samples - 1);
        double fx = f(x);
        double P[n+1];
        legendre_polynomials(n, x, P);
        for (int i = 0; i <= n; i++) {
            coefficients[i] += fx * P[i] * (2.0 / num_samples);
        }
    }
}

double reconstruct_function(double x, int n, double *coefficients) {
    double P[n+1];
    legendre_polynomials(n, x, P);
    double result = 0.0;
    for (int i = 0; i <= n; i++) {
        result += coefficients[i] * P[i];
    }
    return result;
}

double example_function(double x) {
    return cos(M_PI * x);
}

int main() {
    int n = 10; // max order
    double coefficients[n+1];

    legendre_expansion(example_function, n, coefficients);

    int n_of_reconstructions = 100; 
    for (int i = 0; i <= n_of_reconstructions; i++) {
        double x = -1.0 + 2.0 * i / n_of_reconstructions;
        double y = reconstruct_function(x, n, coefficients);
        printf("x = %f, y = %f\n", x, y);
    }

    return 0;
}
