#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>// std::setprecision
#include "grid.h"

using namespace std;

#define M 100

double u1[M][M];
double u2[M][M];
double *u_new = (double *) u1;
double *u_old = (double *) u2;

void switch_arrays() {
    double *t;
    t = u_new;
    u_new = u_old;
    u_old = t;
}

int main() {
//    double e = 1e-5;
    int grid[100][100];
    getGrid(grid);
    double k = 0.5; // коэффициент диффузии
    double e = 1e-3; // точность
    double h = 0.01;
    double m = 100;
    double t = (h * h) / (4 * k);

    int i, j;

//    double *ptr;
//    ptr=(double*)u1;
//    (&ptr[100])[1] = 104;
//    std::cout << *(ptr)   << '\n';//out 0
//    std::cout << (&ptr[100])[1] << '\n';//out 1

    double a = 1 - (4 * t * k / (h * h));
    double a_bord = 1 - (2 * t * k / (h * h));
    double b = t * (k / (h * h) - 1 / (2 * h));
    double c = t * (k / (h * h) + 1 / (2 * h));
    double d = t * (k / (h * h));

    //cout << a << endl;
    //cout << b << endl;
    //cout << c << endl;
    //cout << d << endl;

    for (i = 0; i < M; i++) {
        for (j = 0; j < M; j++) {
            if (i == 0) {
                u_old[i*M+j] = 1;
            } else {
                u_old[i*M+j] = 0;
            }
        }
    }

    double u;
    int n = 0;
    double max_error = 1;
    //while(max_error > 1e-5) {
    while (n < 10000) {
        n++;
        double epoch_error = 0;
        for (i = 0; i < M; i++) {
            for (j = 0; j < M; j++) {
                if (grid[i][j] == 0) {
                    continue;
                }

                u = a * u_old[i*M+j];

                if (i > 0 && grid[i - 1][j] > 0) {
                    u += c * u_old[(i - 1)*M+j];
                }
                else {
                    u += b * u_old[i*M+j];
                    //u += a_bord * u_old[i*M + j];
                }

                if (i < M - 1 && grid[i + 1][j] > 0) {
                    u += b * u_old[(i + 1)*M+j];
                }
                else {
                    u += c * u_old[i*M+j];
                    //u += a_bord * u_old[i*M + j];
                }

                if (j > 0  && grid[i][j - 1] > 0) {
                    u += d * u_old[i*M+(j - 1)];
                }
                else {
                    u += d * u_old[i*M+j];
                    //u += a_bord * u_old[i*M + j];
                }

                if (j < M - 1  && grid[i][j + 1] > 0) {
                    u += d * u_old[i*M+(j + 1)];
                }
                else {
                    u += d * u_old[i*M+j];
                    //u += a_bord * u_old[i*M + j];
                }

                double cell_error = abs(u_old[i*M+j] - u);
                if (cell_error > epoch_error){
                    epoch_error = cell_error;
                }

                u_new[i*M+j] = u;
            }
        }

        max_error = epoch_error;

        if (n % 100 == 0) {
            cout << "epoch " << n << " error " << max_error << "\n";
        }

        switch_arrays();
    }

    cout << "epoch " << n << " error " << max_error << "\n";

    ofstream gridfile;
    gridfile.open("./result.txt");
    gridfile << std::setprecision(4);
    for (j = 99; j >= 0; j--) {
        for (i = 0; i < 100; i++) {
            //int res = u_old[i*M+j] == 0 ? 0 : 1;
            gridfile << u_old[i*M+j] << " ";
        }
        gridfile << "\n";
    }
    gridfile.close();
    return 0;
}