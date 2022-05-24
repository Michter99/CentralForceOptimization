#pragma warning(disable : 6262)
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

#define DIMENSIONES 2      // Nd
#define SONDAS 100          // Np
#define ITERACIONES 10      // Nt
#define ALFA 1
#define BETA 2
#define GAMMA 3

using namespace std;

double ackley(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double a = 20.0f;
    double b = 0.2f;
    double c = (double)2 * M_PI;
    double d = DIMENSIONES;
    double sumatoria1 = 0;
    double sumatoria2 = 0;
    for (int i = 0; i < d; i++) {
        sumatoria1 += pow(posiciones[p][i][j], 2);
    }
    for (int i = 0; i < d; i++) {
        sumatoria2 += cos(c * posiciones[p][i][j]);
    }
    return -((double)(-a * (exp(-b * sqrt(((sumatoria1 / d)))))) - (exp((sumatoria2 / d))) + a + (double)exp(1));
}

double griewank(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double sumatoria = 0;
    double multiplicacion = 1;
    for (int i = 0; i < DIMENSIONES; i++) {
        sumatoria += pow(posiciones[p][i][j], 2);
    }
    for (int i = 1; i <= DIMENSIONES; i++) {
        multiplicacion *= cos(posiciones[p][i - 1][j] / sqrt(i));
    }
    return -((sumatoria / 4000) + multiplicacion + 1);
}

double rastrigin(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double sumatoria = 0;
    for (int i = 0; i < DIMENSIONES; i++) {
        sumatoria += pow(posiciones[p][i][j], 2) - 10 * cos(2 * M_PI * posiciones[p][i][j]);
    }
    return -(10 * DIMENSIONES + sumatoria);
}

double rotatedHyperEllipsoid(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double sumatoria = 0;
    for (int i = 0; i < DIMENSIONES; i++) {
        for (int k = 0; k < i + 1; k++) {
            sumatoria += pow(posiciones[p][k][j], 2);
        }
    }
    return -(sumatoria);
}

double sphere(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double sumatoria = 0;
    for (int i = 0; i < DIMENSIONES; i++) {
        sumatoria += pow(posiciones[p][i][j], 2);
    }
    return -sumatoria;
}

double dixonPrice(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double sumatoria = 0;
    for (int i = 1; i < DIMENSIONES; i++) {
        sumatoria += (i + 1) * pow((2 * pow(posiciones[p][i][j], 2) - posiciones[p][i - 1][j]), 2);
    }
    return -(pow(posiciones[p][0][j] - 1, 2) + sumatoria);
}

double sumSquares(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double sumatoria = 0;
    for (int i = 0; i < DIMENSIONES; i++) {
        sumatoria += (i + 1) * pow(posiciones[p][i][j], 2);
    }
    return -(sumatoria);
}

double trid(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double sumatoria1 = 0;
    double sumatoria2 = 0;
    for (int i = 0; i < DIMENSIONES; i++) {
        sumatoria1 += pow(posiciones[p][i][j] - 1, 2);
    }
    for (int i = 1; i < DIMENSIONES; i++) {
        sumatoria2 += posiciones[p][i][j] * posiciones[p][i - 1][j];
    }
    return -(sumatoria1 - sumatoria2);
}

double rosenbrock(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double sumatoria = 0;
    for (int i = 0; i < DIMENSIONES - 1; i++) {
        sumatoria += (100 * pow(posiciones[p][i + 1][j] - pow(posiciones[p][i][j], 2), 2) + pow(posiciones[p][i][j] - 1, 2));
    }
    return -(sumatoria);
}

double levy(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    double w[DIMENSIONES];
    double sumatoria = 0;
    for (int i = 0; i < DIMENSIONES; i++) {
        w[i] = 1 + ((posiciones[p][i][j] - 1) / 4);
    }
    for (int i = 0; i < DIMENSIONES - 1; i++) {
        sumatoria += pow(w[i] - 1, 2) * (1 + 10 * pow(sin(M_PI * w[i] + 1), 2));
    }
    return -(pow(sin(M_PI * w[0]), 2) + sumatoria + (pow(w[DIMENSIONES - 1] - 1, 2) * (1 + pow(sin(2 * M_PI * w[DIMENSIONES - 1]), 2))));
}

double test(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], int p, int j) {
    return -(sin(posiciones[p][0][j]) * exp(0.1 * posiciones[p][0][j]));
}

int funcionHeaviside(double val) {
    if (val >= 0)
        return 1;
    return 0;
}

void inicializarVariables(double posiciones[SONDAS][DIMENSIONES][ITERACIONES], double aceleraciones[SONDAS][DIMENSIONES][ITERACIONES], double xMin[DIMENSIONES], double xMax[DIMENSIONES]) {
    double deltaXi = 0;

    for (int p = 0; p < SONDAS; p++) {
        for (int i = 0; i < DIMENSIONES; i++) {
            deltaXi = (xMax[i] - xMin[i]) / (SONDAS - 1);
            posiciones[p][i][0] = xMin[i] + (p - 1) * deltaXi;
        }
    }

    for (int p = 0; p < SONDAS; p++) {
        for (int i = 0; i < DIMENSIONES; i++) {
            aceleraciones[p][i][0] = 0;
        }
    }
}

void centralForceOptimization(double xMin[DIMENSIONES], double xMax[DIMENSIONES], double(*func)(double[SONDAS][DIMENSIONES][ITERACIONES], int, int)) {
    double posiciones[SONDAS][DIMENSIONES][ITERACIONES];    // R
    double aceleraciones[SONDAS][DIMENSIONES][ITERACIONES]; // A
    double masas[SONDAS][ITERACIONES];                      // M
    double mejorMasa = -INFINITY;
    int mejorIteracion = 0;
    int mejorSonda = 0;
    double sumSQ = 0;
    double denominador = 0;
    double numerador = 0;

    // Inicializar posiciones y aceleraciones
    inicializarVariables(posiciones, aceleraciones, xMin, xMax);

    // Calcular masas iniciales
    for (int p = 0; p < SONDAS; p++) {
        masas[p][0] = func(posiciones, p, 0); // cambiar función
        if (masas[p][0] >= mejorMasa) {
            mejorMasa = masas[p][0];
            mejorIteracion = 0;
            mejorSonda = p;
        }
    }

    // Mientras que j < Nt
    for (int j = 1; j < ITERACIONES; j++) {

        // Actualizar la posición de la sonda
        for (int p = 0; p < SONDAS; p++) {
            for (int i = 0; i < DIMENSIONES; i++) {
                posiciones[p][i][j] = posiciones[p][i][j - 1] + (aceleraciones[p][i][j - 1] / 2);
            }
        }

        // Verificar que cada posición sea válida y corregir si aplica
        for (int p = 0; p < SONDAS; p++) {
            for (int i = 0; i < DIMENSIONES; i++) {
                if (posiciones[p][i][j] < xMin[i])
                    posiciones[p][i][j] = xMin[i] + (((posiciones[p][i][j - 1] - xMin[i]) / 2));
                if (posiciones[p][i][j] > xMax[i])
                    posiciones[p][i][j] = xMax[i] - (((xMax[i] - posiciones[p][i][j - 1]) / 2));
            }
        }

        // Actualizar masas con nuevas posiciones
        for (int p = 0; p < SONDAS; p++) {
            masas[p][j] = func(posiciones, p, j); // cambiar función
            if (masas[p][j] > mejorMasa) { // si se encuentra una mejor masa (evaluación) se actualiza la mejor masa, iteración y sonda
                mejorMasa = masas[p][j];
                mejorIteracion = j;
                mejorSonda = p;
            }
        }

        // Actualizar aceleraciones con los valores de posiciones y masas actualizados
        for (int p = 0; p < SONDAS; p++) {
            for (int i = 0; i < DIMENSIONES; i++) {
                aceleraciones[p][i][j] = 0;
                for (int k = 0; k < SONDAS; k++) {
                    if (k != p) {
                        sumSQ = 0;
                        for (int L = 0; L < DIMENSIONES; L++)
                            sumSQ += pow(posiciones[k][L][j] - posiciones[p][L][j], 2);
                        if (sumSQ != 0) {
                            denominador = sqrt(sumSQ);
                            numerador = funcionHeaviside(masas[k][j] - masas[p][j]) * (masas[k][j] - masas[p][j]);
                            aceleraciones[p][i][j] += ((posiciones[k][i][j] - posiciones[p][i][j]) * (pow(numerador, ALFA) / pow(denominador, BETA)));
                        }
                    }
                }
            }
        }
    }

    cout << "Mejor iteración: " << mejorIteracion << endl << endl;
    for (int i = 0; i < DIMENSIONES; i++)
        cout << "x" << i + 1 << ": " << posiciones[mejorSonda][i][mejorIteracion] << endl;
    cout << endl << "f(xi) = " << -mejorMasa << endl;
}


int main() {
    double min[DIMENSIONES];
    double max[DIMENSIONES];


    /**** ackley ****/
    //fill_n(min, DIMENSIONES, -32.768);
    //fill_n(max, DIMENSIONES, 32.768);

    /**** griewank ****/
    //fill_n(min, DIMENSIONES, -600);
    //fill_n(max, DIMENSIONES, 600);

    /**** rastrigin ****/
    //fill_n(min, DIMENSIONES, -5.12);
    //fill_n(max, DIMENSIONES, 5.12);

    /**** rotatedHyperEllipsoid ****/
    //fill_n(min, DIMENSIONES, -65.536);
    //fill_n(max, DIMENSIONES, 65.536);

    /**** sphere ****/
    //fill_n(min, DIMENSIONES, -5.12);
    //fill_n(max, DIMENSIONES, 5.12);

    /**** dixonPrice ****/
    //fill_n(min, DIMENSIONES, -10);
    //fill_n(max, DIMENSIONES, 10);

    /**** sumSquares ****/
    //fill_n(min, DIMENSIONES, -10);
    //fill_n(max, DIMENSIONES, 10);

    /**** trid ****/
    //fill_n(min, DIMENSIONES, -pow(DIMENSIONES, 2));
    //fill_n(max, DIMENSIONES, pow(DIMENSIONES, 2));

    /**** rosenbrock ****/
    //fill_n(min, DIMENSIONES, -5);
    //fill_n(max, DIMENSIONES, 10);

    /**** levy ****/
    fill_n(min, DIMENSIONES, -10);
    fill_n(max, DIMENSIONES, 10);


    // Pasar la dirección de la funcion a llamar (opciones comentadas arriba)
    centralForceOptimization(min, max, &levy);
}
