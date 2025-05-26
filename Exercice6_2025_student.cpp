#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
typedef vector<complex<double>> vec_cmplx;

// Fonction resolvant le systeme d'equations A * solution = rhs
// où A est une matrice tridiagonale
template<class T>
void triangular_solve(vector<T> const& diag,  vector<T> const& lower, vector<T> const& upper,
                 vector<T> const& rhs, vector<T>& solution)
{
    vector<T> new_diag = diag;
    vector<T> new_rhs = rhs;

    // forward elimination
    for (size_t i(1); i < diag.size(); ++i) {
        T pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution.resize(diag.size());

    // solve last equation
    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    // backward substitution
    for (int i = diag.size() - 2; i >= 0; --i) {
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];
    }
}

// TODO Potentiel V(x) :
double V(double x, double xL, double xR, double xa, double xb, double V0, double om0) {
    if (x <= xa)
        return 0.5 * om0 * om0 * pow((x - xa) / (1. - xa / xL), 2);
    else if (x <= xb)
        return V0 * pow(sin(M_PI * (x - xa) / (xb - xa)), 2);
    else
        return 0.5 * om0 * om0 * pow((x - xb) / (1. - xb / xR), 2);
}

// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule dans un intervalle [x_i, x_j]
//  - E:    calcule son energie,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.

// TODO: calculer la probabilite de trouver la particule dans un intervalle [x_i, x_j]
double prob(vec_cmplx const& psi, int i, int j, double dx) {
    double sum = 0.5 * norm(psi[i]) + 0.5 * norm(psi[j]);
    for (int k = i + 1; k < j; ++k)
        sum += norm(psi[k]);
    return dx * sum;
}


// TODO calculer l'energie
double E(vec_cmplx const& psi, vector<double> const& x, double dx, double hbar, double m,
         double xL, double xR, double xa, double xb, double V0, double om0)
{
    double energy = 0.0;

    complex<double> d2psi = (psi[2] - 2. * psi[1] + psi[0]) / (dx * dx);
    complex<double> Hpsi = -0.5 * hbar * hbar / m * d2psi + V(x[1], xL, xR, xa, xb, V0, om0) * psi[1];
    energy += 0.5 * real(conj(psi[1]) * Hpsi);

    for (size_t i = 2; i < x.size() - 2; ++i)
    {
        complex<double> d2psi = (psi[i + 1] - 2. * psi[i] + psi[i - 1]) / (dx * dx);
        complex<double> Hpsi = -0.5 * hbar * hbar / m * d2psi + V(x[i], xL, xR, xa, xb, V0, om0) * psi[i];
        energy += real(conj(psi[i]) * Hpsi);
    }

    size_t mmm = x.size() - 2;
    complex<double> d2psi = (psi[mmm + 1] - 2. * psi[mmm] + psi[mmm - 1]) / (dx * dx);
    complex<double> Hpsi = -0.5 * hbar * hbar / m * d2psi + V(x[mmm], xL, xR, xa, xb, V0, om0) * psi[mmm];
    energy += 0.5 * real(conj(psi[mmm]) * Hpsi);

    return energy * dx;
}


// TODO calculer xmoyenne
double xmoy(vec_cmplx const& psi, vector<double> const& x, double dx) {
    complex<double> total = 0.;
    total += 0.5 * conj(psi[0]) * x[0] * psi[0];
    total += 0.5 * conj(psi.back()) * x.back() * psi.back();
    for (size_t i = 1; i < psi.size() - 1; ++i)
        total += conj(psi[i]) * x[i] * psi[i];
    return real(total) * dx;
}


// TODO calculer x.^2 moyenne
double x2moy(vec_cmplx const& psi, vector<double> const& x, double dx) {
    complex<double> total = 0.;
    total += 0.5 * conj(psi[0]) * x[0] * x[0] * psi[0];
    total += 0.5 * conj(psi.back()) * x.back() * x.back() * psi.back();
    for (size_t i = 1; i < psi.size() - 1; ++i)
        total += conj(psi[i]) * x[i] * x[i]* psi[i];
    return real(total) * dx;
}


// TODO calculer p moyenne
double pmoy(vec_cmplx const& psi, double dx, double hbar) {
    complex<double> sum = 0.;

    // Point gauche : dérivée avant
    complex<double> dpsi0 = (psi[1] - psi[0]) / dx;
    sum += conj(psi[0]) * dpsi0 * 0.5;

    // Centre : dérivée centrée
    for (size_t i = 1; i < psi.size() - 1; ++i) {
        complex<double> dpsi = (psi[i + 1] - psi[i - 1]) / (2. * dx);
        sum += conj(psi[i]) * dpsi;
    }

    // Point droit : dérivée arrière
    complex<double> dpsif = (psi.back() - psi[psi.size() - 2]) / dx;
    sum += conj(psi.back()) * dpsif * 0.5;

    return real(-hbar * complex<double>(0, 1) * dx * sum);
}


// TODO calculer p.^2 moyenne
double p2moy(vec_cmplx const& psi, double dx, double hbar) {
    complex<double> sum = 0.;

    // Point gauche : dérivée nulle imposée (psi = 0 aux bords)
    // donc ∂²ψ/∂x² = 0

    for (size_t i = 1; i < psi.size() - 1; ++i) {
        complex<double> d2psi = (psi[i + 1] - 2. * psi[i] + psi[i - 1]) / (dx * dx);
        sum += conj(psi[i]) * d2psi;
    }

    return real(-hbar * hbar * dx * sum);
}


// TODO calculer la normalization
vec_cmplx normalize(vec_cmplx const& psi, double const& dx) {
    double norm = 0;

    norm += 0.5 * std::norm(psi[0]) + 0.5 * std::norm(psi.back());
    for (size_t i = 1; i < psi.size() - 1; ++i)
        norm += std::norm(psi[i]);

    norm *= dx;

    vec_cmplx psi_norm(psi.size());
    for (size_t i = 0; i < psi.size(); ++i)
        psi_norm[i] = psi[i] / sqrt(norm);

    return psi_norm;
}

int main(int argc, char** argv)
{
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
    const double PI = 3.1415926535897932384626433832795028841971e0;

    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
        inputPath = argv[1];

    ConfigFile configFile(
      inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for (int i(2); i < argc;
         ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Parametres physiques :
    double hbar = 1.;
    double m = 1.;
    double tfin = configFile.get<double>("tfin");
    double xL = configFile.get<double>("xL");
    double xR = configFile.get<double>("xR");
    double xa = configFile.get<double>("xa");
    double xb = configFile.get<double>("xb");
    double V0 = configFile.get<double>("V0");
    double om0 = configFile.get<double>("om0");
    double n  = configFile.get<int>("n"); // Read mode number as integer, convert to double

    double x0 = configFile.get<double>("x0");
    double sigma0 = configFile.get<double>("sigma_norm") * (xR - xL);

    int Nsteps = configFile.get<int>("Nsteps");
    int Nintervals = configFile.get<int>("Nintervals");

    // TODO: initialiser le paquet d'onde, equation (4.116) du cours
    double k0 = 2 * PI * n / (xR-xL); 

    int Npoints = Nintervals + 1;
    double dx = (xR - xL) / Nintervals;
    double dt = tfin / Nsteps;

    const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    for (int i(0); i < Npoints; ++i)
        x[i] = xL + i * dx;

    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints);

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0 = floor((0 - xL)/(xR-xL)*Npoints); //chosen xR*0.5 since top of potential is at half x domain
  
    // TODO initialize psi
	for (int i = 0; i < Npoints; ++i) {
        double gauss = exp(-pow((x[i] - x0), 2) / (2 * sigma0 * sigma0));
        psi[i] = gauss * exp(complex_i * k0 * x[i]);
    }
    // Modifications des valeurs aux bords :
    psi[0] = complex<double>(0., 0.);
    psi[Npoints - 1] = complex<double>(0., 0.);
    
    // Normalisation :
    psi = normalize(psi, dx);

    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a = complex_i * hbar * dt / (4.*m*dx*dx); // Coefficient complexe a de l'equation (4.100)

    // TODO: calculer les éléments des matrices A, B et H.
    // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales
    // supérieures et inférieures
    for (int i(0); i < Npoints; ++i) // Boucle sur les points de maillage
    {
        double Vx = V(x[i], xL, xR, xa, xb, V0, om0);
        dH[i] = Vx;
        dA[i] = 1. + 2.0 * a + complex_i * dt * Vx / (2. * hbar);
        dB[i] = 1. - 2.0 * a - complex_i * dt * Vx / (2. * hbar);
    }
    for (int i(0); i < Nintervals; ++i) // Boucle sur les intervalles
    {
        aH[i] = cH[i] = -hbar * hbar / (2. * m * dx * dx);
        
        aA[i] = cA[i] = -a;
        aB[i] = cB[i] = a;
    }

    // Conditions aux limites: psi nulle aux deux bords
    // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites

	dA[0] = dA[Npoints - 1] = 1.;
    dB[0] = dB[Npoints - 1] = 1.;
    aA[0] = aA[Npoints - 2] = 0.;
    aB[0] = aB[Npoints - 2] = 0.;



    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " << V(x[i], xL, xR, xa, xb, V0, om0) << endl;
    fichier_potentiel.close();

    ofstream fichier_psi((output + "_psi2.out").c_str());
    fichier_psi.precision(6);

    ofstream fichier_observables((output + "_obs.out").c_str());
    fichier_observables.precision(15);

    // t0 writing
    for (int i(0); i < Npoints; ++i){
        fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << endl;
        }

    // Ecriture des observables :
    // TODO: introduire les arguments des fonctions prob, E, xmoy, x2moy, pmoy et p2moy
    //       en accord avec la façon dont vous les aurez programmés plus haut
    fichier_observables << t << " "
                        << prob(psi, 0, Npoints / 2, dx) << " "
                        << prob(psi, Npoints / 2, Npoints - 1, dx) << " "
                        << E(psi, x, dx, hbar, m, xL, xR, xa, xb, V0, om0) << " "
                        << xmoy(psi, x, dx) << " "
                        << x2moy(psi, x, dx) << " "
                        << pmoy(psi, dx, hbar) << " "
                        << p2moy(psi, dx, hbar) << endl;
    // Boucle temporelle :    
    while (t < tfin) {

        // Multiplication psi_tmp = B * psi :
        vec_cmplx psi_tmp(Npoints, 0.);
        for (int i(0); i < Npoints; ++i)
            psi_tmp[i] = dB[i] * psi[i];
        for (int i(0); i < Nintervals; ++i) {
            psi_tmp[i] += cB[i] * psi[i + 1];
            psi_tmp[i + 1] += aB[i] * psi[i];
        }

        // Resolution de A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        
        psi[0] = complex<double>(0., 0.);
        psi[Npoints - 1] = complex<double>(0., 0.);
        
        t += dt;

        // t0 writing
        for (int i(0); i < Npoints; ++i){
            fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << endl;
            }

        // Ecriture des observables :
	// TODO: introduire les arguments des fonctions prob, E, xmoy, x2moy, pmoy et p2moy
	//       en accord avec la façon dont vous les aurez programmés plus haut
        fichier_observables << t << " "
                            << prob(psi, 0, Npoints / 2, dx) << " "
                            << prob(psi, Npoints / 2, Npoints - 1, dx) << " "
                            << E(psi, x, dx, hbar, m, xL, xR, xa, xb, V0, om0) << " "
                            << xmoy(psi, x, dx) << " "
                            << x2moy(psi, x, dx) << " "
                            << pmoy(psi, dx, hbar) << " "
                            << p2moy(psi, dx, hbar) << endl;
    } // Fin de la boucle temporelle





    fichier_observables.close();
    fichier_psi.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}
