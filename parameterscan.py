import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os
"""
# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = 'Ex2_2025_student'  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/home/chatelin/Desktop/MyFiles/1quantique_physNum"
os.chdir(repertoire)


CONFIG_FILE = os.path.join(os.path.dirname(__file__), "configuration.in.example")

input_filename = 'configuration.in.example'  # Name of the input file

def lire_configuration():
    config_path = os.path.join(os.path.dirname(__file__), "configuration.in.example")
    configuration = {}
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Le fichier {config_path} n'existe pas.")
    
    with open(config_path, "r", encoding="utf-8") as fichier:
        for ligne in fichier:
            ligne = ligne.strip()
            if ligne and "=" in ligne and not ligne.startswith("#"):
                cle, valeur = ligne.split("=", 1)
                configuration[cle.strip()] = valeur.strip()
    
    return configuration

def ecrire_configuration(nouvelles_valeurs):
    Écrit les nouvelles valeurs dans le fichier de configuration.
    if not os.path.exists(CONFIG_FILE):
        raise FileNotFoundError(f"Le fichier {CONFIG_FILE} n'existe pas.")

    lignes_modifiees = []
    
    with open(CONFIG_FILE, "r", encoding="utf-8") as fichier:
        for ligne in fichier:
            ligne_strippée = ligne.strip()
            if ligne_strippée and "=" in ligne_strippée and not ligne_strippée.startswith("#"):
                cle, _ = ligne_strippée.split("=", 1)
                cle = cle.strip()
                if cle in nouvelles_valeurs:
                    ligne = f"{cle} = {nouvelles_valeurs[cle]}\n"
            lignes_modifiees.append(ligne)

    with open(CONFIG_FILE, "w", encoding="utf-8") as fichier:
        fichier.writelines(lignes_modifiees)

Omega = 0.0
kappa = 0.0
m = 0.0
L = 0.0
B1 = 0.0
B0 = 0.0
mu = 0.0
theta0 = 0.0
thetadot0 = 0.0
sampling = 0.0
N_excit = 0.0
Nperiod = 0.0
nsteps = 0.0
C = 0.0
alpha = 0.0
beta = 0.0

valeurs = lire_configuration()

def actualise_valeur():
    global Omega, kappa, m, L, B1, B0, mu, theta0, thetadot0, sampling, N_excit, Nperiod, nsteps, C, alpha, beta
    Omega = float(valeurs.get("Omega"))
    kappa = float(valeurs.get("kappa"))
    m = float(valeurs.get("m"))
    L = float(valeurs.get("L"))
    B1 = float(valeurs.get("B1"))
    B0 = float(valeurs.get("B0"))
    mu = float(valeurs.get("mu"))
    theta0 = float(valeurs.get("theta0"))
    thetadot0 = float(valeurs.get("thetadot0"))
    sampling = float(valeurs.get("sampling"))
    N_excit = float(valeurs.get("N_excit"))
    Nperiod = float(valeurs.get("Nperiod"))
    nsteps = float(valeurs.get("nsteps"))
    C = float(valeurs.get("C"))
    alpha = float(valeurs.get("alpha"))
    beta = float(valeurs.get("beta"))

def ecrire_valeur(nom,valeur):
    global valeurs
    valeurs[nom] = valeur
    ecrire_configuration(valeurs)

actualise_valeur()


# Question 1

paramstr = 'nsteps'  # Paramètre à scanner
param = nsteps_values


ecrire_valeur("B1",0)
ecrire_valeur("kappa",0)
ecrire_valeur("theta0",0)
ecrire_valeur("thetadot0",0)


actualise_valeur()

outputs = []  # Liste pour stocker les fichiers de sortie
errors = []  # Liste pour stocker les erreurs

omega_0 = np.sqrt(12*B0*mu/(m*L*L))

T0 = 2 * np.pi / Omega  # Période théorique

tfin = Nperiod * T0

A = np.sqrt(theta0**2 + m * L * L * thetadot0**2 / (12 * B0 * mu))

phi = np.pi/2

if (thetadot0 != 0):
	phi = np.arctan(-theta0/thetadot0*omega_0)
"""

# Paramètres
executable = 'Exercice6_2025_student'  # Nom de l'exécutable compilé
repertoire = r"home/chatelin/Desktop/MyFiles/1quantique_physNum"
os.chdir(repertoire)

CONFIG_FILE = os.path.join(repertoire, "configuration.in")

# Liste des paramètres physiques par défaut
params_physiques = {
    "xL": -1.0,
    "xR": 1.0,
    "V0": 0.0,
    "xa": 0.0,
    "xb": 0.0,
    "omega0": 100.0,
    "x0": -0.5,
    "sigma_norm": 0.04,
    "n": 16,
    "nx": 512,
    "nsteps": 800,
    "tfin": 0.08
}

def ecrire_configuration(valeurs):
    with open(CONFIG_FILE, "w", encoding="utf-8") as f:
        for cle, valeur in valeurs.items():
            f.write(f"{cle} = {valeur}\n")

def lancer_simulation(output_file):
    cmd = f"./{executable} {CONFIG_FILE} output={output_file}"
    subprocess.run(cmd, shell=True)

def charger_resultats(nom_fichier):
    data = np.loadtxt(nom_fichier)
    return data[:, 0], data[:, 1:],  # t, [psi_re, psi_im, observables...]

def comparer_classique(x0, p0, omega0, t):
    xclass = x0 * np.cos(omega0 * t) + p0 / omega0 * np.sin(omega0 * t)
    pclass = p0 * np.cos(omega0 * t) - omega0 * x0 * np.sin(omega0 * t)
    return xclass, pclass

def tracer_observable(t, observable, label, question):
    plt.figure()
    plt.plot(t, observable, label=label)
    plt.xlabel("Temps t")
    plt.ylabel(label)
    plt.grid()
    plt.legend()
    plt.title(f"Q{question} : Evolution de {label} dans le temps")
    plt.show()

def heatmap_psi(t, x, psi_abs, question):
    plt.figure()
    plt.imshow(psi_abs, extent=[x[0], x[-1], t[-1], t[0]], aspect='auto', cmap='viridis')
    plt.colorbar(label='|ψ(x,t)|')
    plt.xlabel("x")
    plt.ylabel("Temps t")
    plt.title(f"Q{question} : Heatmap de |ψ(x,t)|")
    plt.show()

def etude_convergence(param_name, valeurs_testees, fixe_params, observable_fn, question):
    valeurs = []
    observables = []
    for v in valeurs_testees:
        params = fixe_params.copy()
        params[param_name] = v
        ecrire_configuration(params)
        output_file = f"test_{param_name}_{v}.out"
        lancer_simulation(output_file)
        t, data = charger_resultats(output_file)
        obs = observable_fn(t, data)
        valeurs.append(v)
        observables.append(obs)
    plt.figure()
    plt.plot(valeurs, observables, marker='o')
    plt.xlabel(param_name)
    plt.ylabel("Observable")
    plt.grid()
    plt.title(f"Q{question} : Convergence par rapport à {param_name}")
    plt.show()

def calcul_transmission(t, psi, x_vals):
    idx_0 = np.argmin(np.abs(x_vals))
    P_right = np.sum(np.abs(psi[:, idx_0:])**2, axis=1)
    P_left = np.sum(np.abs(psi[:, :idx_0])**2, axis=1)
    return P_left, P_right

def tracer_transmission(t, P_left, P_right, question):
    plt.figure()
    plt.plot(t, P_left, label="P(x<0)")
    plt.plot(t, P_right, label="P(x>0)")
    plt.xlabel("Temps t")
    plt.ylabel("Probabilité")
    plt.title(f"Q{question} : Transmission et réflexion")
    plt.grid()
    plt.legend()
    plt.show()

# Exemple d'utilisation pour la partie 6.3 (oscillateur harmonique)
ecrire_configuration(params_physiques)
lancer_simulation("oscillateur.out")
t, data = charger_resultats("oscillateur.out")

# Extraction des observables et de la fonction d'onde
psi_re = data[:, 0::10]  # exemple: colonne 0, 10, 20,... (adaptez selon la sortie)
psi_im = data[:, 1::10]  # exemple: colonne 1, 11, 21,... (adaptez)
psi_abs = np.sqrt(psi_re**2 + psi_im**2)

x_vals = np.linspace(params_physiques["xL"], params_physiques["xR"], params_physiques["nx"])

# Graphe Q6.3.i : heatmap de |ψ(x,t)|
heatmap_psi(t, x_vals, psi_abs, question="6.3.i")

# à compléter avec les autres graphiques et observables comme :
# tracer_observable(t, x_mean, "<x>", question="6.3.i")
# tracer_observable(t, p_mean, "<p>", question="6.3.i")
# tracer_transmission(t, P_left, P_right, question="6.4.ii")
