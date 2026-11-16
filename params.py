# -*- coding: utf-8 -*-
"""
Paramètres pour le code ballistique
"""
import numpy as np 

params = {
    "v0": 2.0,                     # vitesse initiale (m/s)
    "h0": 10.0,                     # hauteur initiale (m)
    "alpha": 30,                    # angle en deg
    "g": 9.81,                      # gravité (m/s²)
    "t_end": 5.0,                   # durée maximale de simulation
    "nbr_points": 300,              # nb de points pour discrétisation
    "t": np.linspace(0, 5.0, 300),  # vecteur temporel
    "alpha_list": np.linspace(20,60,11), # liste d'angles en degré pour test
    "alpha_2": 20,                      # alpha partie 2 
    "mass" : 0.06,                       # mass
    "area":0.036,                        # surface
    "rho":1.3,                          # masse volumique de l'air
    "cd_list":[0, 0.01, 0.01, 0.01],            # plusieurs coefficient drag
    "cl_list":[0, 0, 0.18, 0.18],               # plusieurs coefficient lift
    "T0_list":[0, 0, 0, 0.01 * 9.81],       # plusieurs coefficient propulsion
    "T0": 0 * 0.01 * 9.81,
    "cd" : 0.18,
    "cl" : 0.01
}
