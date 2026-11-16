# -*- coding: utf-8 -*-
"""
L3 ME

Projet "balistique"

@author: FERAY.M

"""
#import 


import numpy as np 
import matplotlib.pyplot as plt 
from ballistic.colored_messages import *
from ballistic.constantes import *
import scipy.optimize as sc
from scipy.integrate import odeint
from scipy.integrate import solve_ivp 


class AnalyticSolution:

    def __init__(self,params): # params est un dict
        self.v0 = params['v0']
        self.h0 = params['h0']
        self.alpha_deg = params['alpha']
        self.alpha_rad = np.deg2rad(self.alpha_deg)
        self.g = params['g']
        self.t = params['t']
        self.t_end = params['t_end']
        self.nbr_points = params['nbr_points']
        self.alpha_list = params['alpha_list']

        self.v_x = None
        self.v_z = None 
        self.v_norm = None 
        self.impact_values = None 
        self.pmax = None 
        self.alt_max = None 
        self.t_alt_max = None 

    #--------- Run ---------#
    def run(self):
        """
        run des méthodes
        """
        self.set_velocity()
        self.set_trajectory()  
        self.set_impact()
        self.set_alpha_opt()
        self.set_alt_max()
        self.set_pmax(self.alpha_opt) # alpha_opt est déja calculé a set_alpha_opt
        self.print_report()


    #--------- Set ---------#

    def set_alpha_opt(self):
        """
        find alpha opt
        """
        alpha_min_rad,alpha_max_rad = np.deg2rad(self.alpha_list[0]),np.deg2rad(self.alpha_list[-1])

        alphas = np.linspace(alpha_min_rad,alpha_max_rad,len(self.alpha_list)) # balaye les alpha 
        liste_point = [np.abs(self.set_f(alpha)) for alpha in alphas]
        index_min = np.argmin(liste_point)
        f_min = liste_point[index_min]
        alpha_opt = alphas[index_min]

        self.alpha_opt = alpha_opt
        self.liste_point = liste_point

    def set_alt_max(self):
        """
        set alt max
        """
        idx = np.argmin(np.abs(self.v_z))
        self.alt_max = self.z_t[idx]
        self.t_alt_max = self.t[idx]

    def set_alt_max_alpha(self):
        self.alt_max_alpha = []
        for alpha_i in self.alpha_list:
             alpha_i_rad = np.deg2rad(alpha_i)
             z_i = -0.5 * self.g * pow(self.t,2) + self.v0 * np.sin(alpha_i_rad) * self.t + self.h0
             v_z_i = self.v0 * np.sin(alpha_i_rad) - self.g * self.t
             idx = np.argmin(np.abs(v_z_i))
             self.alt_max_alpha.append((alpha_i,z_i[idx]))


    def set_f(self,alpha:float):
        """ la relation 1 """
        L = (2*self.g*self.h0)  / (self.v0**2)
        return np.cos(2 * alpha) + (np.sin(alpha) * (np.cos(2 * alpha) - L))  / (np.sqrt(np.sin(alpha)**2 + L))
    
    def set_impact(self):
        """
        set impact
        """
        t_i = pow(self.g,-1) * (self.v0 * np.sin(self.alpha_rad) + np.sqrt(pow(self.v0 * np.sin(self.alpha_rad),2) + 2 * self.g * self.h0))

        vx_i,vz_i= self.v0 * np.cos(self.alpha_rad), self.v0 * np.sin(self.alpha_rad) - self.g * t_i
        v_norm_i = np.sqrt(vx_i**2 + vz_i**2)

        x_i = self.v0 * np.cos(self.alpha_rad) * t_i # x impact

        théta_i = np.arctan(vz_i * pow(vx_i,-1))
        self.théta_i_deg = np.rad2deg(théta_i)

        self.impact_values = {"t_i": t_i, "x_i": x_i, "angle": self.théta_i_deg, "v": [vx_i, vz_i, v_norm_i]}

    def set_pmax(self,alpha_opt_rad):
        """
        set portée max
        """
        L = (2 * self.g * self.h0) / pow(self.v0,2)
        self.pmax = (pow(self.v0,2)/self.g) * ( 0.5 * np.sin(2*alpha_opt_rad) + np.cos(alpha_opt_rad) * np.sqrt(L + pow(np.sin(alpha_opt_rad),2)) )

    def set_trajectory(self):
        """
        set trajectory 
        """
        self.x_t = self.v0 * np.cos(self.alpha_rad) * self.t
        self.z_t = -0.5 * self.g * pow(self.t,2) + self.v0 * np.sin(self.alpha_rad) * self.t + self.h0
        return self.x_t,self.z_t
    
    def set_velocity(self):
        """
        set velocity 
        """
        self.v_x = self.v0 * np.cos(self.alpha_rad) * np.ones_like(self.t)
        self.v_z = self.v0 * np.sin(self.alpha_rad) - self.g * self.t
        self.v_norm = np.sqrt(self.v_x**2 + self.v_z**2)
        


    #--------- Print ---------#


    def print_report(self):

        set_title("Partie : Modèle analytique")

        set_section("Paramètres d'entrée")
        print("v₀ (m/s)             : %.2f" % self.v0)
        print("h₀ (m)               : %.2f" % self.h0)
        print("alpha (°)            : %.2f" % self.alpha_deg)
        print("g (m/s²)             : %.2f" % self.g)
        print("t_end (s)            : %.2f" % self.t_end)
        print("points de calcul     : %d" % self.nbr_points)
        print("liste alpha          :", self.alpha_list)

        set_line("=", 60, "light_blue")


        set_section("Optimal balistique")

        print("α optimal (°)        : %.2f" % np.rad2deg(self.alpha_opt))
        print("portée maximale (m)  : %.4f" % self.pmax)
        print("hauteur max (m)      : %.2f" % self.alt_max)
        print("temps sommet (s)     : %.3f" % self.t_alt_max)

        set_section("Valeurs à l'impact")

        print("temps d'impact (s)   : %.3f" % self.impact_values["t_i"])
        print("x impact (m)         : %.4f" % self.impact_values["x_i"])
        print("vz impact (m/s)      : %.2f" % self.impact_values["v"][1])
        print("vx impact (m/s)      : %.2f" % self.impact_values["v"][0])
        print("|v| impact (m/s)     : %.2f" % self.impact_values["v"][2])
        print("angle impact (°)     : %.2f" % self.théta_i_deg)

    #--------- Plot ---------#

    def plot_all(self):
        """
        Lance les plots
        """
        self.plot_components()
        self.plot_pmax()
        self.plot_trajectories()
        self.plot_itérations()

    def plot_components(self):
        """
        plot components
        """

        fig, (ax1, ax2) = plt.subplots(2, 1) 
        fig.suptitle('composantes',fontsize=12.)  
        ax1.plot(self.t,self.v_z,label='vz(t)')
        ax1.grid()
        ax1.set_xlabel("t (s)",fontsize=12.)
        ax1.set_ylabel("vz(t) (m/s)",fontsize=12.)
        ax1.fill_between(self.t, self.v_z,color='green', alpha=0.5)

        ax2.plot(self.t,self.v_norm,label='norme v(t)')
        ax2.grid()
        ax2.set_xlabel("t (s)",fontsize=12.)
        ax2.set_ylabel("v(t) (m/s)",fontsize=12.)
        ax2.fill_between(self.t, self.v_norm,color='red', alpha=0.5)

       
        plt.tight_layout()
        savefig("plot_components.png")
        plt.show()

    def plot_itérations(self):
        """
        plot recherche de alpha_opt
        """

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(self.alpha_list,self.liste_point,lw=3)
        ax.set_ylabel(r'Fonction $f(\alpha)$')
        ax.set_xlabel(r'$\alpha$ (°)')
        ax.set_title(r"Recherche de l'angle optimal via $f(\alpha)$")       
        ax.grid(True, which="major", linestyle="-", alpha=0.5)

        savefig("itération.png")
        plt.show()

    def plot_pmax(self):
        """
        plot porté max en foncton de alpha
        """
        fig, ax = plt.subplots(figsize=(8, 4))

        pmax_values = []
        alphas_rad = np.deg2rad(self.alpha_list)

        for alpha_rad in alphas_rad:
            L = (2*self.g*self.h0)  / (self.v0**2)
            p = (self.v0**2 / self.g) * (0.5 * np.sin(2 * alpha_rad)+ np.cos(alpha_rad) * np.sqrt(L + np.sin(alpha_rad)**2))
            pmax_values.append(p)

        # Tracé global
        ax.plot(self.alpha_list, pmax_values, '-o', color='tab:blue', linewidth=3)

        ax.set_title(r'Portée maximale en fonction de $\alpha$', fontsize=14, fontweight="bold")
        ax.set_xlabel('alpha (°)', fontsize=12)
        ax.set_ylabel('p_max (m)', fontsize=12)

        ax.grid(True, which="major", linestyle="-", alpha=0.5)



        plt.tight_layout()
        savefig("pmax_ana.png")
        plt.show()

    def plot_trajectories(self):

        """
        Trace toutes les trajectoires en fonctions des alphas dans alpha_liste.
        """
   
        fig, ax = plt.subplots(figsize=(8, 4))

        for alpha_i in self.alpha_list:
            x = self.v0 * np.cos(np.deg2rad(alpha_i)) * self.t
            z = -0.5 * self.g * pow(self.t,2) + self.v0 * np.sin(np.deg2rad(alpha_i)) * self.t + self.h0
            ax.plot(x, z, label=rf"$\alpha = {alpha_i}^\circ$")
            

        ax.set_title(r'Trajectoires en fonction de $\alpha$', fontsize=14, fontweight="bold")
        ax.set_xlabel('x[m]', fontsize=12)
        ax.set_ylabel('z[m]', fontsize=12)

        ax.grid(True, which="major", linestyle="-", alpha=0.5)
        ax.set_ylim(0, 40)
        ax.set_xlim(0, 60)
        ax.legend(title="Angles")


        plt.tight_layout()
        savefig("trajectoire_modèle_1.png")
        plt.show()


class NumericalSolution: 

    def __init__(self,params):
        self.v0 = params['v0']
        self.h = params['h0']
        self.alpha_deg = params['alpha_2']
        self.alpha_rad = np.deg2rad(self.alpha_deg)
        self.g = params['g']
        self.npt = params['nbr_points']
        self.cd = params['cd'] 
        self.cl = params['cl'] 
        self.mass = params['mass']
        self.rho = params['rho'] 
        self.S = params['area']
        self.T0 = params['T0']
        self.ct = self.T0/self.mass
        
        self.t = params['t']
        self.B = (self.rho * self.S) / (2*self.mass)

        params_ana = params.copy()
        params_ana['alpha'] = params['alpha_2']
        self.analytical_model = AnalyticSolution(params_ana)
        self.analytical_x, self.analytical_z = self.analytical_model.set_trajectory()

        #self.last_case = None # instruction 
        #self.x, self.z, self.v_x, self.v_z = None, None, None, None # modèle sans drag odeint 
        #self.x_d, self.z_d, self.v_x_d, self.v_z_d = None, None, None, None # modèle avec drag odeint 
        #self.x_no_drag_ivp, self.z_no_drag_ivp, self.v_x_no_drag_ivp, self.v_z_no_drag_ivp = None, None, None, None # modèle sans drag ivp
        #self.x_drag_ivp, self.z_drag_ivp, self.v_x_drag_ivp, self.v_z_drag_ivp = None, None, None, None # modèle avec drag ivp 

        #self.impact_without_drag = None  # dictionnaire impact sans drag 
        #self.impact_with_drag = None # dictionnaire impact avec drag 

        # on ne peut pas déclarer les self à cause de hasattr 

    #--------- Run ---------#

    def run_all(self, instruction):
        """
        instruction peut venir de PyQt5 (str) OU être déjà un SimulationMode.
        """
        self.last_case = instruction   

        if instruction == 'no_drag':
            self.solve_trajectory_without_drag()
            self.set_impact()
            self.set_validation()
            self.print_report()

        elif instruction == 'drag':
            self.solve_trajectory_with_drag()
            self.set_impact()
            self.print_report()

        elif instruction == 'no_drag_ivp':
            self.solve_trajectory_no_drag_ivp()
            self.set_impact()
            self.print_report()
            self.set_validation()

        elif instruction == 'drag_ivp':
            self.solve_trajectory_drag_ivp()
            self.set_impact()
            self.print_report()

        else:
            raise ValueError("Mode non géré.")

        self.plot_trajectory()

    #--------- ODE ---------#

    def ode_without_drag(self, y, t):
        """
        ODE to solve
        dy/dt = f(t, y)
        y = [x, z, vx, vz]
        """
        dy = np.zeros(4)
        dy[0] = y[2]  # dx / dt = v_x
        dy[1] = y[3]  # dz / dt = v_z
        dy[2] = 0     # dv_x / dt = 0
        dy[3] = -self.g    # dv_z / dt = - g
        return dy
    
    def ode_with_drag(self, y, t):
        """
        ODE to solve
        dy/dt = f(t, y)
        y = [x, z, vx, vz]
        cz =clcosθ−cdsinθ, cx =−cdcosθ−clsinθ # déballer les attributs pour plus de lisibilité pour la version finale 
        """
        theta = np.arctan2(y[3],y[2]) # arctan(vz/vx)
        cz = self.cl * np.cos(theta) - self.cd * np.sin(theta)
        cx = - self.cd * np.cos(theta) - self.cl * np.sin(theta)
        dy = np.zeros(4)
        dy[0] = y[2]  # dx / dt = v_x
        dy[1] = y[3]  # dz / dt = v_z
        v2 = y[2]**2 + y[3]**2
        dy[2] = self.B * v2 * cx + self.ct * np.cos(theta)
        dy[3] = -self.g + self.B * v2 * cz + self.ct * np.sin(theta)
        return dy
    
    def ode_without_drag_ivp(self,t, y): # on inverse
        """
        ODE to solve
        dy/dt = f(t, y)
        y = [x, z, vx, vz]
        """
        dy = np.zeros(4)
        dy[0] = y[2]  # dx / dt = v_x
        dy[1] = y[3]  # dz / dt = v_z
        dy[2] = 0     # dv_x / dt = 0
        dy[3] = -self.g    # dv_z / dt = - g
        return dy
    
    def ode_with_drag_ivp(self, t,y): # inversé
        """
        ODE to solve
        dy/dt = f(t, y)
        y = [x, z, vx, vz]
        cz =clcosθ−cdsinθ, cx =−cdcosθ−clsinθ # déballer les attributs pour plus de lisibilité pour la version finale 
        """
        theta = np.arctan2(y[3],y[2]) # arctan(vz/vx)
        cz = self.cl * np.cos(theta) - self.cd * np.sin(theta)
        cx = - self.cd * np.cos(theta) - self.cl * np.sin(theta)
        dy = np.zeros(4)
        dy[0] = y[2]  # dx / dt = v_x
        dy[1] = y[3]  # dz / dt = v_z
        v2 = y[2]**2 + y[3]**2
        dy[2] = self.B * v2 * cx + self.ct * np.cos(theta)
        dy[3] = -self.g + self.B * v2 * cz + self.ct * np.sin(theta)
        return dy
    


        #--------- Solve ---------#



    def solve_trajectory_without_drag(self):
        """
        Calcul la trajectoire odeint sans drag
        """
        # initial condition
        y_init = [0, self.h, self.v0 * np.cos(self.alpha_rad), self.v0 * np.sin(self.alpha_rad)]
        y = odeint(self.ode_without_drag, y_init, self.t, full_output=False)       # résolution de l'ode
        self.x, self.z, self.v_x, self.v_z = y[:, 0], y[:, 1], y[:, 2], y[:, 3]

    def solve_trajectory_with_drag(self):
        """
        Calcul la trajectoire odeint avec drag
        """
        y_init = [0, self.h, self.v0 * np.cos(self.alpha_rad), self.v0 * np.sin(self.alpha_rad)]
        y = odeint(self.ode_with_drag, y_init, self.t, full_output=False)       # résolution de l'ode
        self.x_d, self.z_d, self.v_x_d, self.v_z_d = y[:, 0], y[:, 1], y[:, 2], y[:, 3]


    def solve_trajectory_no_drag_ivp(self):
        """
        Calcul la trajectoire solve_ivp sans drag
        """
        y_init = [0, self.h, self.v0 * np.cos(self.alpha_rad), self.v0 * np.sin(self.alpha_rad)]
        t_span = (self.t[0], self.t[-1])
        t_eval = np.linspace(self.t[0],self.t[-1],100)
        #def ground(t,y):
            #return y[1]
        #ground.terminal = True
        #ground.direction = -1

        sol = solve_ivp(self.ode_without_drag_ivp, t_span, y_init, method='LSODA', t_eval=t_eval, dense_output=False, events=None, vectorized=False, args=None)
        self.t_ivp = sol.t
        self.x_no_drag_ivp, self.z_no_drag_ivp, self.v_x_no_drag_ivp, self.v_z_no_drag_ivp = sol.y[0], sol.y[1], sol.y[2], sol.y[3]    


    def solve_trajectory_drag_ivp(self):
        """
        Calcul la trajectoire solve_ivp avec drag
        """
        y_init = [0, self.h, self.v0 * np.cos(self.alpha_rad), self.v0 * np.sin(self.alpha_rad)]
        t_span = (self.t[0], self.t[-1])
        t_eval = np.linspace(self.t[0],self.t[-1],100)

        #def ground(t,y):
            #return y[1]
        #ground.terminal = True
        #ground.direction = -1

        sol = solve_ivp(self.ode_with_drag_ivp, t_span, y_init, method='LSODA', t_eval=t_eval, dense_output=False, events=None, vectorized=False, args=None)
        self.t_ivp = sol.t
        self.x_drag_ivp, self.z_drag_ivp, self.v_x_drag_ivp, self.v_z_drag_ivp = sol.y[0], sol.y[1], sol.y[2], sol.y[3]    
    

        #--------- Set ---------#

    def set_impact(self):
        """
        set l'impact selon le cas 
        """
        if hasattr(self,'x'):
            self.impact_without_drag = self.set_impact_values(self.t,self.x,self.z,self.v_x,self.v_z)
        if hasattr(self,'x_d'):
            self.impact_with_drag = self.set_impact_values(self.t,self.x_d,self.z_d,self.v_x_d,self.v_z_d)
        if hasattr(self,'x_no_drag_ivp'):
            self.impact_without_drag = self.set_impact_values(self.t_ivp,self.x_no_drag_ivp,self.z_no_drag_ivp,self.v_x_no_drag_ivp,self.v_z_no_drag_ivp)
        if hasattr(self,'x_drag_ivp'):
            self.impact_with_drag = self.set_impact_values(self.t_ivp,self.x_drag_ivp,self.z_drag_ivp,self.v_x_drag_ivp,self.v_z_drag_ivp)

    def set_impact_last_case(self):
        """
        set l'impact selon le cas (version last_case)
        """
        case = self.last_case
        if case == "no_drag":
            self.impact_without_drag = self.set_impact_values(self.t,self.x, self.z, self.v_x, self.v_z)
        elif case == "drag":
            self.impact_with_drag = self.set_impact_values(self.t,self.x_d, self.z_d, self.v_x_d, self.v_z_d)
        elif case == "no_drag_ivp":
            self.impact_without_drag = self.set_impact_values(self.t,self.x_no_drag_ivp, self.z_no_drag_ivp, self.v_x_no_drag_ivp, self.v_z_no_drag_ivp)
        elif case == "drag_ivp":
            self.impact_with_drag = self.set_impact_values(self.t,self.x_drag_ivp, self.z_drag_ivp, self.v_x_drag_ivp, self.v_z_drag_ivp)
            
    def set_impact_values(self,t,x,z,v_x,v_z): 
        """
        partie à modifier par les étudiants

        méthode: trouver le temps d'impact t_i tel que v_z est juste au dessus de 0.
        v_z(t_{i+1}) < 0. Puis remplir v_x, v_z, v, theta_i, x_i à cet instant d'impact.
        """
        def interpo(a, n, u):
            return u[n] + a * (u[n + 1] - u[n])

        n = 0

        while n < len(z) - 1  and z[n+1]>0: 
            n += 1 
        
        if n == len(z) - 1 and z[n]>0:
            raise ValueError("aucune valeur de z et nulle ou négative cad qu'on ne touche jamais le sol")


        a = - z[n]/(z[n+1] - z[n])

        t_i = interpo(a, n, t)
        x_i = interpo(a, n, x)
        v_x = interpo(a, n, v_x)
        v_z = interpo(a, n, v_z)
        v   = np.sqrt(v_x**2 + v_z**2)
        theta_i = np.arctan2(v_z, v_x)   
        
        return {"t_i": t_i, "x_i": x_i, "angle": np.rad2deg(theta_i), "v": [v_x, v_z, v]}

    def set_validation(self):
        """
        Lance la validation analytique
        """
        if hasattr(self,'x'):
            self.validation(self.t,self.x,self.z)
        if hasattr(self,'x_no_drag_ivp'):
            self.validation(self.t_ivp,self.x_no_drag_ivp,self.z_no_drag_ivp)

    def validation(self, t, x, z):
        """
        On va regarder si la méthode est valide par rapport à la solution de référence
        """
        epsilon = 1e-8

        set_section("Validation")
        set_msg(("analytical solution at t = %f" % t[-1]))

        x_ref = np.interp(t[-1], self.t, self.analytical_x)
        z_ref = np.interp(t[-1], self.t, self.analytical_z)
        # --------------------------

        set_msg(("x, z                       : %f  %f" % (x_ref, z_ref)))
        set_msg(("numerical solution at the same time:"))
        set_msg(("x, z                       : %f  %f" % (x[-1], z[-1])))

        err_x = abs(x_ref - x[-1])
        err_z = abs(z_ref - z[-1])
        set_msg(f"err_x = {err_x:.16e}")
        set_msg(f"err_z = {err_z:.16e}")

        if err_x < epsilon and err_z < epsilon:
            set_msg(("Validation OK : solution numérique proche de l’analytique."))
        else:
            raise ValueError("Validation échouée : erreur trop grande.")
        
        #--------- Print ---------#

    def _print_aero(self):
        set_info("Paramètres aérodynamiques", color="yellow")
        set_msg(f"  cd          : {self.cd:.3f}")
        set_msg(f"  cl          : {self.cl:.3f}")
        set_msg(f"  ρ (kg/m³)   : {self.rho:.3f}")
        set_msg(f"  S (m²)      : {self.S:.4f}")
        set_msg(f"  B=ρS/(2m)   : {self.B:.5f}")

        set_info("Paramètres masse / poussée", color="magenta")
        set_msg(f"  m (kg)      : {self.mass:.3f}")
        set_msg(f"  T₀ (N)      : {self.T0:.3f}")
        set_msg(f"  cₜ=T₀/m     : {self.ct:.3f} (m/s²)")


    def _print_common(self):
        set_info("Paramètres communs", color="cyan")
        set_msg(f"  v₀ (m/s) : {self.v0:.2f}")
        set_msg(f"  h₀ (m)   : {self.h:.2f}")
        set_msg(f"  α (°)    : {self.alpha_deg:.2f}")
        set_msg(f"  g (m/s²) : {self.g:.2f}")

        set_info("Paramètres temporels", color="light_blue")
        set_msg(f"  t début (s) : {self.t[0]:.4f}")
        set_msg(f"  t fin (s)   : {self.t[-1]:.4f}")
        set_msg(f"  nb points   : {len(self.t)}")


    def print_report(self):
        set_line("=", 60, "light_blue")
        set_title("Modèle numérique")

        case = getattr(self, "last_case", None)

        if case == "no_drag":
            set_section('Cas sans trainée modèle Odeint')
            self._print_common()
            self._print_impact(getattr(self, "impact_without_drag", None), "Résultats (sans traînée)")

        elif case == "drag":
            set_section('Cas avec trainée modèle Odeint')
            self._print_common()
            self._print_aero()
            self._print_impact(getattr(self, "impact_with_drag", None), "Résultats (avec traînée)")

        elif case == "no_drag_ivp":
            set_section('Cas sans trainée modèle solve ivp')
            self._print_common()
            self._print_impact(getattr(self, "impact_without_drag", None), "Résultats (ivp sans traînée)")

        elif case == "drag_ivp":
            set_section('Cas avec trainée modèle solve_ivp')
            self._print_common()
            self._print_aero()
            self._print_impact(getattr(self, "impact_with_drag", None), "Résultats (ivp avec traînée)")
        else:
            set_alert("Aucune trajectoire calculée (last_case manquant).")
    

    def _print_impact(self, impact_dict, label):
        if not impact_dict:
            return
        set_info(label, color="light_green")
        set_msg("IMPACT")
        set_msg(f"temps jusqu'à l'impact : {impact_dict['t_i']:.4f} s")
        set_msg(f"portée                 : {impact_dict['x_i']:.4f} m")
        set_msg(f"angle d'impact         : {impact_dict['angle']:.2f} °")
        set_msg(f"vitesse d'impact       : {impact_dict['v'][2]:.4f} m/s")


        #--------- Plot ---------#

    def plot_trajectory(self):

        fig, ax = plt.subplots(figsize=(8, 4))

        if hasattr(self, "x"):
            ax.plot(self.x,self.z,label=f'alpha = {self.alpha_deg}° pmax = {np.round( self.impact_without_drag["x_i"],2)} m')

        if hasattr(self, "x_d"):
            ax.plot(self.x_d,self.z_d,label=f'alpha = {self.alpha_deg}° pmax = { np.round(self.impact_with_drag["x_i"],2)} m')

        if hasattr(self,'x_drag_ivp'):
            ax.plot(self.x_drag_ivp, self.z_drag_ivp,label=f'alpha = {self.alpha_deg} ° pmax = {np.round(self.impact_with_drag["x_i"],2)} m')

        if hasattr(self,'x_no_drag_ivp'):
            ax.plot(self.x_no_drag_ivp, self.z_no_drag_ivp,label=f'alpha = {self.alpha_deg} ° pmax = {np.round(self.impact_without_drag["x_i"],2)} m')

        ax.set_title(r'Trajectoire ', fontsize=14, fontweight="bold")
        ax.set_xlabel('x[m]', fontsize=12)
        ax.set_ylabel('z[m]', fontsize=12)

        ax.grid(True, which="major", linestyle="-", alpha=0.5)
        ax.set_ylim(0, 30)
        ax.set_xlim(0, 70)
        ax.legend(title="Angles")


        plt.tight_layout()
        savefig("trajectoire_num.png")
        plt.show()

    def plot_trajectory_last_case(self):
        """
        Affiche la trajectoire (last_case)
        """
        fig, ax = plt.subplots(figsize=(8, 4))
        case = self.last_case

        if case == "no_drag":
            ax.plot(self.x, self.z, label=f"odeint sans traînée p={self.impact_without_drag['x_i']:.2f} m")

        elif case == "drag":
            ax.plot(self.x_d, self.z_d, label=f"odeint avec traînée p={self.impact_with_drag['x_i']:.2f} m")

        elif case == "no_drag_ivp":
            ax.plot(self.x_no_drag_ivp, self.z_no_drag_ivp, label=f"ivp sans traînée p={self.impact_without_drag['x_i']:.2f} m")

        elif case == "drag_ivp":
            ax.plot(self.x_drag_ivp, self.z_drag_ivp, label=f"ivp avec traînée p={self.impact_with_drag['x_i']:.2f} m")

        ax.grid(True)
        ax.legend()
        savefig("trajectoires_last_case.png")
        plt.show()


        #--------- Class méthode ---------#

    @classmethod
    def plot_trajectories(cls, base_params):
        """
        Affiche les trajectoires sur un même graphique
        """
        cd_list = base_params["cd_list"]
        cl_list = base_params["cl_list"]
        T_list  = base_params["T0_list"]      

        fig, ax = plt.subplots()

        for cd_i, cl_i, T_i in zip(cd_list, cl_list, T_list):
            params_i = base_params.copy()
            params_i['cd'] = cd_i
            params_i['cl'] = cl_i
            params_i['T0'] = T_i

            model_i = cls(params_i)
            model_i.solve_trajectory_with_drag()
            model_i.set_impact()      
            x_i = model_i.x_d
            z_i = model_i.z_d     
            label = rf"$c_d={cd_i:.2f},\ c_l={cl_i:.2f},\ T_0={T_i:.2f}\,\mathrm{{N}}$"
            ax.plot(x_i,z_i,lw=3,label=label)

        ax.set_title(r'Trajectoires modèle avec drag', fontsize=14, fontweight="bold")
        ax.set_xlabel('x[m]', fontsize=12)
        ax.set_ylabel('z[m]', fontsize=12)

        ax.grid(True, which="major", linestyle="-", alpha=0.5)
        ax.set_ylim(0, 35)
        ax.set_xlim(0, 90)
        ax.legend(fontsize=10, loc="upper right",title="Paramètres", frameon=True)


        plt.tight_layout()
        savefig("trajectoirs_modèle_num.png")
        plt.show()
        

    @classmethod
    def isocontours(cls,base_params):
        """
        Affiche les isocontours
        """
        alpha = np.linspace(20,60,11)
        cd = np.linspace(0,0.1,11)
        P = np.zeros((len(cd), len(alpha)))


        for i, cd_i in enumerate(cd): 
            for j, alpha_j in enumerate(alpha):
                params_i = base_params.copy()
                params_i['cd'] = cd_i
                params_i['cl'] = 0.18
                params_i['alpha_2'] = alpha_j
                params_i["T0"] = 0

                model_i = NumericalSolution(params_i)
                model_i.solve_trajectory_with_drag()
                model_i.set_impact()

                p = model_i.impact_with_drag['x_i']
                P[i, j] = p


        fig, ax = plt.subplots(figsize=(7, 5))

        A, C = np.meshgrid(alpha,cd)

        cont = ax.contourf(A, C, P, levels=20, cmap='jet')
        cbar = fig.colorbar(cont, ax=ax)
        cbar.set_label("Portée p (m)")

        ax.set_xlabel(r"$\alpha\ (\degree)$")
        ax.set_ylabel(r"$c_d$")
        ax.set_title(f"Isocontours de la portée en fonction de $\\alpha$ et $c_d$ T0 ={base_params['T0']:2f} N")
        ax.grid(True, which="major", linestyle="-", alpha=0.5)
        savefig("isocontours.png")
        plt.tight_layout()
        plt.show()




def savefig(name, folder="output", dpi=300):
    from pathlib import Path
    out_dir = Path(__file__).resolve().parent / folder
    out_dir.mkdir(exist_ok=True)
    plt.savefig(out_dir / name, dpi=dpi, bbox_inches="tight")