import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

### This script is used to create the plots of the simulation part. More precisely, it displays :
#       * The Global Diagnoses of the simulation
#       * The complexity of the scheme
#       * The comparison between the solution and the numerical solution
#
# Note that this script is self sufficient : it runs the simulation by itself and plots the results afterward.
# It reproduces all the plots of the simulation part of the report.  



def gaussian_solution_periodic(x, t, U_tot, L, c, sigma):
    x_mod = np.mod(x - c*t, L)
    if x_mod > L/2:
        x_mod -= L
    u = U_tot * np.exp(-(x_mod**2) / (sigma**2))
    return u

def wave_packet_periodic(x, t, U_tot, L, c, sigma):
    x_mod = np.mod(x - c*t, L)
    if x_mod > L/2:
        x_mod -= L
    u = U_tot * np.exp(-(x_mod**2) / (sigma**2)) * np.cos(2 * np.pi * 12 * x_mod / L)
    return u

def get_param(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        param_dict = {line.split('=')[0].strip(): eval(line.split('=')[1].strip()) for line in lines if not ('type' in line)}
        param_dict.update({line.split('=')[0].strip(): line.split('=')[1].strip() for line in lines if 'type' in line})
    return param_dict

def load_files(file_patterns, SIZE = ["32", "64", "128"]):
    return {pattern: [np.loadtxt(f"{pattern}_{size}.txt", delimiter=' ') for size in SIZE] for pattern in file_patterns}

def plot_diagnostics(params, datasets, linestyles):

    colors = ['tab:red', 'tab:blue', 'tab:green']
    labels = [r'$I_{h}^{n}$', r'$E_{h}^{n}$', r'$R_{h}^{n}$']
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(15, 5))

    for i, (param_key, param) in enumerate(params.items()):
        c, L, dx = param["c"], param["L"], param["dx"]
        sigma= L/16
        T, U, I, E, R = datasets["data"][i][:, 0], datasets["data"][i][:, 1:], datasets["diagnose"][i][:, 0], datasets["diagnose"][i][:, 1], datasets["diagnose"][i][:, 2]
        t = c*T/L

        ax1.plot(t[t<=1], I[t<=1], color=colors[0], linestyle=linestyles[i], label=r"$\dfrac{h}{\sigma} = $" + f"{dx / sigma:.3f}")
        ax1.grid()
        ax1.set_title(labels[0])
        ax1.set_xlabel(r'$\dfrac{ct}{L}$')
        ax1.set_ylabel("y")
        ax1.legend()

        ax2.plot(t[t<=1], E[t<=1], color=colors[1], linestyle=linestyles[i], label=r"$\dfrac{h}{\sigma} = $" + f"{dx / sigma:.3f}")
        ax2.grid()
        ax2.set_xlabel(r'$\dfrac{ct}{L}$')
        ax2.set_title(labels[1])
        ax2.legend()


        ax3.plot(t[t<=1], R[t<=1], color=colors[2], linestyle=linestyles[i], label=r"$\dfrac{h}{\sigma} = $" + f"{dx / sigma:.3f}")
        ax3.grid()
        ax3.set_xlabel(r'$\dfrac{ct}{L}$')
        ax3.set_title(labels[2])
        ax3.legend()

    plt.tight_layout()
    plt.show()


def plot_comlextity(params, datasets, alpha, intercept):
    global_error = []
    h_sigma = []
    for i, (param_key, param) in enumerate(params.items()):
        c, L, dx = param["c"], param["L"], param["dx"]
        sigma= L/16
        T, R = datasets["data"][i][:, 0], datasets["diagnose"][i][:, 2]
        t = c*T/L

        global_error.append(R[t==1/2][0])
        h_sigma.append(dx/sigma)

    global_error = np.array(global_error)
    h_sigma = np.array(h_sigma)
    fig, ax = plt.subplots(1,1)
    ax.loglog(h_sigma[global_error!=0], global_error[global_error!=0], "k-.")
    ax.loglog(h_sigma[global_error!=0], global_error[global_error!=0], "k^")
    ax.loglog(h_sigma[global_error!=0], (intercept*h_sigma[global_error!=0])**alpha, color= "blue", linestyle = ":"      , label = r'$\mathcal{O}\left(\dfrac{h}{\sigma}^{%.1d}\right)$'%alpha)
    ax.set_xlabel(r'$\dfrac{h}{\sigma}$')
    ax.set_ylabel(r'$R_{h}^{n}$')
    ax.grid()
    ax.legend()
    plt.show()


def plot_comparison(params, datasets, initial_type = 'GAUSSIAN'):


    labels = [r'$I_{h}^{n}$', r'$E_{h}^{n}$', r'$R_{h}^{n}$']

    for i, (param_key, param) in enumerate(params.items()):
        a, c, L, dx, u, N = param["a"], param["c"], param["L"], param["dx"], param["U"], param["N"]
        sigma= L/16
        T, U = datasets["data"][i][:, 0], datasets["data"][i][:, 1:]
    
        grid = datasets["grid"][i]
        t = c*T/L

        if param['Grid type'] == 'NON_UNIF':
            U = U / (1 - a * np.cos(2*np.pi * grid/L))
            grid = grid - a* L/(2*np.pi) * np.sin(2*np.pi*grid/L)
        
        xsi = np.linspace(-1/2, 1/2, 4*N)
        if initial_type == 'GAUSSIAN':
            U_sol = np.array([[gaussian_solution_periodic(x, time, u, L, c, sigma) for x in xsi] for time in T])
        else :
            U_sol = np.array([[wave_packet_periodic(x, time, u, L, c, sigma) for x in xsi] for time in T])
        

        U1 =  np.array([*U[t==1/4][0], *U[t==1/4][0]])
        U_sol1 = np.array([*U_sol[t==1/4][0], *U_sol[t==1/4][0]])
        U2 = np.array([*U[t==1/2][0], *U[t==1/2][0]])
        U_sol2 = np.array([*U_sol[t==1/2][0], *U_sol[t==1/2][0]])

        grid = np.array([*grid, *(grid + L)])



        fig, ax  = plt.subplots(1, 1)

        ax.plot(np.linspace(-1/2, 3/2, 8*N), U_sol1/u, 'k', label = 'solution')
        ax.plot(np.linspace(-1/2, 3/2, 8*N), U_sol2/u, 'k')
        ax.plot(grid, U2/u, 'tab:orange', label = r'$\dfrac{ct}{L} = 1/2$')#, linestyle='--')
        ax.plot(grid, U1/u, 'tab:red', label = r'$\dfrac{ct}{L} = 1/4$')#,  linestyle='--')
        ax.set_xlabel(r'$\dfrac{x}{L}$')
        ax.set_ylabel(r'$\dfrac{u(x, t)}{U}$')
        ax.grid()

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper right', fancybox=True,framealpha=0.3, facecolor='white', fontsize = 'large', edgecolor='black')

        plt.tight_layout()
        # plt.savefig(f"./figures/{params[param_key]['scheme_type']}/data_{dx/sigma:.3f}_{params[param_key]['scheme_type']}.pdf")
        # plt.savefig(f"./figures/NON_UNIF/data_{dx/sigma:.3f}_{params[param_key]['scheme_type']}.pdf")
        plt.show()

    return  


def plot_wave_packet(grid, U, xsi_2L, U_sol1, U_sol2, u, L, T, t = [1/4, 1/2]):
    """Plot the wave packet and save the figure."""
    fig, ax = plt.subplots(1, 1)

    U1 = np.array([*U[c*T/L == t[0]][0], *U[c*T/L == t[0]][0]])
    U2 = np.array([*U[c*T/L == t[1]][0], *U[c*T/L == t[1]][0]])
    grid = np.array([*grid, *(grid + L)])
    ax.plot(xsi_2L/L, U_sol1/u, 'forestgreen', linestyle='--', linewidth=1.2, label='Analytical solution')
    ax.plot(xsi_2L/L, U_sol2/u, 'forestgreen', linestyle='--', linewidth=1.2)
    ax.plot(grid/L, U2/u, 'k', label=r'$\dfrac{ct}{L} = %s$'%t[1])
    ax.plot(grid/L, U1/u, 'dimgrey', label=r'$\dfrac{ct}{L} = %s$'%t[0])
    ax.set_xlabel(r'$\dfrac{x}{L}$')
    ax.set_ylabel(r'$\dfrac{u_{i}^{n}}{U}$')
    ax.grid()
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":


    # UNIFORM GRID
    alpha = [4, 8, 8, 6]
    INITIAL_TYPE = 'GAUSSIAN' 
    for i, scheme in enumerate(['E2', 'E4', 'I4', 'ED']):

        # os.chdir("..")
        os.system("make clean")
        os.system("make")
        os.system(f"make run N=32 SCHEME_TYPE={scheme} GRID_TYPE=UNIF  INITIAL_TYPE={INITIAL_TYPE}")
        os.system(f"make run N=64 SCHEME_TYPE={scheme} GRID_TYPE=UNIF  INITIAL_TYPE={INITIAL_TYPE}")
        os.system(f"make run N=128 SCHEME_TYPE={scheme} GRID_TYPE=UNIF INITIAL_TYPE={INITIAL_TYPE}")
        os.chdir("./post_processing/data")

        file_patterns = ["data", "diagnose", "grid"]
        datasets = load_files(file_patterns)
        param_files = ["param_32.txt", "param_64.txt", "param_128.txt"]
        params = {f"param_{size}": get_param(f"param_{size}.txt") for size in ["32", "64", "128"]}

        linestyles = ["-", "--", "-."]

        intercept = 1
        plot_diagnostics(params, datasets, linestyles)
        plot_comlextity(params, datasets, alpha[i], intercept)
        plot_comparison(params, datasets)
   
        os.chdir("../..")
        os.system("make clean")
        os.chdir("./post_processing")



    # NON UNIFORM GRID - GAUSSIAN

    INITIAL_TYPE = 'GAUSSIAN'
    GRID_TYPE = 'NON_UNIF'
    time_sol = [1/2, 1]

    for i, scheme in enumerate(['E2', 'E4', 'I4', 'ED']):

        # os.chdir("..")
        os.system("make clean")
        os.system("make")
        os.system(f"make run N=128 SCHEME_TYPE={scheme} GRID_TYPE={GRID_TYPE} INITIAL_TYPE={INITIAL_TYPE}")
        os.chdir("./post_processing/data")

        file_patterns = ["data", "diagnose", "grid"]
        datasets = load_files(file_patterns, SIZE = ["128"])
        param_files = ["param_128.txt"]
        params = {f"param_{size}": get_param(f"param_{size}.txt") for size in ["128"]}
        
        L = float(params['param_128']['L'])
        c = float(params['param_128']['c'])
        u = float(params['param_128']['U'])
        a = float(params['param_128']['a'])
        sigma = L/16

        xsi_2L = np.linspace(-L/2, 3*L/2, 400)  # For smooth analytical plot
        if INITIAL_TYPE == 'GAUSSIAN':
            U_sol1 = np.array([gaussian_solution_periodic(x, L*time_sol[0]/(c), u, L, c, sigma) for x in xsi_2L])
            U_sol2 = np.array([gaussian_solution_periodic(x, L*time_sol[1]/(c), u, L, c, sigma) for x in xsi_2L])
        else :
            U_sol1 = np.array([wave_packet_periodic(x, L*time_sol[0]/(c), u, L, c, sigma) for x in xsi_2L])
            U_sol2 = np.array([wave_packet_periodic(x, L*time_sol[1]/(c), u, L, c, sigma) for x in xsi_2L])
   
        grid = datasets["grid"][0]
        time = datasets["data"][0][:, 0]
        U = datasets["data"][0][:, 1:]

        scheme_type = params['param_128']['Scheme type']

        plot_wave_packet(grid, U, xsi_2L, np.empty_like(U_sol1), np.empty_like(U_sol2), u, L, time, t = time_sol)

        if params['param_128']['Grid type'] == 'NON_UNIF':
            U = U / (1 - a * np.cos(2*np.pi * grid/L))
            grid = grid - a* L/(2*np.pi) * np.sin(2*np.pi*grid/L)
       
        plot_wave_packet(grid, U, xsi_2L, U_sol1, U_sol2, u, L, time, t = time_sol)


        os.chdir("../..")
        os.system("make clean")
        os.chdir("./post_processing")

    # NON UNIFORM GRID - WAVEPACKET
    
    INITIAL_TYPE = 'WAVEPACKET'
    GRID_TYPE = 'NON_UNIF'
    time_sol = [1/4, 1/2]

    for i, scheme in enumerate(['E2', 'E4', 'I4', 'ED']):

        # os.chdir("..")
        os.system("make clean")
        os.system("make")
        os.system(f"make run N=128 SCHEME_TYPE={scheme} GRID_TYPE={GRID_TYPE} INITIAL_TYPE={INITIAL_TYPE}")
        os.chdir("./post_processing/data")

        file_patterns = ["data", "diagnose", "grid"]
        datasets = load_files(file_patterns, SIZE = ["128"])
        param_files = ["param_128.txt"]
        params = {f"param_{size}": get_param(f"param_{size}.txt") for size in ["128"]}
        
        L = float(params['param_128']['L'])
        c = float(params['param_128']['c'])
        u = float(params['param_128']['U'])
        a = float(params['param_128']['a'])
        sigma = L/16

        xsi_2L = np.linspace(-L/2, 3*L/2, 400)  # For smooth analytical plot
        if INITIAL_TYPE == 'GAUSSIAN':
            U_sol1 = np.array([gaussian_solution_periodic(x, L*time_sol[0]/(c), u, L, c, sigma) for x in xsi_2L])
            U_sol2 = np.array([gaussian_solution_periodic(x, L*time_sol[1]/(c), u, L, c, sigma) for x in xsi_2L])
        else :
            U_sol1 = np.array([wave_packet_periodic(x, L*time_sol[0]/(c), u, L, c, sigma) for x in xsi_2L])
            U_sol2 = np.array([wave_packet_periodic(x, L*time_sol[1]/(c), u, L, c, sigma) for x in xsi_2L])
   
        grid = datasets["grid"][0]
        time = datasets["data"][0][:, 0]
        U = datasets["data"][0][:, 1:]

        scheme_type = params['param_128']['Scheme type']

        if params['param_128']['Grid type'] == 'NON_UNIF':
            U = U / (1 - a * np.cos(2*np.pi * grid/L))
            grid = grid - a* L/(2*np.pi) * np.sin(2*np.pi*grid/L)
       
        plot_wave_packet(grid, U, xsi_2L, U_sol1, U_sol2, u, L, time, t = time_sol)



        os.chdir("../..")
        os.system("make clean")
        os.chdir("./post_processing")


    
