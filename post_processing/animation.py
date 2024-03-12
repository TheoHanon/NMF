import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 
import argparse

### This script is used to animate the results of the simulation.


def load_simulation_data(file_base, N):
    
    data = np.loadtxt(f"{file_base}/data_{N}.txt", delimiter=' ')
    diagnose = np.loadtxt(f"{file_base}/diagnose_{N}.txt", delimiter=' ')
    grid = np.loadtxt(f"{file_base}/grid_{N}.txt")
    with open(f"{file_base}/param_{N}.txt", 'r') as file:
        params = dict(line.strip().split('=') for line in file)
    return data, diagnose, grid, params


def animate(grid, U, T, c, L, u):

    fig, ax = plt.subplots()
    line, = ax.plot(grid/L, U[0]/u)
    ax.set_xlabel(r'$\dfrac{x}{L}$')
    ax.set_ylabel(r'$\dfrac{u(x, t)}{U}$')


    def update(frame):
        line.set_xdata(grid/L)
        line.set_ydata(U[frame]/u)
        ax.set_title(r'$\dfrac{ct}{L}$' + f'= {c*T[frame]/L:.3f}')
        return line,

    ani = animation.FuncAnimation(fig, update, frames=len(T), interval=10)
    plt.show()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-N", default=128)
    args = parser.parse_args()
    file_base = "./post_processing/data"  # Update this with the actual path
    N = args.N # Example, use the desired N value
    data, diagnose, grid, params = load_simulation_data(file_base, N)

    T, V = data[:, 0], data[:, 1:]
    
    L = float(params['L'])
    c = float(params['c'])
    u = float(params['U'])
    sigma = L/16

    animate(grid, V, T, c, L, u)












