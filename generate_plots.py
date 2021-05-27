import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
import os

def main():
    if len(sys.argv) != 2:
        print("Wrong number of parameters! Please provide only the temperature T [K] as an argument.")
        sys.exit()

    T = float(sys.argv[1])

    filename = 'data_T=' + str(T) + '.pkl'

    if os.path.exists(filename):
        with open(filename, 'rb') as file:
            data = pickle.load(file)
    else:
        print("Data for this temperature hasn't been generated yet! Run the generate_data.py script.")
        sys.exit()

    # print("data: ", data[:])
    taus = data[-1][1]
    # "data" is an array of triples. We want to take the 0th entry from each triple, because these
    # correspond to times "t"
    ts = [triple[0] for triple in data]
    # print("ts: ", ts)
    # print("taus: ", taus)
    X, Y = np.meshgrid(taus, ts)
    Z = [triple[2] for triple in data]

    # Fill the 'blanks' in the Z set
    target_length = len(Z[-1])
    for i in range(len(Z)):
        if len(Z[i]) < target_length:
            for j in range(target_length - len(Z[i])):
                Z[i].append(np.nan)
            # print("Z: ", Z[i])
    Z = np.array(Z)

    # print("X: ", X)
    # print("Y: ", Y)
    # print("Z: ", Z)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    scatter = ax.plot_wireframe(X, Y, Z)
    ax.set_title(r'$T\ =\ ' + str(T) + '$')
    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel(r'$t$')
    ax.set_zlabel(r'$g_{av}$')

    plt.show()

if __name__ == "__main__":
    main()
