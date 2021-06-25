import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from coherence_gain_not_normalized import calc_g_av

def main():
    if len(sys.argv) != 2:
        print("Wrong number of parameters! Please provide only the temperature T [K] as an argument.")
        sys.exit()

    T = float(sys.argv[1])

    t = 1
    tau_step = 0.0001

    taus = list(np.arange(0, 0.01, tau_step)) # List of all "tau"s
    results_first_part = []
    results_second_part = []
    results_third_part = []
    results_final = []

    for tau in taus:
        print("t: ", format(t, '.5f'), " || tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'), end="\r")
        first_part, second_part, third_part, final_result = calc_g_av(t, tau, T)
        results_first_part.append(first_part)
        results_second_part.append(second_part)
        results_third_part.append(third_part)
        results_final.append(final_result)
    print("t: ", format(t, '.5f'), " || tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'))

    plt.plot(taus, results_first_part, ".-", label=r'$|p_+g_+|$')
    plt.plot(taus, results_second_part, ".-", label=r'$|p_-g_-|$')
    plt.plot(taus, results_third_part, ".-", label=r'$|D(t)|$')
    plt.plot(taus, results_final, ".-", label=r'$g_{av}$')
    plt.title(r'$t\ =\ ' + str(t) + '$')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()

    # Calculate the points half between |p+g+| and |p_g_|
    results_middle = []
    for i in range(len(results_first_part)):
        larger = max(results_first_part[i], results_second_part[i])
        smaller = min(results_first_part[i], results_second_part[i])
        half_distance = (larger - smaller)/2
        middle = smaller + half_distance
        results_middle.append(middle)

    plt.plot(taus, results_middle, ".-", label=r'$ half\ distance\ between\ |p_+g_+| \ and\ |p_-g_-| $')
    plt.legend()
    filename = 'T=' + str(T) + '_t=' + str(t) + '_detailed.pdf'
    plt.savefig(filename)
    # plt.show()

    # print("results_first_part: ", results_first_part)
    # print("results_second_part: ", results_second_part)
    # print("results_third_part: ", results_third_part)
    # print("results_final: ", results_final)
    # print("results_middle: ", results_middle)

    filename = 'data_T=' + str(T) + '_t=' + str(t) + '_detailed.csv'
    with open(filename, 'w') as file:
        save_list('|p+g+|', results_first_part, file)
        save_list('|p-g-|', results_second_part, file)
        save_list('|D(t)|', results_third_part, file)
        save_list('g_av', results_final, file)
        save_list('middle', results_middle, file)


def save_list(row_name, list, file):
    file.write(row_name + ',')
    for entry in list:
        file.write(str(entry) + ",")
    file.write('\n')

if __name__ == "__main__":
    main()
