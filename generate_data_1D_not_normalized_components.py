import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from coherence_gain_not_normalized_components import calc_g_av_components

def main():
    if len(sys.argv) != 2:
        print("Wrong number of parameters! Please provide only the temperature T [K] as an argument.")
        sys.exit()

    T = float(sys.argv[1])

    t = 1
    tau_step = 0.0001

    taus = list(np.arange(0, 0.01, tau_step)) # List of all "tau"s
    results_W_0_real = []
    results_W_1_real = []
    results_W_2_real = []
    results_W_3_real = []

    results_W_0_imag = []
    results_W_1_imag = []
    results_W_2_imag = []
    results_W_3_imag = []

    results_pure_phase_real = []
    results_pure_phase_imag = []

    for tau in taus:
        print("t: ", format(t, '.5f'), " || tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'), end="\r")
        W_0, phase_W_1, phase_W_2, W_3, pure_phase = calc_g_av_components(t, tau, T)
        results_W_0_real.append(W_0.real)
        results_W_0_imag.append(W_0.imag)
        results_W_1_real.append(phase_W_1.real)
        results_W_1_imag.append(phase_W_1.imag)
        results_W_2_real.append(phase_W_2.real)
        results_W_2_imag.append(phase_W_2.imag)
        results_W_3_real.append(W_3.real)
        results_W_3_imag.append(W_3.imag)

        results_pure_phase_real.append(pure_phase.real)
        results_pure_phase_imag.append(pure_phase.imag)

    print("t: ", format(t, '.5f'), " || tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'))


    plt.plot(taus, results_W_0_real, ".-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(taus, results_W_1_real, ".-", label=r'$e^{-iE\tau/\hbar}\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(taus, results_W_2_real, ".-", label=r'$\langle e^{iE\tau/\hbar}W(\tau)W(t-\tau)\rangle$')
    plt.plot(taus, results_W_3_real, ".-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(taus, results_pure_phase_real, ".-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$t\ =\ ' + str(t) + '\ (real)$')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    plt.legend()
    filename = 'T=' + str(T) + '_t=' + str(t) + '_components_real.pdf'
    plt.savefig(filename)
    # plt.show()

    print("results_W_0_real: ", results_W_0_real)
    print("results_W_1_real: ", results_W_1_real)
    print("results_W_2_real: ", results_W_2_real)
    print("results_W_3_real: ", results_W_3_real)
    print("results_pure_phase_real: ", results_pure_phase_real)

    filename = 'data_T=' + str(T) + '_t=' + str(t) + '_components_real.csv'
    with open(filename, 'w') as file:
        save_list('W_0_real', results_W_0_real, file)
        save_list('W_1_real', results_W_1_real, file)
        save_list('W_2_real', results_W_2_real, file)
        save_list('W_3_real', results_W_3_real, file)
        save_list('phase_real', results_pure_phase_real, file)


    plt.clf()
    plt.plot(taus, results_W_0_imag, ".-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(taus, results_W_1_imag, ".-", label=r'$e^{-iE\tau/\hbar}\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(taus, results_W_2_imag, ".-", label=r'$\langle e^{iE\tau/\hbar}W(\tau)W(t-\tau)\rangle$')
    plt.plot(taus, results_W_3_imag, ".-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(taus, results_pure_phase_imag, ".-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$t\ =\ ' + str(t) + '\ (imag)$')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    plt.legend()
    filename = 'T=' + str(T) + '_t=' + str(t) + '_components_imag.pdf'
    plt.savefig(filename)
    # plt.show()

    print("results_W_0_imag: ", results_W_0_imag)
    print("results_W_1_imag: ", results_W_1_imag)
    print("results_W_2_imag: ", results_W_2_imag)
    print("results_W_3_imag: ", results_W_3_imag)
    print("results_pure_phase_imag: ", results_pure_phase_imag)

    filename = 'data_T=' + str(T) + '_t=' + str(t) + '_components_imag.csv'
    with open(filename, 'w') as file:
        save_list('W_0_imag', results_W_0_imag, file)
        save_list('W_1_imag', results_W_1_imag, file)
        save_list('W_2_imag', results_W_2_imag, file)
        save_list('W_3_imag', results_W_3_imag, file)
        save_list('phase_imag', results_pure_phase_imag, file)


def save_list(row_name, list, file):
    file.write(row_name + ',')
    for entry in list:
        file.write(str(entry) + ",")
    file.write('\n')

if __name__ == "__main__":
    main()
