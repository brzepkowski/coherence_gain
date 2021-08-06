import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from coherence_gain_t_finite import calc_g_av

def main():
    if len(sys.argv) != 5:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- the temperature 'T' [K],")
        print("- time 't' [ps],")
        print("- 'tau_min' [ps],")
        print("- 'tau_max' [ps].")
        sys.exit()

    T = float(sys.argv[1])
    t = float(sys.argv[2])
    tau_min = float(sys.argv[3])
    tau_max = float(sys.argv[4])

    if tau_min <= 0.0001:
        print("tau_min <= 0.0001. It might cause numerical problems. Provide larger tau_min.")
        sys.exit()
    if tau_min >= tau_max:
        print("tau_min >= tau_max. Provide smaller tau_min or larger tau_max.")
        sys.exit()
    if t <= tau_max:
        print("t <= tau_max. Provide larger value of 't'.")
        sys.exit()

    tau_step = 0.0001

    taus = list(np.arange(tau_min, tau_max + tau_step, tau_step)) # List of all "tau"s
    results_D_plus = []
    results_D_minus = []
    results_D_t = []
    results_p_plus = []
    results_p_minus = []
    results_g_av = []

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
        print("tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'), end="\r")
        # D_plus, D_minus, D_t, p_plus, p_minus, g_av = calc_g_av(t, tau, T)
        W_0, phase_W_1, phase_W_2, W_3, pure_phase, D_plus, D_minus, D_t, p_plus, p_minus, g_av = calc_g_av(t, tau, T)
        results_D_plus.append(D_plus)
        results_D_minus.append(D_minus)
        results_D_t.append(D_t)

        results_p_plus.append(p_plus)
        results_p_minus.append(p_minus)

        results_g_av.append(g_av)

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
    print("tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'))

    plt.figure(figsize=(25,4))
    plt.plot(taus, results_D_plus, "-", label=r'$D_+$')
    plt.plot(taus, results_D_minus, "-", label=r'$D_-$')
    plt.plot(taus, results_D_t, "-", label=r'$D_t$')
    plt.title(r'$t\ =\ ' + str(t) + '$')
    plt.ylabel(r'$D$')
    plt.xlabel(r'$\tau$')
    plt.legend()
    plt.grid()
    filename = 'Ds_T=' + str(T) + '_t=' + str(t) + '_tau_min=' + str(tau_min) + '_tau_max=' + str(tau_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    plt.figure(figsize=(25,4))
    plt.plot(taus, results_p_plus, "-", label=r'$p_+$')
    plt.plot(taus, results_p_minus, "-", label=r'$p_-$')
    plt.title(r'$t\ =\ ' + str(t) + '$')
    plt.ylabel(r'$p$')
    plt.xlabel(r'$\tau$')
    plt.legend()
    plt.grid()
    filename = 'probabilities_T=' + str(T) + '_t=' + str(t) + '_tau_min=' + str(tau_min) + '_tau_max=' + str(tau_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    plt.figure(figsize=(25,4))
    plt.plot(taus, results_g_av, "-")
    plt.title(r'$t\ =\ ' + str(t) + '$')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    filename = 'g_av_T=' + str(T) + '_t=' + str(t) + '_tau_min=' + str(tau_min) + '_tau_max=' + str(tau_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    plt.figure(figsize=(25,4))
    plt.plot(taus, results_W_0_real, "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(taus, results_W_1_real, "-", label=r'$e^{-iE\tau/\hbar}\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(taus, results_W_2_real, "-", label=r'$\langle e^{iE\tau/\hbar}W(\tau)W(t-\tau)\rangle$')
    plt.plot(taus, results_W_3_real, "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(taus, results_pure_phase_real, "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$t\ =\ ' + str(t) + '\ (real)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    plt.legend()
    filename = 'components_real_T=' + str(T) + '_t=' + str(t) + '_tau_min=' + str(tau_min) + '_tau_max=' + str(tau_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    plt.figure(figsize=(25,4))
    plt.plot(taus, results_W_0_imag, "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(taus, results_W_1_imag, "-", label=r'$e^{-iE\tau/\hbar}\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(taus, results_W_2_imag, "-", label=r'$\langle e^{iE\tau/\hbar}W(\tau)W(t-\tau)\rangle$')
    plt.plot(taus, results_W_3_imag, "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(taus, results_pure_phase_imag, "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$t\ =\ ' + str(t) + '\ (imag)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    plt.legend()
    filename = 'components_imag_T=' + str(T) + '_t=' + str(t) + '_tau_min=' + str(tau_min) + '_tau_max=' + str(tau_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

if __name__ == "__main__":
    main()
