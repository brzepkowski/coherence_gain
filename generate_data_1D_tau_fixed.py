import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from coherence_gain_t_finite_detailed import calc_g_av

def main():
    if len(sys.argv) != 5:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- the temperature 'T' [K],")
        print("- time 'tau' [ps],")
        print("- 't_min' [ps],")
        print("- 't_max' [ps].")
        sys.exit()

    T = float(sys.argv[1])
    tau = float(sys.argv[2])
    t_min = float(sys.argv[3])
    t_max = float(sys.argv[4])

    if t_min - tau <= 0.0001:
        print("t_min - tau <= 0.0001. It might cause numerical problems. Provide larger t_min.")
        sys.exit()
    if t_min >= t_max:
        print("t_min >= t_max. Provide smaller t_min or larger t_max.")
        sys.exit()
    if tau >= t_min:
        print("tau >= t_min. Provide larger value of t_min.")
        sys.exit()

    t_step = 0.0001

    ts = list(np.arange(t_min, t_max + t_step, t_step)) # List of all "ts"
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

    results_phase_W_0_real = []
    results_phase_W_1_real = []
    results_phase_W_2_real = []
    results_phase_W_3_real = []

    results_phase_W_0_imag = []
    results_phase_W_1_imag = []
    results_phase_W_2_imag = []
    results_phase_W_3_imag = []

    results_pure_phase_real = []
    results_pure_phase_imag = []

    for t in ts:
        print("t: ", format(t, '.5f'), " / ", format(ts[-1], '.5f'), end="\r")
        W_0, W_1, phase_W_1, W_2, phase_W_2, W_3, pure_phase, D_plus, D_minus, D_t, p_plus, p_minus, g_av = calc_g_av(t, tau, T)

        results_D_plus.append([tau, t, D_plus])
        results_D_minus.append([tau, t, D_minus])
        results_D_t.append([tau, t, D_t])

        results_p_plus.append([tau, t, p_plus])
        results_p_minus.append([tau, t, p_minus])

        results_g_av.append([tau, t, g_av])

        results_phase_W_0_real.append([tau, t, W_0.real])
        results_phase_W_0_imag.append([tau, t, W_0.imag])
        results_phase_W_1_real.append([tau, t, phase_W_1.real])
        results_phase_W_1_imag.append([tau, t, phase_W_1.imag])
        results_phase_W_2_real.append([tau, t, phase_W_2.real])
        results_phase_W_2_imag.append([tau, t, phase_W_2.imag])
        results_phase_W_3_real.append([tau, t, W_3.real])
        results_phase_W_3_imag.append([tau, t, W_3.imag])

        results_pure_phase_real.append([tau, t, pure_phase.real])
        results_pure_phase_imag.append([tau, t, pure_phase.imag])

        results_W_0_real.append([tau, t, W_0.real])
        results_W_0_imag.append([tau, t, W_0.imag])
        results_W_1_real.append([tau, t, W_1.real])
        results_W_1_imag.append([tau, t, W_1.imag])
        results_W_2_real.append([tau, t, W_2.real])
        results_W_2_imag.append([tau, t, W_2.imag])
        results_W_3_real.append([tau, t, W_3.real])
        results_W_3_imag.append([tau, t, W_3.imag])

        """
        print("=== REAL (phase) ===")
        print("W_0: ", W_0.real)
        print("phase_W_1: ", phase_W_1.real)
        print("phase_W_2: ", phase_W_2.real)
        print("W_3: ", W_3.real)

        print("=== IMAG (phase) ===")
        print("W_0: ", W_0.imag)
        print("phase_W_1: ", phase_W_1.imag)
        print("phase_W_2: ", phase_W_2.imag)
        print("W_3: ", W_3.imag)

        print("=== REAL ===")
        print("W_0: ", W_0.real)
        print("W_1: ", W_1.real)
        print("W_2: ", W_2.real)
        print("W_3: ", W_3.real)

        print("=== IMAG ===")
        print("W_0: ", W_0.imag)
        print("W_1: ", W_1.imag)
        print("W_2: ", W_2.imag)
        print("W_3: ", W_3.imag)
        """
    print("t: ", format(t, '.5f'), " / ", format(ts[-1], '.5f'))

    data = {}

    data["results_D_plus"] = results_D_plus
    data["results_D_minus"] = results_D_minus
    data["results_D_t"] = results_D_t
    data["results_p_plus"] = results_p_plus
    data["results_p_minus"] = results_p_minus
    data["results_g_av"] = results_g_av
    data["results_phase_W_0_real"] = results_phase_W_0_real
    data["results_phase_W_0_imag"] = results_phase_W_0_imag
    data["results_phase_W_1_real"] = results_phase_W_1_real
    data["results_phase_W_1_imag"] = results_phase_W_1_imag
    data["results_phase_W_2_real"] = results_phase_W_2_real
    data["results_phase_W_2_imag"] = results_phase_W_2_imag
    data["results_phase_W_3_real"] = results_phase_W_3_real
    data["results_phase_W_3_imag"] = results_phase_W_3_imag
    data["results_pure_phase_real"] = results_pure_phase_real
    data["results_pure_phase_imag"] = results_pure_phase_imag
    data["results_W_0_real"] = results_W_0_real
    data["results_W_0_imag"] = results_W_0_imag
    data["results_W_1_real"] = results_W_1_real
    data["results_W_1_imag"] = results_W_1_imag
    data["results_W_2_real"] = results_W_2_real
    data["results_W_2_imag"] = results_W_2_imag
    data["results_W_3_real"] = results_W_3_real
    data["results_W_3_imag"] = results_W_3_imag

    filename = 'data_T=' + str(T) + '_tau=' + str(tau) + '_t_min=' + str(t_min) + '_t_max=' + str(t_max) + '.pickle'
    with open(filename, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    sys.exit()

    ############################################################################

    """
    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_D_plus, "-", label=r'$D_+$')
    plt.plot(ts, results_D_minus, "-", label=r'$D_-$')
    plt.plot(ts, results_D_t, "-", label=r'$D_t$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '$')
    plt.ylabel(r'$D$')
    plt.xlabel(r'$\tau$')
    plt.legend()
    plt.grid()
    filename = 'Ds_T=' + str(T) + '_tau=' + str(tau) + '_t_min=' + str(t_min) + '_t_max=' + str(t_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_p_plus, "-", label=r'$p_+$')
    plt.plot(ts, results_p_minus, "-", label=r'$p_-$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '$')
    plt.ylabel(r'$p$')
    plt.xlabel(r'$\tau$')
    plt.legend()
    plt.grid()
    filename = 'probabilities_T=' + str(T) + '_tau=' + str(tau) + '_t_min=' + str(t_min) + '_t_max=' + str(t_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_g_av, "-")
    plt.title(r'$\tau\ =\ ' + str(tau) + '$')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    filename = 'g_av_T=' + str(T) + '_tau=' + str(tau) + '_t_min=' + str(t_min) + '_t_max=' + str(t_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_phase_W_0_real, "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_phase_W_1_real, "-", label=r'$e^{-iE\tau/\hbar} \langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_phase_W_2_real, "-", label=r'$e^{iE\tau/\hbar} \langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_phase_W_3_real, "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_pure_phase_real, "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (real)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    plt.legend()
    filename = 'components_real_T=' + str(T) + '_tau=' + str(tau) + '_t_min=' + str(t_min) + '_t_max=' + str(t_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_phase_W_0_imag, "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_phase_W_1_imag, "-", label=r'$e^{-iE\tau/\hbar} \langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_phase_W_2_imag, "-", label=r'$e^{iE\tau/\hbar} \langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_phase_W_3_imag, "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_pure_phase_imag, "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (imag)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    plt.legend()
    filename = 'components_imag_T=' + str(T) + '_tau=' + str(tau) + '_t_min=' + str(t_min) + '_t_max=' + str(t_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_W_0_real, "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_W_1_real, "-", label=r'$\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_W_2_real, "-", label=r'$\langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_W_3_real, "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_pure_phase_real, "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (real)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    plt.legend()
    filename = 'components_no_phase_real_T=' + str(T) + '_tau=' + str(tau) + '_t_min=' + str(t_min) + '_t_max=' + str(t_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_W_0_imag, "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_W_1_imag, "-", label=r'$\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_W_2_imag, "-", label=r'$\langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_W_3_imag, "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_pure_phase_imag, "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (imag)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    plt.legend()
    filename = 'components_no_phase_imag_T=' + str(T) + '_tau=' + str(tau) + '_t_min=' + str(t_min) + '_t_max=' + str(t_max) + '.pdf'
    plt.savefig(filename)
    # plt.show()
    """

if __name__ == "__main__":
    main()
