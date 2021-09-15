from matplotlib import pyplot as plt
import numpy as np
import glob, os
import sys

def string_to_complex(str):
    str = str.replace("(","").replace(")","").split(",")
    real_part = float(str[0])
    imag_part = float(str[1])
    return complex(real_part, imag_part)

def main():
    if len(sys.argv) != 3:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- the temperature 'T' [K],")
        print("- prefix ('MIN' or 'MAX').")
        sys.exit()

    T_string = sys.argv[1]
    T = float(T_string)
    prefix = sys.argv[2].lower()

    taus = []

    results_D_plus_all_taus = {}
    results_D_minus_all_taus = {}
    results_D_t_all_taus = {}
    results_p_plus_all_taus = {}
    results_p_minus_all_taus = {}
    results_g_av_all_taus = {}

    results_W_0_real_all_taus = {}
    results_W_1_real_all_taus = {}
    results_W_2_real_all_taus = {}
    results_W_3_real_all_taus = {}

    results_W_0_imag_all_taus = {}
    results_W_1_imag_all_taus = {}
    results_W_2_imag_all_taus = {}
    results_W_3_imag_all_taus = {}

    results_W_0_phased_real_all_taus = {}
    results_W_1_phased_real_all_taus = {}
    results_W_2_phased_real_all_taus = {}
    results_W_3_phased_real_all_taus = {}

    results_W_0_phased_imag_all_taus = {}
    results_W_1_phased_imag_all_taus = {}
    results_W_2_phased_imag_all_taus = {}
    results_W_3_phased_imag_all_taus = {}

    results_phase_real_all_taus = {}
    results_phase_imag_all_taus = {}

    os.chdir("./")
    for filename in glob.glob(prefix + "_g_av_*.dat"):
        if '_T=' + str(T) + '_' in filename and '.dat' in filename:
            tau = float(filename.split('_')[4].replace('tau=', '').replace('.dat', ''))
            taus.append(tau)

            input_file = open(filename, "r")
            content = input_file.readlines()

            ts = []

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

            results_W_0_phased_real = []
            results_W_1_phased_real = []
            results_W_2_phased_real = []
            results_W_3_phased_real = []

            results_W_0_phased_imag = []
            results_W_1_phased_imag = []
            results_W_2_phased_imag = []
            results_W_3_phased_imag = []

            results_phase_real = []
            results_phase_imag = []

            for line in content:
                # Read values from file
                splitted_line = line.split()
                t = float(splitted_line[0])
                W_0 = string_to_complex(splitted_line[1])
                W_1 = string_to_complex(splitted_line[2])
                W_1_phased = string_to_complex(splitted_line[3])
                W_2 = string_to_complex(splitted_line[4])
                W_2_phased = string_to_complex(splitted_line[5])
                W_3 = string_to_complex(splitted_line[6])
                phase = string_to_complex(splitted_line[7])
                D_plus = float(splitted_line[8])
                D_minus = float(splitted_line[9])
                D_t = float(splitted_line[10])
                p_plus = float(splitted_line[11])
                p_minus = float(splitted_line[12])
                g_av = float(splitted_line[13])

                # Save results to lists
                ts.append(t-tau)

                results_D_plus.append(D_plus)
                results_D_minus.append(D_minus)
                results_D_t.append(D_t)

                results_p_plus.append(p_plus)
                results_p_minus.append(p_minus)

                results_g_av.append(g_av)

                results_W_0_phased_real.append(W_0.real)
                results_W_0_phased_imag.append(W_0.imag)
                results_W_1_phased_real.append(W_1_phased.real)
                results_W_1_phased_imag.append(W_1_phased.imag)
                results_W_2_phased_real.append(W_2_phased.real)
                results_W_2_phased_imag.append(W_2_phased.imag)
                results_W_3_phased_real.append(W_3.real)
                results_W_3_phased_imag.append(W_3.imag)

                results_phase_real.append(phase.real)
                results_phase_imag.append(phase.imag)

                results_W_0_real.append(W_0.real)
                results_W_0_imag.append(W_0.imag)
                results_W_1_real.append(W_1.real)
                results_W_1_imag.append(W_1.imag)
                results_W_2_real.append(W_2.real)
                results_W_2_imag.append(W_2.imag)
                results_W_3_real.append(W_3.real)
                results_W_3_imag.append(W_3.imag)
            input_file.close()

            results_D_plus_all_taus[tau] = results_D_plus
            results_D_minus_all_taus[tau] = results_D_minus
            results_D_t_all_taus[tau] = results_D_t
            results_p_plus_all_taus[tau] = results_p_plus
            results_p_minus_all_taus[tau] = results_p_minus
            results_g_av_all_taus[tau] = results_g_av
            results_W_0_real_all_taus[tau] = results_W_0_real
            results_W_1_real_all_taus[tau] = results_W_1_real
            results_W_2_real_all_taus[tau] = results_W_2_real
            results_W_3_real_all_taus[tau] = results_W_3_real
            results_W_0_imag_all_taus[tau] = results_W_0_imag
            results_W_1_imag_all_taus[tau] = results_W_1_imag
            results_W_2_imag_all_taus[tau] = results_W_2_imag
            results_W_3_imag_all_taus[tau] = results_W_3_imag
            results_W_0_phased_real_all_taus[tau] = results_W_0_phased_real
            results_W_1_phased_real_all_taus[tau] = results_W_1_phased_real
            results_W_2_phased_real_all_taus[tau] = results_W_2_phased_real
            results_W_3_phased_real_all_taus[tau] = results_W_3_phased_real
            results_W_0_phased_imag_all_taus[tau] = results_W_0_phased_imag
            results_W_1_phased_imag_all_taus[tau] = results_W_1_phased_imag
            results_W_2_phased_imag_all_taus[tau] = results_W_2_phased_imag
            results_W_3_phased_imag_all_taus[tau] = results_W_3_phased_imag
            results_phase_real_all_taus[tau] = results_phase_real
            results_phase_imag_all_taus[tau] = results_phase_imag

    ############################################################################
    # Generate plots
    ############################################################################

    # Sort all taus before plotting
    taus = sorted(taus)

    """
    plt.plot(ts, results_D_plus[:,-1], "-", label=r'$D_+$')
    plt.plot(ts, results_D_minus[:,-1], "-", label=r'$D_-$')
    plt.plot(ts, results_D_t[:,-1], "-", label=r'$D_t$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '$')
    plt.ylabel(r'$D$')
    plt.xlabel(r'$t$')
    plt.legend()
    plt.grid()
    # filename = 'Ds_T=' + str(T) + '_tau=' + str(tau) +'.pdf'
    # plt.savefig(filename)
    plt.show()

    # Clear the plot
    plt.clf()

    plt.plot(ts, results_p_plus[:,-1], "-", label=r'$p_+$')
    plt.plot(ts, results_p_minus[:,-1], "-", label=r'$p_-$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '$')
    plt.ylabel(r'$p$')
    plt.xlabel(r'$t$')
    plt.legend()
    plt.grid()
    # filename = 'probabilities_T=' + str(T) + '_tau=' + '.pdf'
    # plt.savefig(filename)
    plt.show()

    # Clear the plot
    plt.clf()
    """

    for tau in taus:
        plt.plot(ts, results_g_av_all_taus[tau], "-", label=r'$\tau=' + str(tau) + '$')
    plt.title('Average coherence gain (' + prefix.upper() + ')')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t-\tau$')
    # plt.ylim([0.23, 1.02])
    plt.grid()
    plt.legend()
    filename = prefix + '_g_av_T=' + str(T) + '.pdf'
    plt.savefig(filename)
    plt.show()
    """
    # Clear the plot
    plt.clf()

    plt.plot(ts, results_W_0_phased_real[:,-1], "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_W_1_phased_real[:,-1], "-", label=r'$e^{-iE\tau/\hbar} \langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_W_2_phased_real[:,-1], "-", label=r'$e^{iE\tau/\hbar} \langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_W_3_phased_real[:,-1], "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_phase_real[:,-1], "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (real)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t$')
    plt.grid()
    plt.legend()
    # filename = 'components_real_T=' + str(T) + '_tau=' + str(tau) + '.pdf'
    # plt.savefig(filename)
    plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_W_0_phased_imag[:,-1], "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_W_1_phased_imag[:,-1], "-", label=r'$e^{-iE\tau/\hbar} \langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_W_2_phased_imag[:,-1], "-", label=r'$e^{iE\tau/\hbar} \langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_W_3_phased_imag[:,-1], "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_phase_imag[:,-1], "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (imag)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t$')
    plt.grid()
    plt.legend()
    # filename = 'components_imag_T=' + str(T) + '_tau=' + str(tau) + '.pdf'
    # plt.savefig(filename)
    plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_W_0_real[:,-1], "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_W_1_real[:,-1], "-", label=r'$\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_W_2_real[:,-1], "-", label=r'$\langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_W_3_real[:,-1], "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_phase_real[:,-1], "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (real)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t$')
    plt.grid()
    plt.legend()
    # filename = 'components_no_phase_real_T=' + str(T) + '_tau=' + str(tau) + '.pdf'
    # plt.savefig(filename)
    plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_W_0_imag[:,-1], "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_W_1_imag[:,-1], "-", label=r'$\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_W_2_imag[:,-1], "-", label=r'$\langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_W_3_imag[:,-1], "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_phase_imag[:,-1], "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (imag)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t$')
    plt.grid()
    plt.legend()
    # filename = 'components_no_phase_imag_T=' + str(T) + '_tau=' + str(tau) + '.pdf'
    # plt.savefig(filename)
    plt.show()
    """

if __name__ == '__main__':
    main()
