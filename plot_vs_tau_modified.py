from matplotlib import pyplot as plt
import numpy as np
import sys, os, glob

def main():
    if len(sys.argv) != 5:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- tau_time_start, ")
        print("- t_time, ")
        print("- T_temperature,")
        print("- MODE: 0 - simple plot, 1 - extended plot.")
        sys.exit()

    tau_time_string = sys.argv[1]
    t_time_string = sys.argv[2]
    T_temp_string = sys.argv[3]
    extended_mode = sys.argv[4]

    if extended_mode != "0" and extended_mode != "1":
        print("Wrong type of mode! PLease provide '0' or '1'.")
        sys.exit()
    else:
        if extended_mode == "0":
            extended_mode = False
        else:
            extended_mode = True

    print("tau_time: ", tau_time_string)
    print("t_time: ", t_time_string)
    print("T_temp: ", T_temp_string)

    os.chdir("./")
    input_file_prefix = "g_av_vs_tau_T=" + T_temp_string + "_" + t_time_string + "_" + tau_time_string + "*"
    input_file_postfix = "_modified.dat"
    for filename in glob.glob(input_file_prefix):
        if input_file_postfix in filename:
            taus = []
            g_avs = []
            g_avs_no_phase = []
            g_avs_negative_one = []
            g_avs_i = []
            g_avs_negative_i = []

            W_taus = []
            W_0s_real = []
            W_0s_imag = []
            W_1s_real = []
            W_1s_imag = []
            W_1s_phased_real = []
            W_1s_phased_imag = []
            W_2s_real = []
            W_2s_imag = []
            W_2s_phased_real = []
            W_2s_phased_imag = []
            W_3s_real = []
            W_3s_imag = []

            Ds_plus = []
            Ds_plus_no_phase = []
            Ds_minus = []
            Ds_minus_no_phase = []
            Ds_t_minus_tau = []

            ps_plus = []
            ps_minus = []
            ps_plus_no_phase = []
            ps_minus_no_phase = []

            pure_phases_real = []
            pure_phases_imag = []

            print("filename: ", filename)
            input_file = open(filename, "r")

            content = input_file.readlines()
            for line in content:
                if line != "\n":
                    splitted_line = line.split()

                    ############################################################
                    tau = float(splitted_line[0])
                    g_av = float(splitted_line[1])
                    W_tau = splitted_line[2].replace('(', '').replace(')', '').split(',')
                    W_tau = (float(W_tau[0]), float(W_tau[1]))
                    W_tau = complex(W_tau[0], W_tau[1])
                    W_0 = splitted_line[3].replace('(', '').replace(')', '').split(',')
                    W_0 = (float(W_0[0]), float(W_0[1]))
                    W_0 = complex(W_0[0], W_0[1])
                    W_1 = splitted_line[4].replace('(', '').replace(')', '').split(',')
                    W_1 = (float(W_1[0]), float(W_1[1]))
                    W_1 = complex(W_1[0], W_1[1])
                    W_1_phased = splitted_line[5].replace('(', '').replace(')', '').split(',')
                    W_1_phased = (float(W_1_phased[0]), float(W_1_phased[1]))
                    W_1_phased = complex(W_1_phased[0], W_1_phased[1])
                    W_2 = splitted_line[6].replace('(', '').replace(')', '').split(',')
                    W_2 = (float(W_2[0]), float(W_2[1]))
                    W_2 = complex(W_2[0], W_2[1])
                    W_2_phased = splitted_line[7].replace('(', '').replace(')', '').split(',')
                    W_2_phased = (float(W_2_phased[0]), float(W_2_phased[1]))
                    W_2_phased = complex(W_2_phased[0], W_2_phased[1])
                    W_3 = splitted_line[8].replace('(', '').replace(')', '').split(',')
                    W_3 = (float(W_3[0]), float(W_3[1]))
                    W_3 = complex(W_3[0], W_3[1])
                    D_plus = float(splitted_line[9])
                    D_minus = float(splitted_line[10])
                    D_t_minus_tau = float(splitted_line[11])
                    p_plus = float(splitted_line[12])
                    p_minus = float(splitted_line[13])
                    p_plus_no_phase = float(splitted_line[14])
                    p_minus_no_phase = float(splitted_line[15])
                    pure_phase = splitted_line[16].replace('(', '').replace(')', '').split(',')
                    pure_phase = (float(pure_phase[0]), float(pure_phase[1]))
                    pure_phase = complex(pure_phase[0], pure_phase[1])

                    ############################################################
                    # No phase at all (which is equal to exp(iE tau / hbar) = 1)
                    ############################################################
                    denominator_plus_no_phase = np.abs(p_plus_no_phase)
                    denominator_minus_no_phase = np.abs(p_minus_no_phase)

                    D_plus_no_phase = np.abs(0.25*(W_0 + W_1 + W_2 + W_3))/denominator_plus_no_phase
                    D_minus_no_phase = np.abs(0.25*(-W_0 + W_1 + W_2 - W_3))/denominator_minus_no_phase

                    g_plus_no_phase = (D_plus_no_phase - D_t_minus_tau)
                    g_minus_no_phase = (D_minus_no_phase - D_t_minus_tau)

                    g_av_no_phase = (p_plus_no_phase*g_plus_no_phase) + (p_minus_no_phase*g_minus_no_phase)

                    ############################################################
                    # exp(iE tau / hbar) = -1
                    ############################################################
                    p_plus_negative_one = 0.5*(1 -(W_tau).real)
                    p_minus_negative_one = 0.5*(1 + (W_tau).real)
                    denominator_plus_negative_one = np.abs(p_plus_negative_one)
                    denominator_minus_negative_one = np.abs(p_minus_negative_one)

                    D_plus_negative_one = np.abs(0.25*(W_0 - W_1 - W_2 + W_3))/denominator_plus_negative_one
                    D_minus_negative_one = np.abs(0.25*(-W_0 - W_1 - W_2 - W_3))/denominator_minus_negative_one

                    g_plus_negative_one = (D_plus_negative_one - D_t_minus_tau)
                    g_minus_negative_one = (D_minus_negative_one - D_t_minus_tau)

                    g_av_negative_one = (p_plus_negative_one*g_plus_negative_one) + (p_minus_negative_one*g_minus_negative_one)

                    ############################################################
                    # exp(iE tau / hbar) = i
                    ############################################################
                    p_plus_i = 0.5*(1 + (1j*W_tau).real)
                    p_minus_i = 0.5*(1 - (1j*W_tau).real)

                    denominator_plus_i = np.abs(p_plus_i)
                    denominator_minus_i = np.abs(p_minus_i)

                    D_plus_i = np.abs(0.25*(W_0 - 1j*W_1 + 1j*W_2 + W_3))/denominator_plus_i
                    D_minus_i = np.abs(0.25*(-W_0 - 1j*W_1 + 1j*W_2 - W_3))/denominator_minus_i

                    g_plus_i = (D_plus_i - D_t_minus_tau)
                    g_minus_i = (D_minus_i - D_t_minus_tau)

                    g_av_i = (p_plus_i*g_plus_i) + (p_minus_i*g_minus_i)

                    ############################################################
                    # exp(iE tau / hbar) = -i
                    ############################################################
                    p_plus_negative_i = 0.5*(1 + (-1j*W_tau).real)
                    p_minus_negative_i = 0.5*(1 - (-1j*W_tau).real)

                    denominator_plus_negative_i = np.abs(p_plus_negative_i)
                    denominator_minus_negative_i = np.abs(p_minus_negative_i)

                    D_plus_negative_i = np.abs(0.25*(W_0 + 1j*W_1 - 1j*W_2 + W_3))/denominator_plus_negative_i
                    D_minus_negative_i = np.abs(0.25*(-W_0 + 1j*W_1 - 1j*W_2 - W_3))/denominator_minus_negative_i

                    g_plus_negative_i = (D_plus_negative_i - D_t_minus_tau)
                    g_minus_negative_i = (D_minus_negative_i - D_t_minus_tau)

                    g_av_negative_i = (p_plus_negative_i*g_plus_negative_i) + (p_minus_negative_i*g_minus_negative_i)
                    ############################################################

                    taus.append(tau)
                    g_avs.append(g_av)
                    g_avs_no_phase.append(g_av_no_phase)
                    g_avs_negative_one.append(g_av_negative_one)
                    g_avs_i.append(g_av_i)
                    g_avs_negative_i.append(g_av_negative_i)

                    W_0s_real.append(W_0.real)
                    W_0s_imag.append(W_0.imag)
                    W_1s_real.append(W_1.real)
                    W_1s_imag.append(W_1.imag)
                    W_1s_phased_real.append(W_1_phased.real)
                    W_1s_phased_imag.append(W_1_phased.imag)
                    W_2s_real.append(W_2.real)
                    W_2s_imag.append(W_2.imag)
                    W_2s_phased_real.append(W_2_phased.real)
                    W_2s_phased_imag.append(W_2_phased.imag)
                    W_3s_real.append(W_3.real)
                    W_3s_imag.append(W_3.imag)

                    Ds_plus.append(D_plus)
                    Ds_minus.append(D_minus)
                    Ds_t_minus_tau.append(D_t_minus_tau)

                    Ds_plus_no_phase.append(D_plus_no_phase)
                    Ds_minus_no_phase.append(D_minus_no_phase)

                    ps_plus.append(p_plus)
                    ps_minus.append(p_minus)
                    ps_plus_no_phase.append(p_plus_no_phase)
                    ps_minus_no_phase.append(p_minus_no_phase)

                    pure_phases_real.append(pure_phase.real)
                    pure_phases_imag.append(pure_phase.imag)
                else:
                    break
            input_file.close()

            if extended_mode:
                fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(12.8,10))
                axs[0][0].plot(taus, g_avs, "-", label=r'$g_{av}$')
                axs[1][0].plot(taus, ps_plus, "-", label=r'$p_+$')
                axs[1][0].plot(taus, ps_minus, "-", label=r'$p_-$')
                axs[2][0].plot(taus, Ds_plus, "-", label=r'$D_+$')
                axs[2][0].plot(taus, Ds_minus, "-", label=r'$D_-$')
                axs[3][0].plot(taus, pure_phases_real, "-", label=r'$Re(exp(iE\tau)/\hbar))$')
                axs[3][0].plot(taus, pure_phases_imag, "-", label=r'$Im(exp(iE\tau)/\hbar))$')

                # plt.title(r'$\tau\ =\ ' + str(tau) + '$')
                axs[0][0].grid()
                axs[0][0].set_ylabel(r'$g_{av}$')
                axs[1][0].grid()
                axs[1][0].legend()
                axs[1][0].set_ylabel(r'$probability$')
                axs[2][0].grid()
                axs[2][0].legend()
                axs[2][0].set_ylabel(r'$D$')
                axs[3][0].grid()
                axs[3][0].legend()
                axs[3][0].set_xlabel(r'$\tau$')

                # Second column which doesn't include phase
                axs[0][1].plot(taus, g_avs, "-", label=r'$g_{av}$')
                axs[0][1].plot(taus, g_avs_no_phase, "--", label=r'$g_{av}\ (e^{iE\tau/\hbar} = +1)$')
                axs[0][1].plot(taus, g_avs_negative_one, "--", label=r'$g_{av}\ (e^{iE\tau/\hbar} = -1)$')
                axs[0][1].plot(taus, g_avs_i, "--", label=r'$g_{av}\ (e^{iE\tau/\hbar} = +i)$')
                axs[0][1].plot(taus, g_avs_negative_i, "--", label=r'$g_{av}\ (e^{iE\tau/\hbar} = -i)$')

                axs[1][1].plot(taus, ps_plus, "-", label=r'$p_+$')
                axs[1][1].plot(taus, ps_minus, "-", label=r'$p_-$')
                axs[1][1].plot(taus, ps_plus_no_phase, "--", label=r'$p_+\ (no\ phase)$')
                axs[1][1].plot(taus, ps_minus_no_phase, "--", label=r'$p_-\ (no\ phase)$')

                axs[2][1].plot(taus, Ds_plus, "-", label=r'$D_+$')
                axs[2][1].plot(taus, Ds_minus, "-", label=r'$D_-$')
                axs[2][1].plot(taus, Ds_plus_no_phase, "--", label=r'$D_+\ (no\ phase)$')
                axs[2][1].plot(taus, Ds_minus_no_phase, "--", label=r'$D_-\ (no\ phase)$')

                axs[3][1].plot(taus, pure_phases_real, "-", label=r'$Re(exp(iE\tau)/\hbar))$')
                axs[3][1].plot(taus, pure_phases_imag, "-", label=r'$Im(exp(iE\tau)/\hbar))$')

                axs[0][1].grid()
                axs[0][1].legend()
                axs[0][1].set_ylabel(r'$g_{av}$')
                axs[1][1].grid()
                axs[1][1].legend()
                axs[1][1].set_ylabel(r'$probability$')
                axs[2][1].grid()
                axs[2][1].legend()
                axs[2][1].set_ylabel(r'$D$')
                axs[3][1].grid()
                axs[3][1].legend()
                axs[3][1].set_xlabel(r'$\tau$')


                # filename = 'g_av_vs_tau_T=' + str(T_temp_string) + "_" + t_time_string + "_" + tau_time_string + "_" + "{:.2f}".format(taus[-1]) + '.pdf'
                # plt.savefig(filename)



                plt.show()
                plt.clf()
                sys.exit()
            else:
                plt.plot(taus, g_avs, "-")
                plt.ylabel(r'$g_{av}$')
                plt.xlabel(r'$\tau$')
                plt.grid()
                filename = 'g_av_vs_tau_T=' + str(T_temp_string) + "_" + t_time_string + "_" + tau_time_string + "_" + "{:.2f}".format(taus[-1]) + '.pdf'
                plt.savefig(filename)
                # plt.show()
                plt.clf()


if __name__ == '__main__':
    main()
