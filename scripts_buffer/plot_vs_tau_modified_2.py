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
            W_0s = []
            W_0s_real = []
            W_0s_imag = []

            W_1s = []
            W_1s_real = []
            W_1s_imag = []

            W_1s_phased = []
            W_1s_phased_real = []
            W_1s_phased_imag = []

            W_2s = []
            W_2s_real = []
            W_2s_imag = []

            W_2s_phased = []
            W_2s_phased_real = []
            W_2s_phased_imag = []

            W_3s = []
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

                    W_taus.append(W_tau)
                    W_0s.append(W_0)
                    W_0s_real.append(W_0.real)
                    W_0s_imag.append(W_0.imag)
                    W_1s.append(W_1)
                    W_1s_real.append(W_1.real)
                    W_1s_imag.append(W_1.imag)
                    W_1s_phased.append(W_1_phased)
                    W_1s_phased_real.append(W_1_phased.real)
                    W_1s_phased_imag.append(W_1_phased.imag)
                    W_2s.append(W_2)
                    W_2s_real.append(W_2.real)
                    W_2s_imag.append(W_2.imag)
                    W_2s_phased.append(W_2_phased)
                    W_2s_phased_real.append(W_2_phased.real)
                    W_2s_phased_imag.append(W_2_phased.imag)
                    W_3s.append(W_3)
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

            # Cast all data to ndarrays
            taus = np.array(taus)
            g_avs = np.array(g_avs)
            g_avs_no_phase = np.array(g_avs_no_phase)
            g_avs_negative_one = np.array(g_avs_negative_one)
            g_avs_i = np.array(g_avs_i)
            g_avs_negative_i = np.array(g_avs_negative_i)

            W_taus = np.array(W_taus)
            W_0s = np.array(W_0s)
            W_0s_real = np.array(W_0s_real)
            W_0s_imag = np.array(W_0s_imag)
            W_1s = np.array(W_1s)
            W_1s_real = np.array(W_1s_real)
            W_1s_imag = np.array(W_1s_imag)
            W_1s_phased = np.array(W_1s_phased)
            W_1s_phased_real = np.array(W_1s_phased_real)
            W_1s_phased_imag = np.array(W_1s_phased_imag)
            W_2s = np.array(W_2s)
            W_2s_real = np.array(W_2s_real)
            W_2s_imag = np.array(W_2s_imag)
            W_2s_phased = np.array(W_2s_phased)
            W_2s_phased_real = np.array(W_2s_phased_real)
            W_2s_phased_imag = np.array(W_2s_phased_imag)
            W_3s = np.array(W_3s)
            W_3s_real = np.array(W_3s_real)
            W_3s_imag = np.array(W_3s_imag)

            Ds_plus = np.array(Ds_plus)
            Ds_plus_no_phase = np.array(Ds_plus_no_phase)
            Ds_minus = np.array(Ds_minus)
            Ds_minus_no_phase = np.array(Ds_minus_no_phase)
            Ds_t_minus_tau = np.array(Ds_t_minus_tau)

            ps_plus = np.array(ps_plus)
            ps_minus = np.array(ps_minus)
            ps_plus_no_phase = np.array(ps_plus_no_phase)
            ps_minus_no_phase = np.array(ps_minus_no_phase)

            pure_phases_real = np.array(pure_phases_real)
            pure_phases_imag = np.array(pure_phases_imag)
            ####################################################################

            # max_W_1_phased_imag = np.max(W_1s_phased_imag)
            # min_W_2_phased_imag = np.min(W_2s_phased_imag)
            # min_pure_phase_imag = np.min(pure_phases_imag)
            # max_pure_phase_imag = np.max(pure_phases_imag)
            #
            # print("W_1 MAX (imag):")
            # for i in range(len(W_1s_phased_imag)):
            #     if np.isclose(W_1s_phased_imag[i], max_W_1_phased_imag, atol=1e-5):
            #         print("tau: ", taus[i])
            #
            # print("W_2 MIN (imag):")
            # for i in range(len(W_2s_phased_imag)):
            #     if np.isclose(W_2s_phased_imag[i], min_W_2_phased_imag, atol=1e-5):
            #         print("tau: ", taus[i])
            #
            # print("pure phase MIN (imag):")
            # for i in range(len(pure_phases_imag)):
            #     if np.isclose(pure_phases_imag[i], min_pure_phase_imag, atol=1e-5):
            #         print("tau: ", taus[i])
            #
            # print("pure phase MAX (imag):")
            # for i in range(len(pure_phases_imag)):
            #     if np.isclose(pure_phases_imag[i], max_pure_phase_imag, atol=1e-5):
            #         print("tau: ", taus[i])

            print("W_0s: ", W_0s[:10])
            print("W_1s: ", W_1s[:10])
            print("W_2s: ", W_2s[:10])
            print("W_3s: ", W_3s[:10])
            print("W_taus: ", W_taus[:10])
            # sys.exit()

            if extended_mode:
                fig, axs = plt.subplots(nrows=6, ncols=2, figsize=(25,14))

                # First column
                axs[0][0].plot(taus, W_0s_real, "-", label=r'$Re(W_0)$')
                axs[0][0].plot(taus, W_0s_imag, "-", label=r'$Im(W_0)$')

                axs[1][0].plot(taus, W_1s_real, "-", label=r'$Re(W_1)$')
                axs[1][0].plot(taus, W_1s_imag, "-", label=r'$Im(W_1)$')
                axs[1][0].plot(taus, W_1s_phased_real, "-", label=r'$Re(e^{(...)}W_1)$')
                axs[1][0].plot(taus, W_1s_phased_imag, "-", label=r'$Im(e^{(...)}W_1)$')

                axs[2][0].plot(taus, W_2s_real, "-", label=r'$Re(W_2)$')
                axs[2][0].plot(taus, W_2s_imag, "-", label=r'$Im(W_2)$')
                axs[2][0].plot(taus, W_2s_phased_real, "-", label=r'$Re(e^{(...)}W_2)$')
                axs[2][0].plot(taus, W_2s_phased_imag, "-", label=r'$Im(e^{(...)}W_2)$')

                axs[3][0].plot(taus, W_3s_real, "-", label=r'$Re(W_3)$')
                axs[3][0].plot(taus, W_3s_imag, "-", label=r'$Im(W_3)$')

                axs[4][0].plot(taus, pure_phases_real, "-", label=r'$Re(exp(iE\tau)/\hbar))$')
                axs[4][0].plot(taus, pure_phases_imag, "-", label=r'$Im(exp(iE\tau)/\hbar))$')

                axs[0][0].grid()
                axs[0][0].legend()
                axs[1][0].grid()
                axs[1][0].legend()
                axs[2][0].grid()
                axs[2][0].legend()
                axs[3][0].grid()
                axs[3][0].legend()
                axs[4][0].grid()
                axs[4][0].legend()
                axs[4][0].set_xlabel(r'$\tau$')

                # Second column
                # print("Re(W_0 + W_3): ", W_0s_real + W_3s_real)
                axs[0][1].plot(taus, [x.real for x in W_0s + W_3s], "-", label=r'$Re(...),\ simple$')
                axs[0][1].plot(taus, [x.imag for x in W_0s + W_3s], "-", label=r'$Im(...),\ simple$')
                axs[0][1].plot(taus, W_0s_real + W_3s_real, "-", label=r'$Re(...),\ splitted$')
                axs[0][1].plot(taus, W_0s_imag + W_3s_imag, "-", label=r'$Im(...),\ splitted$')
                axs[0][1].set_title(r"$W_0 + W_3$")

                axs[1][1].plot(taus, [x.real for x in W_1s_phased + W_2s_phased], "-", label=r'$Re(...),\ simple$')
                axs[1][1].plot(taus, [x.imag for x in W_1s_phased + W_2s_phased], "-", label=r'$Im(...),\ simple$')
                axs[1][1].plot(taus, W_1s_phased_real + W_2s_phased_real, "-", label=r'$Re(...)$')
                axs[1][1].plot(taus, W_1s_phased_imag + W_2s_phased_imag, "-", label=r'$Im(...)$')
                axs[1][1].set_title(r"$e^{-i...}W_1 + e^{i...}W_2$")

                axs[2][1].plot(taus, W_0s_real + W_1s_phased_real + W_2s_phased_real + W_3s_real, "-", label=r'$Re(...)$')
                axs[2][1].plot(taus, W_0s_imag + W_1s_phased_imag + W_2s_phased_imag + W_3s_imag, "-", label=r'$Im(...)$')
                axs[2][1].set_title(r"$W_0 + e^{-i...}W_1 + e^{i...}W_2 + W_3$")

                Ws_sum_plus = W_0s + W_1s_phased + W_2s_phased + W_3s
                Ws_sum_plus_real = W_0s_real + W_1s_phased_real + W_2s_phased_real + W_3s_real
                Ws_sum_plus_imag = W_0s_imag + W_1s_phased_imag + W_2s_phased_imag + W_3s_imag

                Ws_sum_minus = -W_0s + W_1s_phased + W_2s_phased + -W_3s
                Ws_sum_minus_real = -W_0s_real + W_1s_phased_real + W_2s_phased_real + -W_3s_real
                Ws_sum_minus_imag = -W_0s_imag + W_1s_phased_imag + W_2s_phased_imag + -W_3s_imag

                axs[3][1].plot(taus, 0.25*np.abs(Ws_sum_plus), "-", label="+ simple")
                axs[3][1].plot(taus, 0.25*np.sqrt(Ws_sum_plus_real**2 + Ws_sum_plus_imag**2), "--", label="+ splitted")
                axs[3][1].plot(taus, 0.25*np.abs(Ws_sum_minus), "-", label="- simple")
                axs[3][1].plot(taus, 0.25*np.sqrt(Ws_sum_minus_real**2 + Ws_sum_minus_imag**2), "--", label="- splitted")
                axs[3][1].plot(taus, Ds_t_minus_tau, "-", label=r"$D_{t-\tau}$")
                ################################################################
                axs[3][1].plot(taus, 0.25*np.abs(Ws_sum_plus) + 0.25*np.abs(Ws_sum_minus), "-")
                ################################################################
                axs[3][1].set_title(r"$1/4|W_0 + e^{-i...}W_1 + e^{i...}W_2 + W_3|$")

                axs[4][1].plot(taus, 0.25*np.abs(Ws_sum_plus) / np.abs(ps_plus), "-", label="+ simple")
                axs[4][1].plot(taus, 0.25*np.sqrt(Ws_sum_plus_real**2 + Ws_sum_plus_imag**2) / np.abs(ps_plus), "--", label="+ splitted")
                axs[4][1].plot(taus, 0.25*np.abs(Ws_sum_minus) / np.abs(ps_minus), "-", label="- simple")
                axs[4][1].plot(taus, Ds_t_minus_tau, "-", label=r"$D_{t-\tau}$")
                axs[4][1].set_title(r"$1/4|W_0 + e^{-i...}W_1 + e^{i...}W_2 + W_3| / |p_+|$")

                # axs[5][1].plot(taus, ps_plus, "-", label=r'$p_+$')
                # axs[5][1].plot(taus, pure_phases_real, "--", label=r'$Re(exp(iE\tau)/\hbar))$')
                # axs[5][1].plot(taus, pure_phases_imag, "--", label=r'$Im(exp(iE\tau)/\hbar))$')
                axs[5][1].plot(taus, g_avs, "-", label=r'$g_{av}$')
                axs[5][1].plot(taus, 0.25*np.abs(Ws_sum_plus) + 0.25*np.abs(Ws_sum_minus) - Ds_t_minus_tau, "--", label=r'$g_{av}^{NO\ PROB.}$')


                axs[0][1].grid()
                axs[0][1].legend()
                axs[1][1].grid()
                axs[1][1].legend()
                axs[2][1].grid()
                axs[2][1].legend()
                axs[3][1].grid()
                axs[3][1].legend()
                axs[4][1].grid()
                axs[4][1].legend()
                axs[5][1].grid()
                axs[5][1].legend()

                # filename = 'g_av_vs_tau_T=' + str(T_temp_string) + "_" + t_time_string + "_" + tau_time_string + "_" + "{:.2f}".format(taus[-1]) + '_modified.pdf'
                # plt.savefig(filename)

                plt.show()
                plt.clf()
                sys.exit()
            else:
                plt.plot(taus, g_avs, "-")
                plt.ylabel(r'$g_{av}$')
                plt.xlabel(r'$\tau$')
                plt.grid()
                filename = 'g_av_vs_tau_T=' + str(T_temp_string) + "_" + t_time_string + "_" + tau_time_string + "_" + "{:.2f}".format(taus[-1]) + '_modified.pdf'
                plt.savefig(filename)
                # plt.show()
                plt.clf()


if __name__ == '__main__':
    main()
