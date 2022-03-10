from matplotlib import pyplot as plt
import numpy as np
import sys, os, glob

def main():
    if len(sys.argv) != 4:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- tau_MIN,")
        print("- tau_MAX,")
        print("- T_temperature.")
        sys.exit()

    tau_min_string = sys.argv[1]
    tau_max_string = sys.argv[2]
    T_temp_string = sys.argv[3]
    tau_min = float(tau_min_string)
    tau_max = float(tau_max_string)

    print("tau_min_string: ", tau_min_string)
    print("tau_max_string: ", tau_max_string)
    print("T_temp: ", T_temp_string)

    xs = []
    g_avs = []
    labels = []
    W_0s = []
    W_1s = []
    W_1s_phased = []
    W_2s = []
    W_2s_phased = []
    W_3s = []
    Ds_plus = []
    Ds_minus = []
    Ds_t_minus_tau = []
    ps_plus = []
    ps_minus = []
    pure_phases = []

    filename_postfix = ".dat"
    os.chdir("./")

    ############################################################################
    # tau_min
    ############################################################################
    filename_tau_min_prefix = "g_av_vs_t_MIN_T=" + T_temp_string + "_" + tau_min_string + "_" + tau_min_string + "_*"
    filename_tau_max_prefix = "g_av_vs_t_MAX_T=" + T_temp_string + "_" + tau_max_string + "_" + tau_max_string + "_*"

    for filename_prefix in [filename_tau_min_prefix, filename_tau_max_prefix]:
        print("filename_prefix: ", filename_prefix)
        # Get all target filenames and sort them before final plotting
        for filename in glob.glob(filename_prefix):
            if filename_postfix in filename:
                print("filename: ", filename)
                input_file = open(filename, "r")

                xs_partial = []
                g_avs_partial = []
                W_0s_partial = []
                W_1s_partial = []
                W_1s_phased_partial = []
                W_2s_partial = []
                W_2s_phased_partial = []
                W_3s_partial = []
                Ds_plus_partial = []
                Ds_minus_partial = []
                Ds_t_minus_tau_partial = []
                ps_plus_partial = []
                ps_minus_partial = []
                pure_phases_partial = []

                content = input_file.readlines()
                for line in content:
                    if line != "\n":
                        splitted_line = line.split()
                        t_time = float(splitted_line[0])
                        t_minus_tau = t_time - tau_min
                        g_av = float(splitted_line[1])
                        # print("splitted_line[2]: ", splitted_line[2])
                        W_0 = splitted_line[2].replace("(", "").replace(")", "").split(",")
                        W_0 = complex(float(W_0[0]), float(W_0[1]))
                        # print("W_0: ", W_0)
                        # print("splitted_line[3]: ", splitted_line[3])
                        W_1 = splitted_line[3].replace("(", "").replace(")", "").split(",")
                        W_1 = complex(float(W_1[0]), float(W_1[1]))
                        # print("W_1: ", W_1)
                        W_1_phased = splitted_line[4].replace("(", "").replace(")", "").split(",")
                        W_1_phased = complex(float(W_1_phased[0]), float(W_1_phased[1]))
                        W_2 = splitted_line[5].replace("(", "").replace(")", "").split(",")
                        W_2 = complex(float(W_2[0]), float(W_2[1]))
                        W_2_phased = splitted_line[6].replace("(", "").replace(")", "").split(",")
                        W_2_phased = complex(float(W_2_phased[0]), float(W_2_phased[1]))
                        W_3 = splitted_line[7].replace("(", "").replace(")", "").split(",")
                        W_3 = complex(float(W_3[0]), float(W_3[1]))
                        D_plus = float(splitted_line[8])
                        D_minus = float(splitted_line[9])
                        D_t_minus_tau = float(splitted_line[10])
                        p_plus = float(splitted_line[11])
                        p_minus = float(splitted_line[12])
                        pure_phase = splitted_line[13].replace("(", "").replace(")", "").split(",")
                        pure_phase = complex(float(pure_phase[0]), float(pure_phase[1]))


                        epsilon = 1e-4
                        if len(g_avs_partial) == 0 or (abs(g_av - g_avs_partial[-1]) < epsilon and abs(D_plus - Ds_plus_partial[-1]) < epsilon and abs(D_minus - Ds_minus_partial[-1]) < epsilon):
                            xs_partial.append(t_minus_tau)
                            g_avs_partial.append(g_av)
                            W_0s_partial.append(W_0)
                            W_1s_partial.append(W_1)
                            W_1s_phased_partial.append(W_1_phased)
                            W_2s_partial.append(W_2)
                            W_2s_phased_partial.append(W_2_phased)
                            W_3s_partial.append(W_3)
                            Ds_plus_partial.append(D_plus)
                            Ds_minus_partial.append(D_minus)
                            Ds_t_minus_tau_partial.append(D_t_minus_tau)
                            ps_plus_partial.append(p_plus)
                            ps_minus_partial.append(p_minus)
                            pure_phases_partial.append(pure_phase)
                    else:
                        break
                input_file.close()
                xs.append(xs_partial)
                g_avs.append(g_avs_partial)
                W_0s.append(W_0s_partial)
                W_1s.append(W_1s_partial)
                W_1s_phased.append(W_1s_phased_partial)
                W_2s.append(W_2s_partial)
                W_2s_phased.append(W_2s_phased_partial)
                W_3s.append(W_3s_partial)
                Ds_plus.append(Ds_plus_partial)
                Ds_minus.append(Ds_minus_partial)
                Ds_t_minus_tau.append(Ds_t_minus_tau_partial)
                ps_plus.append(ps_plus_partial)
                ps_minus.append(ps_minus_partial)
                pure_phases.append(pure_phases_partial)
                # labels.append(r"$\tau_{MIN}=" + tau_min_string + "$")
                break

    labels.append(r"$\tau_{MIN}=" + tau_min_string)
    labels.append(r"$\tau_{MAX}=" + tau_max_string)

    ############################################################################
    # Plot g_av for both tau_min and tau_max
    ############################################################################

    plt.plot(xs[0], g_avs[0], "-", label=labels[0] + "$")
    plt.plot(xs[1], g_avs[1], "-", label=labels[1] + "$")
    # plt.title()
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t - \tau$')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    filename_plot = 'g_av_vs_t_MIN_MAX_T=' + T_temp_string + "_" + tau_min_string + "_" + tau_max_string + ".pdf"
    # plt.savefig(filename_plot)
    plt.show()
    plt.clf()

    ############################################################################
    # Plot Ds for both tau_min and tau_max
    ############################################################################

    plt.plot(xs[0], Ds_plus[0], "-", label=labels[0] + ", D_+$")
    plt.plot(xs[0], Ds_minus[0], "-", label=labels[0] + ", D_-$")
    plt.plot(xs[0], Ds_t_minus_tau[0], "-", label=labels[0] + ", D$")
    plt.plot(xs[1], Ds_plus[1], "-", label=labels[1] + ", D_+$")
    plt.plot(xs[1], Ds_minus[1], "-", label=labels[1] + ", D_-$")
    plt.plot(xs[1], Ds_t_minus_tau[1], "-", label=labels[1] + ", D$")
    # plt.title()
    plt.ylabel(r'$D$')
    plt.xlabel(r'$t - \tau$')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    filename_plot = 'g_av_vs_t_MIN_MAX_T=' + T_temp_string + "_" + tau_min_string + "_" + tau_max_string + ".pdf"
    # plt.savefig(filename_plot)
    plt.show()

if __name__ == '__main__':
    main()
