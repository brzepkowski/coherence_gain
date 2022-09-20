from matplotlib import pyplot as plt
import numpy as np
import sys, os, glob
import math

def smooth_data(taus_raw, Ds_plus_raw, Ds_minus_raw, Ds_t_raw, epsilon=5e-3):
    Ds_plus = []
    Ds_minus = []
    Ds_t = []
    taus = []
    # multipliers = np.linspace(0, 2, len(taus_raw))
    multipliers = np.linspace(0, 3, len(taus_raw))
    for i in range(len(taus_raw)):
        tau = taus_raw[i]
        D_plus = Ds_plus_raw[i]
        D_minus = Ds_minus_raw[i]
        D_t = Ds_t_raw[i]
        multiplier = math.exp(-multipliers[i])
        if len(Ds_plus) == 0 or (abs(D_plus - Ds_plus[-1]) < multiplier*epsilon
            and abs(D_minus - Ds_minus[-1]) < multiplier*epsilon
            and abs(D_t - Ds_t[-1]) < multiplier*epsilon):
            taus.append(tau)
            Ds_plus.append(D_plus)
            Ds_minus.append(D_minus)
            Ds_t.append(D_t)

    return taus, Ds_plus, Ds_minus, Ds_t


def main():
    # tau_min_0_string = "0.423879981"
    # tau_max_0_string = "0.424919993"
    tau_min_0_string = "4.0012598"
    tau_max_0_string = "4.00021982"
    T_temp_0_string = "34.0"


    # tau_min_1_string = "0.423879981"
    # tau_max_1_string = "0.424919993"
    tau_min_1_string = "4.0012598"
    tau_max_1_string = "4.00021982"
    T_temp_1_string = "70.0"

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
    Ds_t = []
    ps_plus = []
    ps_minus = []
    pure_phases = []

    for (tau_min_string, tau_max_string, T_temp_string) in [(tau_min_0_string, tau_max_0_string, T_temp_0_string), (tau_min_1_string, tau_max_1_string, T_temp_1_string)]:
        tau_min = float(tau_min_string)
        tau_max = float(tau_max_string)

        filename_postfix = ".dat"
        os.chdir("./")

        # Below we assume, that the values on x axis (the t_time) starts at 0
        filename_tau_min_prefix = "g_av_vs_t_MIN_T=" + T_temp_string + "_" + tau_min_string + "_0.00_*"
        filename_tau_max_prefix = "g_av_vs_t_MAX_T=" + T_temp_string + "_" + tau_max_string + "_0.00_*"

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
                    Ds_t_partial = []
                    ps_plus_partial = []
                    ps_minus_partial = []
                    pure_phases_partial = []

                    content = input_file.readlines()
                    for line in content:
                        if line != "\n":
                            splitted_line = line.split()
                            t_time = float(splitted_line[0])
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
                            D_t = float(splitted_line[10])
                            p_plus = float(splitted_line[11])
                            p_minus = float(splitted_line[12])
                            pure_phase = splitted_line[13].replace("(", "").replace(")", "").split(",")
                            pure_phase = complex(float(pure_phase[0]), float(pure_phase[1]))

                            if t_time >= 7.999:
                                break

                            xs_partial.append(t_time)
                            g_avs_partial.append(g_av)
                            W_0s_partial.append(W_0)
                            W_1s_partial.append(W_1)
                            W_1s_phased_partial.append(W_1_phased)
                            W_2s_partial.append(W_2)
                            W_2s_phased_partial.append(W_2_phased)
                            W_3s_partial.append(W_3)
                            Ds_plus_partial.append(D_plus)
                            Ds_minus_partial.append(D_minus)
                            Ds_t_partial.append(D_t)
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
                    Ds_t.append(Ds_t_partial)
                    ps_plus.append(ps_plus_partial)
                    ps_minus.append(ps_minus_partial)
                    pure_phases.append(pure_phases_partial)
                    # labels.append(r"$\tau_{MIN}=" + tau_min_string + "$")
                    break

        labels.append(r"$\tau_{MIN}=" + tau_min_string + "\ T=" + T_temp_string)
        labels.append(r"$\tau_{MAX}=" + tau_max_string + "\ T=" + T_temp_string)

    fig = plt.figure(figsize=[12, 6])
    font_size = 25 # Changes the size of all fonts in the plot
    tick_size = 25 # Changes the size of all labels on axes in the plot
    plt.rc('font', size=font_size)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    ############################################################################
    # Plot Ds for both tau_min and tau_max
    ############################################################################

    plt.subplot(1, 2, 1)
    # tau_min values
    xs[0], Ds_plus[0], Ds_minus[0], Ds_t[0] = smooth_data(xs[0], Ds_plus[0], Ds_minus[0], Ds_t[0])
    plt.plot(xs[0], Ds_plus[0], "--", label=labels[0] + ", D_+$", color='C0')
    plt.plot(xs[0], Ds_minus[0], "--", label=labels[0] + ", D_-$", color='C0')
    plt.plot(xs[0], Ds_t[0], ":", label=labels[0] + ", D$", color='C2')

    # tau_max values
    xs[1], Ds_plus[1], Ds_minus[1], Ds_t[1] = smooth_data(xs[1], Ds_plus[1], Ds_minus[1], Ds_t[1], epsilon=1e-2)
    plt.plot(xs[1], Ds_plus[1], "-", label=labels[1] + ", D_+$", color='C1')
    # plt.plot(xs[1], Ds_minus[1], "--", label=labels[1] + ", D_-$", color='C1') # This line overlaps with the D_+ plotted one line above
    # plt.plot(xs[1], Ds_t[1], ":", label=labels[1] + ", D$", color='C2') # This line overlaps with the D line plotted 7 lines above

    # plt.annotate(r'$34\ K$', xy =(3.3, 0.9))
    plt.annotate(r'$34\ K$', xy =(6.6, 0.9))
    # plt.title()
    plt.ylabel(r'$D$', fontsize=font_size)
    plt.ylim([0, 1])
    plt.xlabel(r'$t\ [ps]$', fontsize=font_size)
    plt.xlim([min(xs[0]), max(xs[0])])
    plt.tick_params(axis='x', pad=15, labelsize=tick_size)
    plt.tick_params(axis='y', labelsize=tick_size)
    # plt.locator_params(axis='x', nbins=4)
    plt.locator_params(axis='x', nbins=8)
    # plt.legend()
    # plt.grid()

    plt.subplot(1, 2, 2)
    # tau_min values
    xs[2], Ds_plus[2], Ds_minus[2], Ds_t[2] = smooth_data(xs[2], Ds_plus[2], Ds_minus[2], Ds_t[2], epsilon=1e-2)
    plt.plot(xs[2], Ds_plus[2], "--", label=labels[2] + ", D_+$", color='C0')
    plt.plot(xs[2], Ds_minus[2], "--", label=labels[2] + ", D_-$", color='C0')
    plt.plot(xs[2], Ds_t[2], ":", label=labels[2] + ", D$", color='C2')

    # tau_max values
    xs[3], Ds_plus[3], Ds_minus[3], Ds_t[3] = smooth_data(xs[3], Ds_plus[3], Ds_minus[3], Ds_t[3], epsilon=1e-2)
    plt.plot(xs[3], Ds_plus[3], "-", label=labels[3] + ", D_+$", color='C1')
    # plt.plot(xs[3], Ds_minus[3], "--", label=labels[3] + ", D_-$", color='C1') # This line overlaps with the D_+ plotted one line above
    # plt.plot(xs[3], Ds_t[3], ":", label=labels[3] + ", D$", color='C2') # This line overlaps with the D line plotted 7 lines above
    # plt.annotate(r'$70\ K$', xy =(3.3, 0.9))
    plt.annotate(r'$70\ K$', xy =(6.6, 0.9))
    # plt.title()
    # plt.ylabel(r'$D$')
    plt.ylim([0, 1])
    plt.yticks([])
    plt.xlabel(r'$t\ [ps]$', fontsize=font_size)
    plt.xlim([min(xs[2]), max(xs[2])])
    plt.tick_params(axis='x', pad=15, labelsize=tick_size)
    # plt.locator_params(axis='x', nbins=4)
    plt.locator_params(axis='x', nbins=8)
    # plt.legend()
    # plt.grid()


    # plt.tight_layout()
    plt.subplots_adjust(left=0.1,
                    bottom=0.2,
                    right=0.95,
                    top=0.95,
                    wspace=0.02)
    filename_plot = "D_vs_t_tau=" + str(round(float(tau_min_0_string),2)) + "_subplots.pdf"
    plt.savefig(filename_plot)
    # plt.show()

if __name__ == '__main__':
    main()
