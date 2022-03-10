from matplotlib import pyplot as plt
import sys, os, glob
import numpy as np
import math
from envelopes import envelopes

def main():
    tau_time_string = "0.00"
    t_time_string = "0.5"
    T_temp_string_1 = "34.0"
    T_temp_string_2 = "70.0"

    all_taus = []
    all_g_avs = []
    all_labels = []

    for T_temp_string in [T_temp_string_1, T_temp_string_2]:
        os.chdir("./")
        input_file_prefix = "g_av_vs_tau_T=" + T_temp_string + "_" + t_time_string + "_" + tau_time_string + "*"
        input_file_postfix = ".dat"
        for filename in glob.glob(input_file_prefix):
            if input_file_postfix in filename:
                g_avs = []
                taus = []
                D_pluses = []
                D_minuses = []
                D_t_minus_taus = []
                p_pluses = []
                p_minuses = []
                pure_phases_real = []
                pure_phases_imag = []

                print("filename: ", filename)
                input_file = open(filename, "r")
                iter = 0

                content = input_file.readlines()
                for line in content:
                    if line != "\n":
                        splitted_line = line.split()
                        tau = float(splitted_line[0])
                        g_av = float(splitted_line[1])
                        D_plus = float(splitted_line[2])
                        D_minus = float(splitted_line[3])
                        D_t_minus_tau = float(splitted_line[4])
                        p_plus = float(splitted_line[5])
                        p_minus = float(splitted_line[6])
                        pure_phase = splitted_line[7].replace('(', '').replace(')', '').split(',')
                        pure_phase = (float(pure_phase[0]), float(pure_phase[1]))
                        pure_phase = complex(pure_phase[0], pure_phase[1])
                        # print("pure_phase: ", pure_phase)

                        taus.append(tau)
                        g_avs.append(g_av)
                        D_pluses.append(D_plus)
                        D_minuses.append(D_minus)
                        D_t_minus_taus.append(D_t_minus_tau)
                        p_pluses.append(p_plus)
                        p_minuses.append(p_minus)
                        pure_phases_real.append(pure_phase.real)
                        pure_phases_imag.append(pure_phase.imag)
                    else:
                        break
                input_file.close()

                all_taus.append(taus)
                all_g_avs.append(g_avs)
                all_labels.append("T = " + T_temp_string + " K, tau = " + tau_time_string)

    # fig = plt.figure(figsize=[10, 12])
    plt.rc('font', size=25)

    # plt.subplot(2, 1, 1)
    # plt.plot(all_taus[0], all_g_avs[0], '-', label=all_labels[0])
    # plt.plot(all_taus[1], all_g_avs[1], '-', label=all_labels[1])


    # Extract envelopes
    min_taus_g_av, min_g_avs, max_taus_g_av, max_g_avs = envelopes(all_taus[0], all_g_avs[0], epsilon_min=1e-7, epsilon_max=1e-7, diff=0.0001)
    plt.plot(max_taus_g_av, max_g_avs, '-', color='C0')
    plt.plot(min_taus_g_av, min_g_avs, '-', color='C0')

    # Extract envelopes
    min_taus_g_av, min_g_avs, max_taus_g_av, max_g_avs = envelopes(all_taus[1], all_g_avs[1], epsilon_min=1e-7, epsilon_max=1e-7, diff=1.0001)
    plt.plot(max_taus_g_av, max_g_avs, '--', color='C1')
    plt.plot(min_taus_g_av, min_g_avs, '--', color='C1')


    plt.annotate('34 K', xy =(0.23, 0.0044))
    plt.annotate('70 K', xy =(0.33, 0.0085))
    # plt.hlines(0, 0.0, 0.2, color='black', linestyle='dashed')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau\ [ps]$')
    plt.xlim([min(all_taus[0]), max(all_taus[0])])
    plt.locator_params(axis='x', nbins=5)
    plt.tick_params(axis='x', pad=15)
    # plt.ylim([-0.1, 0.1])
    plt.locator_params(axis='y', nbins=5)


    plt.tight_layout()
    filename = 'g_av_vs_tau_t=' + t_time_string + '_2_temperatures.pdf'
    plt.savefig(filename, bbox_inches='tight')
    # plt.show()
    plt.clf()


if __name__ == '__main__':
    main()
