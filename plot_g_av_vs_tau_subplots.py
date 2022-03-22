from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys, os, glob
import numpy as np
import math
from envelopes import envelopes

def smooth_data(taus_raw, g_avs_raw):
    epsilon = 7.5e-3
    g_avs = []
    taus = []
    multipliers = np.linspace(0, 1, len(taus_raw))
    # multipliers = np.linspace(0, 3, len(taus_raw))
    for i in range(len(taus_raw)):
        tau = taus_raw[i]
        g_av = g_avs_raw[i]
        multiplier = math.exp(-multipliers[i])
        if len(g_avs) == 0 or (abs(g_av - g_avs[-1]) < multiplier*epsilon):
            taus.append(tau)
            g_avs.append(g_av)

    return taus, g_avs

def smooth_data_2(taus_raw, g_avs_raw, max_index):
    g_avs = []
    taus = []
    for i in range(max_index):
        tau = taus_raw[i]
        g_av = g_avs_raw[i]
        if len(g_avs) == 0 or g_av - g_avs[-1] < 0:
            taus.append(tau)
            g_avs.append(g_av)
    g_avs = g_avs + g_avs_raw[max_index:]
    taus = taus + taus_raw[max_index:]
    return taus, g_avs

def main():
    tau_time_0_string = "0.00"
    tau_time_1_string = "4.00"
    t_time_string = "20.0"
    T_temp_string_1 = "34.0"
    T_temp_string_2 = "70.0"

    all_taus = []
    all_g_avs = []
    all_labels = []

    for tau_time_string in [tau_time_0_string, tau_time_1_string]:
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

                            # Normalize g_av
                            g_av /= (1 - D_t_minus_tau)
                            g_av *= 100

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

    fig = plt.figure(figsize=[10, 12])
    plt.rc('font', size=25)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    ax = plt.subplot(2, 1, 1)
    # plt.plot(all_taus[0], all_g_avs[0], '-', label=all_labels[0], color='C3')
    # plt.plot(all_taus[1], all_g_avs[1], '-', label=all_labels[1], color='C4')

    # Extract envelopes
    min_taus_g_av_0, min_g_avs_0, max_taus_g_av_0, max_g_avs_0 = envelopes(all_taus[0], all_g_avs[0], epsilon_min=1e-4, epsilon_max=1e-8, diff=0.015)
    min_taus_g_av_0, min_g_avs_0 = smooth_data(min_taus_g_av_0, min_g_avs_0)
    min_taus_g_av_0, min_g_avs_0 = smooth_data_2(min_taus_g_av_0, min_g_avs_0, 100)
    max_taus_g_av_0, max_g_avs_0 = smooth_data_2(max_taus_g_av_0, max_g_avs_0, 30)
    plt.plot(max_taus_g_av_0, max_g_avs_0, '-', color='#1f77b4')
    plt.plot(min_taus_g_av_0, min_g_avs_0, '-', color='#1f77b4')

    # Extract envelopes
    min_taus_g_av_1, min_g_avs_1, max_taus_g_av_1, max_g_avs_1 = envelopes(all_taus[1], all_g_avs[1], epsilon_min=1e-4, epsilon_max=1e-8, diff=0.015)
    min_taus_g_av_1, min_g_avs_1 = smooth_data_2(min_taus_g_av_1, min_g_avs_1, 20)
    plt.plot(max_taus_g_av_1, max_g_avs_1, '--', color='#ff7f0e')
    plt.plot(min_taus_g_av_1, min_g_avs_1, '--', color='#ff7f0e')

    # plt.annotate('34 K', xy =(0.85, 5.75))
    # plt.annotate('70 K', xy =(1.45, 4.2))
    plt.annotate('34 K', xy =(1.15, 6.55))
    plt.annotate('70 K', xy =(1.95, 5.85))
    # plt.hlines(0, 0.0, 0.2, color='black', linestyle='dashed')
    plt.ylabel(r'$g_{av}^\prime\ [\%]$')
    plt.xlim([min(all_taus[0]), max(all_taus[0])])
    plt.locator_params(axis='x', nbins=5)
    plt.tick_params(axis='x', pad=15)
    # plt.ylim([-0.1, 0.1])
    plt.locator_params(axis='y', nbins=5)

    # Add inset
    axins = inset_axes(ax, width="50%", height="35%", loc="center right", bbox_to_anchor=(-0.02,-0.12,1,1), bbox_transform=ax.transAxes)
    # axins.plot(all_taus[0], all_g_avs[0], '-', label=all_labels[0], color='C1')
    axins.plot(min_taus_g_av_0, min_g_avs_0, '-', color='C0')
    axins.plot(max_taus_g_av_0, max_g_avs_0, '-', color='C0')
    axins.plot(all_taus[0][:3136], all_g_avs[0][:3136], '-', color='C0') # Fill the initial gap by exact values
    # axins.plot(all_taus[1], all_g_avs[1], '-', label=all_labels[1], color='C0')
    axins.plot(min_taus_g_av_1, min_g_avs_1, '--', color='C1')
    axins.plot(max_taus_g_av_1, max_g_avs_1, '--', color='C1')
    axins.plot(all_taus[1][:3136], all_g_avs[1][:3136], '--', color='C1') # Fill the initial gap by exact values
    axins.set_ylim([-0.15, 0])
    axins.set_xlim([0, 1])

    plt.subplot(2, 1, 2)
    plt.plot(all_taus[2], all_g_avs[2], '-', label=all_labels[2])
    plt.plot(all_taus[3], all_g_avs[3], '--', label=all_labels[3])
    plt.annotate('34 K', xy =(4.007, 7))
    plt.annotate('70 K', xy =(4.0078, 1))
    plt.ylabel(r'$g_{av}^\prime\ [\%]$')
    plt.xlabel(r'$\tau\ [ps]$')
    plt.xlim([min(all_taus[2]), max(all_taus[2])])
    plt.locator_params(axis='x', nbins=5)
    plt.tick_params(axis='x', pad=15)
    plt.ylim([0, 8])
    plt.locator_params(axis='y', nbins=5)

    plt.tight_layout()
    filename = 'g_av_vs_tau_t=' + t_time_string + '_subplots.pdf'
    plt.savefig(filename, bbox_inches='tight')
    # plt.show()
    # plt.clf()


if __name__ == '__main__':
    main()
