from matplotlib import pyplot as plt
from envelopes import envelopes, combine_two_envelopes
import sys, os, glob
import math
import numpy as np

def smooth_data(min_taus_D_plus_raw, min_D_plus_raw, max_taus_D_plus_raw, max_D_plus_raw, min_taus_D_minus_raw, min_D_minus_raw, max_taus_D_minus_raw, max_D_minus_raw):
    epsilon = 1e-1
    min_D_plus = []
    min_taus_D_plus = []
    multipliers = np.linspace(0, 8, len(min_D_plus_raw))
    for i in range(len(min_D_plus_raw)):
        D_plus = min_D_plus_raw[i]
        multiplier = math.exp(-multipliers[i])
        if len(min_D_plus) == 0 or ((abs(D_plus - min_D_plus[-1])/min_D_plus[-1]) < multiplier*epsilon):
            min_D_plus.append(D_plus)
            min_taus_D_plus.append(min_taus_D_plus_raw[i])

    epsilon = 1e-4
    min_D_minus = []
    min_taus_D_minus = []
    multipliers = np.linspace(0, 4, len(min_D_minus_raw))
    for i in range(len(min_D_minus_raw)):
        D_minus = min_D_minus_raw[i]
        multiplier = math.exp(-multipliers[i])
        if len(min_D_minus) == 0 or ((abs(D_minus - min_D_minus[-1])/ min_D_minus[-1]) < multiplier*epsilon):
            min_D_minus.append(D_minus)
            min_taus_D_minus.append(min_taus_D_minus_raw[i])

    epsilon = 1e-3
    max_D_plus = []
    max_taus_D_plus = []
    for i in range(len(max_D_plus_raw)):
        D_plus = max_D_plus_raw[i]
        if len(max_D_plus) == 0 or ((abs(D_plus - max_D_plus[-1])/max_D_plus[-1]) < epsilon):
            max_D_plus.append(D_plus)
            max_taus_D_plus.append(max_taus_D_plus_raw[i])

    max_D_minus = []
    max_taus_D_minus = []
    for i in range(len(max_D_minus_raw)):
        D_minus = max_D_minus_raw[i]
        if len(max_D_minus) == 0 or ((abs(D_minus - max_D_minus[-1])/ max_D_minus[-1]) < epsilon):
            max_D_minus.append(D_minus)
            max_taus_D_minus.append(max_taus_D_minus_raw[i])

    return min_taus_D_plus, min_D_plus, max_taus_D_plus, max_D_plus, min_taus_D_minus, min_D_minus, max_taus_D_minus, max_D_minus

def main():
    tau_time_0_string = "0.00"
    tau_time_1_string = "4.00"
    t_time_string = "20.0"
    T_temp_string_0 = "34.0"
    T_temp_string_1 = "70.0"

    all_taus = []
    all_D_pluses = []
    all_D_minuses = []
    all_D_ts = []
    all_labels = []

    for tau_time_string in [tau_time_0_string, tau_time_1_string]:
        for T_temp_string in [T_temp_string_0, T_temp_string_1]:
            os.chdir("./")
            input_file_prefix = "g_av_vs_tau_T=" + T_temp_string + "_" + t_time_string + "_" + tau_time_string + "*"
            input_file_postfix = ".dat"
            for filename in glob.glob(input_file_prefix):
                if input_file_postfix in filename:
                    g_avs = []
                    taus = []
                    D_pluses = []
                    D_minuses = []
                    D_ts = []
                    p_pluses = []
                    p_minuses = []
                    pure_phases_real = []
                    pure_phases_imag = []

                    print("filename: ", filename)
                    input_file = open(filename, "r")

                    content = input_file.readlines()
                    for line in content:
                        if line != "\n":
                            splitted_line = line.split()
                            tau = float(splitted_line[0])
                            g_av = float(splitted_line[1])
                            D_plus = float(splitted_line[2])
                            D_minus = float(splitted_line[3])
                            D_t = float(splitted_line[4])
                            p_plus = float(splitted_line[5])
                            p_minus = float(splitted_line[6])
                            pure_phase = splitted_line[7].replace('(', '').replace(')', '').split(',')
                            pure_phase = (float(pure_phase[0]), float(pure_phase[1]))
                            pure_phase = complex(pure_phase[0], pure_phase[1])
                            # print("pure_phase: ", pure_phase)

                            # Normalize g_av
                            g_av /= (1 - D_t)
                            g_av *= 100

                            if len(D_ts) == 0 or (not math.isnan(D_plus) and not math.isnan(D_minus) and abs(D_t - D_ts[-1]) < 1e-3):
                                taus.append(tau)
                                g_avs.append(g_av)
                                D_pluses.append(D_plus)
                                D_minuses.append(D_minus)
                                D_ts.append(D_t)
                                p_pluses.append(p_plus)
                                p_minuses.append(p_minus)
                                pure_phases_real.append(pure_phase.real)
                                pure_phases_imag.append(pure_phase.imag)
                        else:
                            break
                    input_file.close()

                    # For plotting purposes get rid of the last entry in the first dataset (to remove
                    # the rightmost entry in the xlabels on the left subplot)
                    if tau_time_string == tau_time_0_string:
                        taus = taus[:-1]
                        D_pluses = D_pluses[:-1]
                        D_minuses = D_minuses[:-1]
                        D_ts = D_ts[:-1]

                    all_taus.append(taus)
                    all_D_pluses.append(D_pluses)
                    all_D_minuses.append(D_minuses)
                    all_D_ts.append(D_ts)
                    # all_labels.append("T = " + T_temp_string + " K, tau = " + tau_time_string)

    fig = plt.figure(figsize=[12, 7.12])
    font_size = 35 # Changes the size of all fonts in the plot
    tick_size = 35 # Changes the size of all labels on axes in the plot
    plt.rc('font', size=font_size)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    plt.subplot(1, 2, 1)
    plt.plot(all_taus[0], all_D_pluses[0], '-', color='black')
    plt.plot(all_taus[0], all_D_minuses[0], '--', color='black')
    plt.plot(all_taus[1], all_D_pluses[1], '-', color='black')
    plt.plot(all_taus[1], all_D_minuses[1], '--', color='black')

    # Extract envelopes
    min_taus_D_plus_raw, min_D_plus_raw, _, _ = envelopes(all_taus[0], all_D_pluses[0], epsilon_min=1e-8, epsilon_max=1e-5, diff=0.005)
    _, _, max_taus_D_plus_raw, max_D_plus_raw = envelopes(all_taus[0], all_D_pluses[0], diff=0.1)

    min_taus_D_minus_raw, min_D_minus_raw, _, _ = envelopes(all_taus[0], all_D_minuses[0], epsilon_min=1e-8, epsilon_max=1e-5, diff=0.005)
    _, _, max_taus_D_minus_raw, max_D_minus_raw = envelopes(all_taus[0], all_D_minuses[0], diff=0.1)

    # Smooth the data to get rid of the 'blops'
    min_taus_D_plus, min_D_plus, max_taus_D_plus, max_D_plus, min_taus_D_minus, min_D_minus, max_taus_D_minus, max_D_minus = smooth_data(min_taus_D_plus_raw, min_D_plus_raw, max_taus_D_plus_raw, max_D_plus_raw, min_taus_D_minus_raw, min_D_minus_raw, max_taus_D_minus_raw, max_D_minus_raw)

    # Combine envelopes obtained for D_pluses and D_minuses
    max_taus_D, max_Ds = combine_two_envelopes(max_taus_D_plus, max_D_plus, max_taus_D_minus, max_D_minus)
    min_taus_D, min_Ds = combine_two_envelopes(min_taus_D_plus, min_D_plus, min_taus_D_minus, min_D_minus)
    plt.plot(max_taus_D, max_Ds, '-', color="tab:blue")
    plt.plot(min_taus_D, min_Ds, '--', color="tab:orange")
    plt.plot(all_taus[0], all_D_ts[0], ':', color="tab:green")
    #---------------------------------------------------------------------------

    # Extract envelopes
    min_taus_D_plus_raw, min_D_plus_raw, max_taus_D_plus_raw, max_D_plus_raw = envelopes(all_taus[1], all_D_pluses[1], diff=0.1)
    min_taus_D_minus_raw, min_D_minus_raw, max_taus_D_minus_raw, max_D_minus_raw = envelopes(all_taus[1], all_D_minuses[1], diff=0.1)

    # Smooth the data to get rid of the 'blops'
    min_taus_D_plus, min_D_plus, max_taus_D_plus, max_D_plus, min_taus_D_minus, min_D_minus, max_taus_D_minus, max_D_minus = smooth_data(min_taus_D_plus_raw, min_D_plus_raw, max_taus_D_plus_raw, max_D_plus_raw, min_taus_D_minus_raw, min_D_minus_raw, max_taus_D_minus_raw, max_D_minus_raw)

    # Combine envelopes obtained for D_pluses and D_minuses
    max_taus_D, max_Ds = combine_two_envelopes(max_taus_D_plus, max_D_plus, max_taus_D_minus, max_D_minus)
    min_taus_D, min_Ds = combine_two_envelopes(min_taus_D_plus, min_D_plus, min_taus_D_minus, min_D_minus)
    plt.plot(max_taus_D, max_Ds, '-', color="tab:red")
    plt.plot(min_taus_D, min_Ds, '--', color="tab:purple")
    plt.plot(all_taus[1], all_D_ts[1], ':', color="tab:brown")


    plt.annotate(r'$34\ K$', xy =(all_taus[0][math.floor(0.75*len(all_taus[0]))], 0.7))
    plt.annotate(r'$70\ K$', xy =(all_taus[0][math.floor(0.75*len(all_taus[0]))], 0.2))
    plt.xlabel(r'$\tau\ [ps]$', fontsize=font_size)
    plt.ylabel(r'$D$', fontsize=font_size)
    plt.xlim([min(all_taus[0]), max(all_taus[0])])
    # plt.locator_params(axis='x', nbins=5)
    plt.tick_params(axis='x', pad=15, labelsize=tick_size)
    plt.tick_params(axis='y', labelsize=tick_size)
    plt.ylim([0, 1])
    # plt.locator_params(axis='y', nbins=5)

    ############################################################################

    plt.subplot(1, 2, 2)
    plt.plot(all_taus[2], all_D_pluses[2], '-')
    plt.plot(all_taus[2], all_D_minuses[2], '--')
    plt.plot(all_taus[2], all_D_ts[2], ':')
    plt.plot(all_taus[3], all_D_pluses[3], '-')
    plt.plot(all_taus[3], all_D_minuses[3], '--')
    plt.plot(all_taus[3], all_D_ts[3], ':')

    plt.annotate(r'$34\ K$', xy =(all_taus[2][math.floor(0.75*len(all_taus[2]))], 0.7))
    plt.annotate(r'$70\ K$', xy =(all_taus[2][math.floor(0.75*len(all_taus[2]))], 0.2))
    plt.xlabel(r'$\tau\ [ps]$', fontsize=font_size)
    plt.xlim([min(all_taus[2]), max(all_taus[2])])
    # plt.locator_params(axis='x', nbins=5)
    plt.tick_params(axis='x', pad=15, labelsize=tick_size)
    plt.ylim([0, 1])
    plt.yticks([])
    # plt.locator_params(axis='y', nbins=5)

    # plt.tight_layout()
    plt.subplots_adjust(left=0.11,
                    bottom=0.18,
                    right=0.98,
                    top=0.95,
                    wspace=0.02)
    filename = 'D_vs_tau_t=' + t_time_string + '_subplots.pdf'
    plt.savefig(filename, bbox_inches='tight')
    # plt.show()
    plt.clf()


if __name__ == '__main__':
    main()
