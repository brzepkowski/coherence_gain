from matplotlib import pyplot as plt
from envelopes import envelopes, combine_two_envelopes
import sys, os, glob
import math
import numpy as np

def convert_to_complex(str):
    result = str.replace('(', '').replace(')', '').split(',')
    result = (float(result[0]), float(result[1]))
    result = complex(result[0], result[1])
    return result

def main():
    tau_time_0_string = "0.00"
    # tau_time_1_string = "4.00"
    t_time_string = "20.0"
    T_temp_string_0 = "34.0"
    # T_temp_string_0 = "70.0"

    all_taus = []
    all_D_pluses = []
    all_D_minuses = []
    all_D_ts = []
    all_labels = []

    # for tau_time_string in [tau_time_0_string, tau_time_1_string]:
    #     for T_temp_string in [T_temp_string_0, T_temp_string_1]:
    for tau_time_string in [tau_time_0_string]:
        for T_temp_string in [T_temp_string_0]:
            os.chdir("./")
            input_file_prefix = "g_av_vs_tau_extended_info_T=" + T_temp_string + "_" + t_time_string + "_" + tau_time_string + "*"
            input_file_postfix = ".dat"
            for filename in glob.glob(input_file_prefix):
                if input_file_postfix in filename:
                    tau_times = []
                    W_0s = []
                    W_1s = []
                    W_1s_phased = []
                    W_2s = []
                    W_2s_phased = []
                    W_3s = []
                    phases = []
                    D_pluses = []
                    D_minuses = []
                    D_ts = []
                    D_avs = []
                    p_pluses = []
                    p_minuses = []
                    g_avs = []
                    lower_bounds = []
                    upper_bounds = []

                    print("filename: ", filename)
                    input_file = open(filename, "r")

                    content = input_file.readlines()
                    for line in content:
                        if line != "\n":
                            splitted_line = line.split()
                            tau_time = float(splitted_line[0])
                            W_0 = convert_to_complex(splitted_line[1])
                            W_1 = convert_to_complex(splitted_line[2])
                            W_1_phased = convert_to_complex(splitted_line[3])
                            W_2 = convert_to_complex(splitted_line[4])
                            W_2_phased = convert_to_complex(splitted_line[5])
                            W_3 = convert_to_complex(splitted_line[6])
                            phase = convert_to_complex(splitted_line[7])
                            D_plus = float(splitted_line[8])
                            D_minus = float(splitted_line[9])
                            D_t = float(splitted_line[10])
                            p_plus = float(splitted_line[11])
                            p_minus = float(splitted_line[12])
                            g_av = float(splitted_line[13])

                            A = 0.25*(W_0 + W_3)
                            B = 0.25*(W_1_phased + W_2_phased)
                            lower_bound = 2*max(np.abs(A), np.abs(B))
                            upper_bound = 2*np.sqrt(np.abs(A)**2 + np.abs(B)**2)

                            # print("A: ", A)
                            # print("B: ", B)
                            # print("|A|: ", np.abs(A))
                            # print("|B|: ", np.abs(B))
                            # print("2|A|: ", 2*np.abs(A))
                            # print("2|B|: ", 2*np.abs(B))
                            # print("upper_bound: ", upper_bound)
                            # print("lower_bound: ", lower_bound)

                            # print("tau_time: ", tau_time)
                            # print("W_0: ", W_0)
                            # print("W_1: ", W_1)
                            # print("W_1_phased: ", W_1_phased)
                            # print("W_2: ", W_2)
                            # print("W_2_phased: ", W_2_phased)
                            # print("W_3: ", W_3)
                            # print("phase: ", phase)
                            # print("D_plus: ", D_plus)
                            # print("D_minus: ", D_minus)
                            # print("D_t: ", D_t)
                            # print("p_plus: ", p_plus)
                            # print("p_minus: ", p_minus)
                            # print("g_av: ", g_av)

                            # Normalize g_av
                            g_av /= (1 - D_t)
                            g_av *= 100

                            if len(D_ts) == 0 or (not math.isnan(D_plus) and not math.isnan(D_minus) and abs(D_t - D_ts[-1]) < 1e-3):
                                tau_times.append(tau_time)
                                W_0s.append(W_0)
                                W_1s.append(W_1)
                                W_1s_phased.append(W_1_phased)
                                W_2s.append(W_2)
                                W_2s_phased.append(W_2_phased)
                                W_3s.append(W_3)
                                phases.append(phase)
                                D_pluses.append(D_plus)
                                D_minuses.append(D_minus)
                                D_ts.append(D_t)
                                D_avs.append(p_plus*D_plus + p_minus*D_minus)
                                p_pluses.append(p_plus)
                                p_minuses.append(p_minus)
                                g_avs.append(g_av)
                                lower_bounds.append(lower_bound)
                                upper_bounds.append(upper_bound)
                        else:
                            break
                    input_file.close()

                    # # For plotting purposes get rid of the last entry in the first dataset (to remove
                    # # the rightmost entry in the xlabels on the left subplot)
                    # if tau_time_string == tau_time_0_string:
                    #     taus = taus[:-1]
                    #     D_pluses = D_pluses[:-1]
                    #     D_minuses = D_minuses[:-1]
                    #     D_ts = D_ts[:-1]

                    # all_taus.append(taus)
                    # all_D_pluses.append(D_pluses)
                    # all_D_minuses.append(D_minuses)
                    # all_D_ts.append(D_ts)
                    # # all_labels.append("T = " + T_temp_string + " K, tau = " + tau_time_string)

                    # plt.plot(tau_times, D_pluses, "-", label=r'$D_+$')
                    # plt.plot(tau_times, D_minuses, "-", label=r'$D_-$')
                    # plt.plot(tau_times, D_ts, "-", label=r'$D_t$')
                    plt.plot(tau_times, D_avs, "-", label=r'$D_{av}$')
                    plt.plot(tau_times, upper_bounds, "--", label='upper')
                    plt.plot(tau_times, lower_bounds, "--", label='lower')
                    plt.legend()
                    plt.show()
                    sys.exit()

    # fig = plt.figure(figsize=[12, 6])
    font_size = 25 # Changes the size of all fonts in the plot
    tick_size = 25 # Changes the size of all labels on axes in the plot
    plt.rc('font', size=font_size)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'


    plt.plot(tau_times, D_ts, "-", label=r'$D_t$')
    plt.plot(tau_times, D_pluses, "-", label=r'$D_+$')
    plt.plot(tau_times, D_minuses, "-", label=r'$D_-$')
    plt.plot(tau_times, upper_bounds, "--", label='upper')
    plt.plot(tau_times, lower_bounds, "--", label='lower')
    plt.legend()
    plt.show()
    sys.exit()

    plt.subplot(1, 2, 1)
    plt.plot(all_taus[0], all_D_pluses[0], '-')
    plt.plot(all_taus[0], all_D_minuses[0], '--')
    plt.plot(all_taus[1], all_D_pluses[1], '-')
    plt.plot(all_taus[1], all_D_minuses[1], '--')

    # Extract envelopes
    min_taus_D_plus_raw, min_D_plus_raw, max_taus_D_plus_raw, max_D_plus_raw = envelopes(all_taus[0], all_D_pluses[0], diff=0.01)
    min_taus_D_minus_raw, min_D_minus_raw, max_taus_D_minus_raw, max_D_minus_raw = envelopes(all_taus[0], all_D_minuses[0], diff=0.01)

    # Smooth the data to get rid of the 'blops'
    min_taus_D_plus, min_D_plus, max_taus_D_plus, max_D_plus, min_taus_D_minus, min_D_minus, max_taus_D_minus, max_D_minus = smooth_data(min_taus_D_plus_raw, min_D_plus_raw, max_taus_D_plus_raw, max_D_plus_raw, min_taus_D_minus_raw, min_D_minus_raw, max_taus_D_minus_raw, max_D_minus_raw)

    # Combine envelopes obtained for D_pluses and D_minuses
    max_taus_D, max_Ds = combine_two_envelopes(max_taus_D_plus, max_D_plus, max_taus_D_minus, max_D_minus)
    min_taus_D, min_Ds = combine_two_envelopes(min_taus_D_plus, min_D_plus, min_taus_D_minus, min_D_minus)
    # plt.plot(max_taus_D, max_Ds, '-')
    # plt.plot(min_taus_D, min_Ds, '--')
    # plt.plot(all_taus[0], all_D_ts[0], ':')
    #---------------------------------------------------------------------------

    # Extract envelopes
    min_taus_D_plus_raw, min_D_plus_raw, max_taus_D_plus_raw, max_D_plus_raw = envelopes(all_taus[1], all_D_pluses[1], diff=0.01)
    min_taus_D_minus_raw, min_D_minus_raw, max_taus_D_minus_raw, max_D_minus_raw = envelopes(all_taus[1], all_D_minuses[1], diff=0.01)

    # Smooth the data to get rid of the 'blops'
    min_taus_D_plus, min_D_plus, max_taus_D_plus, max_D_plus, min_taus_D_minus, min_D_minus, max_taus_D_minus, max_D_minus = smooth_data(min_taus_D_plus_raw, min_D_plus_raw, max_taus_D_plus_raw, max_D_plus_raw, min_taus_D_minus_raw, min_D_minus_raw, max_taus_D_minus_raw, max_D_minus_raw)

    # Combine envelopes obtained for D_pluses and D_minuses
    max_taus_D, max_Ds = combine_two_envelopes(max_taus_D_plus, max_D_plus, max_taus_D_minus, max_D_minus)
    min_taus_D, min_Ds = combine_two_envelopes(min_taus_D_plus, min_D_plus, min_taus_D_minus, min_D_minus)
    # plt.plot(max_taus_D, max_Ds, '-')
    # plt.plot(min_taus_D, min_Ds, '--')
    # plt.plot(all_taus[1], all_D_ts[1], ':')


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
    plt.subplots_adjust(left=0.1,
                    bottom=0.2,
                    right=0.95,
                    top=0.95,
                    wspace=0.02)
    filename = 'D_vs_tau_t=' + t_time_string + '_subplots.pdf'
    # plt.savefig(filename, bbox_inches='tight')
    plt.show()
    plt.clf()


if __name__ == '__main__':
    main()
