from matplotlib import pyplot as plt
import sys, os, glob
from .envelopes import envelopes

def add_subplot(axes):

    tau_time_string = "0.00"
    t_time_string = "20.0"
    T_temp_string = "34.0"

    print("tau_time: ", tau_time_string)
    print("t_time: ", t_time_string)
    print("T_temp: ", T_temp_string)

    os.chdir("../k_first_two/") # python 'remained' in the previous location ("k_0_2564" catalog), so we need to switch it
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
                    if abs(D_t - 1) > 1e-5:
                        g_av /= (1 - D_t)
                        g_av *= 100
                    else:
                        break

                    # print("splitted_line: ", splitted_line)
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

            # Extract envelopes
            min_taus_g_av_0, min_g_avs_0, max_taus_g_av_0, max_g_avs_0 = envelopes(taus, g_avs, epsilon_min=1e-4, epsilon_max=1e-3, diff=0.41)

            plt.sca(axes[1, 0])

            # Plot envelopes
            plt.plot(min_taus_g_av_0, min_g_avs_0, "-", color='#1f77b4')
            plt.plot(max_taus_g_av_0, max_g_avs_0, "-", color='#1f77b4')
            plt.ylabel(r'$g_{av}^\prime\ [\%]$', font="stix", fontsize=35)
            plt.xlabel(r'$\tau\ [ps]$', font="stix", fontsize=35)
            plt.xlim(0, 10)
            plt.xticks([0,1,2,3,4,5,6,7,8,9], font="stix", fontsize=35)
            plt.ylim(-20.0, 35.0)
            plt.yticks([-20, -10, 0, 10, 20, 30], font="stix", fontsize=35)
            # plt.grid()

            plt.axhline(y=0.0, color='black', linestyle='--')

            plt.annotate(r'$k \in \{k_0, k_1 \}$', (9.9, 30.0), horizontalalignment='right')
            # plt.annotate(r'$c)$', (9.6, -18.0), horizontalalignment='right')
            axes[1, 0].text(
                9.7, -16.5, r'$c)$', ha="right", va="center", size=35,
                bbox=dict(boxstyle="square,pad=0.1", facecolor='none'))

            return axes

if __name__ == '__main__':
    main()
