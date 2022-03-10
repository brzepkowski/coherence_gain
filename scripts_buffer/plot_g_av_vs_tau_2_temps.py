from matplotlib import pyplot as plt
import sys, os, glob

def main():
    if len(sys.argv) != 5:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- tau_time_start, ")
        print("- t_time, ")
        print("- T_temperature 1,")
        print("- T_temperature 2.")
        sys.exit()

    tau_time_string = sys.argv[1]
    t_time_string = sys.argv[2]
    T_temp_string_1 = sys.argv[3]
    T_temp_string_2 = sys.argv[4]

    print("tau_time: ", tau_time_string)
    print("t_time: ", t_time_string)
    print("T_temp_1: ", T_temp_string_1)
    print("T_temp_2: ", T_temp_string_2)

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

                        # print("splitted_line: ", splitted_line)
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

                all_g_avs.append(g_avs)
                all_labels.append("T = " + T_temp_string + " K")


            plt.figure(figsize=[10, 7.5])
            plt.rc('font', size=30)
            for i in range(len(all_g_avs)):
                # plt.plot(taus, all_g_avs[i], "-", label=all_labels[i])
                plt.plot(taus, all_g_avs[i], "-", color='black', label=all_labels[i])
            plt.ylabel(r'$g_{av}^\prime\ [\%]$')
            plt.xlabel(r'$\tau$')
            # plt.annotate('34 K', xy =(0.155, -0.085))
            # plt.annotate('70 K', xy =(0.155, 0.055))

            plt.annotate('34 K', xy =(4.0069, 7))
            plt.annotate('70 K', xy =(4.0078, 1))

            plt.xlim([4, 4.01])
            plt.locator_params(axis='x', nbins=5)
            plt.tick_params(axis='x', pad=15)
            plt.ylim([0, 8])
            plt.locator_params(axis='y', nbins=5)
            # plt.hlines(0, 0.0, 0.2, color='black', linestyle='dashed')
            # plt.legend()
            filename = 'g_av_vs_tau_T1=' + T_temp_string_1 + "_T2=" + T_temp_string_2 + "_" + t_time_string + "_" + tau_time_string + "_" + "{:.2f}".format(taus[-1]) + '.pdf'
            plt.savefig(filename, bbox_inches='tight')
            # plt.show()
            plt.clf()


if __name__ == '__main__':
    main()
