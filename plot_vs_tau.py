from matplotlib import pyplot as plt
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
    input_file_postfix = ".dat"
    for filename in glob.glob(input_file_prefix):
        if input_file_postfix in filename:
            g_avs = []
            taus = []
            D_pluses = []
            D_minuses = []
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
                    p_plus = float(splitted_line[4])
                    p_minus = float(splitted_line[5])
                    pure_phase = splitted_line[6].replace('(', '').replace(')', '').split(',')
                    pure_phase = (float(pure_phase[0]), float(pure_phase[1]))
                    pure_phase = complex(pure_phase[0], pure_phase[1])
                    # print("pure_phase: ", pure_phase)

                    # print("splitted_line: ", splitted_line)
                    taus.append(tau)
                    g_avs.append(g_av)
                    D_pluses.append(D_plus)
                    D_minuses.append(D_minus)
                    p_pluses.append(p_plus)
                    p_minuses.append(p_minus)
                    pure_phases_real.append(pure_phase.real)
                    pure_phases_imag.append(pure_phase.imag)
                else:
                    break
            input_file.close()

            if extended_mode:
                fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(6.4,10))
                axs[0].plot(taus, g_avs, "-", label=r'$g_{av}$')
                axs[1].plot(taus, p_pluses, "-", label=r'$p_+$')
                axs[1].plot(taus, p_minuses, "-", label=r'$p_-$')
                axs[2].plot(taus, D_pluses, "-", label=r'$D_+$')
                axs[2].plot(taus, D_minuses, "-", label=r'$D_-$')
                axs[3].plot(taus, pure_phases_real, "-", label=r'$Re(exp(iE\tau)/\hbar))$')
                axs[3].plot(taus, pure_phases_imag, "-", label=r'$Im(exp(iE\tau)/\hbar))$')

                # plt.title(r'$\tau\ =\ ' + str(tau) + '$')
                axs[0].grid()
                axs[0].set_ylabel(r'$g_{av}$')
                axs[1].grid()
                axs[1].legend()
                axs[1].set_ylabel(r'$probability$')
                axs[2].grid()
                axs[2].legend()
                axs[2].set_ylabel(r'$D$')
                axs[3].grid()
                axs[3].legend()
                axs[3].set_xlabel(r'$\tau$')
                filename = 'g_av_vs_tau_T=' + str(T_temp_string) + "_" + t_time_string + "_" + tau_time_string + "_" + "{:.2f}".format(taus[-1]) + '.pdf'
                plt.savefig(filename)
                # plt.show()
                plt.clf()
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
