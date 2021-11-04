from matplotlib import pyplot as plt
import sys, os, glob

def main():
    if len(sys.argv) != 4:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- tau_time_start, ")
        print("- t_time, ")
        print("- T_temperature.")
        sys.exit()

    tau_time_string = sys.argv[1]
    t_time_string = sys.argv[2]
    T_temp_string = sys.argv[3]

    print("tau_time: ", tau_time_string)
    print("t_time: ", t_time_string)
    print("T_temp: ", T_temp_string)

    ys = []
    xs = []

    os.chdir("./")
    input_file_prefix = "g_av_vs_tau_T=" + T_temp_string + "_" + t_time_string + "_" + tau_time_string + "*"
    for filename in glob.glob(input_file_prefix):
        print("filename: ", filename)
        input_file = open(filename, "r")

        content = input_file.readlines()
        for line in content:
            if line != "\n":
                splitted_line = line.split()
                tau = float(splitted_line[0])
                g_av = float(splitted_line[1])
                xs.append(tau)
                ys.append(g_av)
            else:
                break
        input_file.close()

        plt.plot(xs, ys, "-")
        # plt.title(r'$\tau\ =\ ' + str(tau) + '$')
        plt.ylabel(r'$g_{av}$')
        plt.xlabel(r'$\tau$')
        # plt.legend()
        plt.grid()
        filename = 'g_av_vs_tau_T=' + str(T_temp_string) + "_" + t_time_string + "_" + tau_time_string + "_" + "{:.2f}".format(xs[-1]) + '.pdf'
        plt.savefig(filename)
        # plt.show()


if __name__ == '__main__':
    main()
