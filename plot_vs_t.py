from matplotlib import pyplot as plt
import sys, os, glob

def main():
    if len(sys.argv) != 3:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- MIN/MAX,")
        print("- T_temperature.")
        sys.exit()

    min_or_max_type_string = sys.argv[1]
    T_temp_string = sys.argv[2]

    print("T_temp: ", T_temp_string)

    xs = []
    ys = []
    labels = []

    # Get all target filenames and sort them before final plotting
    filenames = []
    os.chdir("./")
    input_file_prefix = "g_av_vs_t_" + min_or_max_type_string + "_T=" + T_temp_string + "*"
    for filename in glob.glob(input_file_prefix):
        filenames.append(filename)
    filenames.sort()

    for filename in filenames:
        print("filename: ", filename)
        input_file = open(filename, "r")

        tau_time = float(filename.split("_")[6])

        xs_partial = []
        ys_partial = []

        content = input_file.readlines()
        for line in content:
            if line != "\n":
                splitted_line = line.split()
                t_minus_tau = float(splitted_line[0]) - tau_time
                g_av = float(splitted_line[1])
                xs_partial.append(t_minus_tau)
                ys_partial.append(g_av)
            else:
                break
        input_file.close()

        xs.append(xs_partial)
        ys.append(ys_partial)
        labels.append(r"$\tau=" + str(tau_time) + "$")

    fig, ax = plt.subplots()
    for i in range(len(labels)):
        ax.plot(xs[i], ys[i], "-", label=labels[i])
    # plt.title()
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t - \tau$')
    plt.legend()
    plt.grid()
    filename = 'g_av_vs_t_' + min_or_max_type_string + '_T=' + str(T_temp_string) + '.pdf'
    plt.savefig(filename)
    # plt.show()


if __name__ == '__main__':
    main()
