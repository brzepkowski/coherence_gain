from matplotlib import pyplot as plt
import numpy as np
import sys, os, glob

def main():
    if len(sys.argv) != 4:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- tau_MIN,")
        print("- tau_MAX,")
        print("- T_temperature.")
        sys.exit()

    tau_min_string = sys.argv[1]
    tau_max_string = sys.argv[2]
    T_temp_string = sys.argv[3]
    tau_min = float(tau_min_string)
    tau_max = float(tau_max_string)

    print("tau_min_string: ", tau_min_string)
    print("tau_max_string: ", tau_max_string)
    print("T_temp: ", T_temp_string)

    # TO DO REST

    xs = []
    g_avs = []
    labels = []

    filename_postfix = ".dat"
    os.chdir("./")

    ############################################################################
    # tau_min
    ############################################################################
    filename_tau_min_prefix = "g_av_vs_t_MIN_T=" + T_temp_string + "_" + tau_min_string + "_" + tau_min_string + "_*"

    print("filename_tau_min_prefix: ", filename_tau_min_prefix)

    # Get all target filenames and sort them before final plotting
    for filename in glob.glob(filename_tau_min_prefix):
        if filename_postfix in filename:
            print("filename: ", filename)
            input_file = open(filename, "r")

            xs_partial = []
            g_avs_partial = []

            content = input_file.readlines()
            for line in content:
                if line != "\n":
                    splitted_line = line.split()
                    t_minus_tau = float(splitted_line[0]) - tau_min
                    g_av = float(splitted_line[1])
                    xs_partial.append(t_minus_tau)
                    g_avs_partial.append(g_av)
                else:
                    break
            input_file.close()
            xs.append(xs_partial)
            g_avs.append(g_avs_partial)
            labels.append(r"$\tau_{MIN}=" + tau_min_string + "$")
            break

    xs_inset = xs[0][:int(0.0015*len(xs[0]))]
    g_avs_inset = g_avs[0][:int(0.0015*len(g_avs[0]))]

    fig, ax = plt.subplots()
    ax.plot(xs[0], g_avs[0], "-", label=labels[0])

    axin1 = ax.inset_axes([0.6, 0.5, 0.35, 0.35])
    axin1.plot(xs_inset, g_avs_inset, "-")

    # Plot line corresponding to t = 2 * tau
    ax.vlines(x = tau_min, ymin = np.min(g_avs[0]), ymax = np.max(g_avs[0]),
           colors = 'purple',
           linestyles = 'dashed',
           label = r'$t = 2\tau_{MIN}$')

    # plt.title()
    ax.set_ylabel(r'$g_{av}$')
    ax.set_xlabel(r'$t - \tau$')
    ax.legend()
    ax.grid()
    filename_plot = 'g_av_vs_t_MIN_T=' + T_temp_string + "_" + tau_min_string + ".pdf"
    plt.savefig(filename_plot)
    plt.clf()

    ############################################################################
    # tau_max
    ############################################################################
    filename_tau_max_prefix = "g_av_vs_t_MAX_T=" + T_temp_string + "_" + tau_max_string + "_" + tau_max_string + "_*"

    print("filename_tau_max_prefix: ", filename_tau_max_prefix)

    # Get all target filenames and sort them before final plotting
    for filename in glob.glob(filename_tau_max_prefix):
        if filename_postfix in filename:
            print("filename: ", filename)
            input_file = open(filename, "r")

            xs_partial = []
            g_avs_partial = []

            content = input_file.readlines()
            for line in content:
                if line != "\n":
                    splitted_line = line.split()
                    t_minus_tau = float(splitted_line[0]) - tau_max
                    g_av = float(splitted_line[1])
                    xs_partial.append(t_minus_tau)
                    g_avs_partial.append(g_av)
                else:
                    break
            input_file.close()
            xs.append(xs_partial)
            g_avs.append(g_avs_partial)
            labels.append(r"$\tau_{MAX}=" + tau_max_string + "$")
            break

    xs_inset = xs[1][:int(0.0015*len(xs[1]))]
    g_avs_inset = g_avs[1][:int(0.0015*len(g_avs[1]))]

    fig, ax = plt.subplots()
    ax.plot(xs[1], g_avs[1], "-", label=labels[1])

    axin1 = ax.inset_axes([0.6, 0.5, 0.35, 0.35])
    axin1.plot(xs_inset, g_avs_inset, "-")

    # Plot line corresponding to t = 2 * tau
    ax.vlines(x = tau_max, ymin = np.min(g_avs[1]), ymax = np.max(g_avs[1]),
           colors = 'purple',
           linestyles = 'dashed',
           label = r'$t = 2\tau_{MAX}$')

    # plt.title()
    ax.set_ylabel(r'$g_{av}$')
    ax.set_xlabel(r'$t - \tau$')
    ax.legend()
    ax.grid()
    filename_plot = 'g_av_vs_t_MAX_T=' + T_temp_string + "_" + tau_max_string + ".pdf"
    plt.savefig(filename_plot)
    plt.clf()

    ############################################################################
    # Plot both g_avs at the same plot
    ############################################################################

    plt.plot(xs[0], g_avs[0], "-", label=labels[0])
    plt.plot(xs[1], g_avs[1], "-", label=labels[1])
    # plt.title()
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t - \tau$')
    plt.legend()
    plt.grid()
    filename_plot = 'g_av_vs_t_MIN_MAX_T=' + T_temp_string + "_" + tau_min_string + "_" + tau_max_string + ".pdf"
    plt.savefig(filename_plot)
    # plt.show()


if __name__ == '__main__':
    main()
