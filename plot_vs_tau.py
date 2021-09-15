from matplotlib import pyplot as plt
import sys

def main():
    if len(sys.argv) != 2:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- the temperature 'T' [K].")
        sys.exit()

    T_string = sys.argv[1]
    T = float(T_string)

    ys = []
    xs = []

    input_file = open("g_av_T=" + T_string + ".dat", "r")
    content = input_file.readlines()
    for line in content:
        splitted_line = line.split()
        tau = float(splitted_line[0])
        g_av = float(splitted_line[1])
        xs.append(tau)
        ys.append(g_av)
    input_file.close()

    plt.plot(xs, ys, "-")
    # plt.title(r'$\tau\ =\ ' + str(tau) + '$')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    # plt.legend()
    plt.grid()
    filename = 'g_av_T=' + str(T) +'.pdf'
    plt.savefig(filename)
    plt.show()


if __name__ == '__main__':
    main()
