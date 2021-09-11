from matplotlib import pyplot as plt

def main():
    ys = []
    xs = []

    input_file = open("g_av.dat", "r")
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
    # filename = 'Ds_T=' + str(T) + '_tau=' + str(tau) +'.pdf'
    # plt.savefig(filename)
    plt.show()


if __name__ == '__main__':
    main()
