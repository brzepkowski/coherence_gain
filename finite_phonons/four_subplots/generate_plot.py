from matplotlib import pyplot as plt
from k_0_2564.plot_g_av_vs_tau import add_subplot as add_subplot_0_2564
from k_1_4102.plot_g_av_vs_tau import add_subplot as add_subplot_1_4102
from k_first_two.plot_g_av_vs_tau import add_subplot as add_subplot_first_two
from k_first_eleven.plot_g_av_vs_tau import add_subplot as add_subplot_first_eleven

def main():
    figure = plt.figure(figsize=[12, 10])
    font_size = 25 # Changes the size of all fonts in the plot
    tick_size = 25 # Changes the size of all labels on axes in the plot
    plt.rc('font', size=font_size)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
    figure, ax = add_subplot_0_2564(figure)
    figure, ax = add_subplot_1_4102(figure, ax)
    figure, ax = add_subplot_first_two(figure, ax)
    figure, ax = add_subplot_first_eleven(figure, ax)

    plt.subplots_adjust(left=0.1,
                    bottom=0.2,
                    right=0.95,
                    top=0.95,
                    wspace=0.02)
    plt.tight_layout()
    plt.savefig("../finite_phonons.pdf", bbox_inches='tight') # ../ is added, because python will think, that it's in
                                                              # the k_first_eleven catalog



if __name__ == '__main__':
    main()
