import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob, os
import sys

def main():
    if len(sys.argv) != 3:
        print("Wrong number of parameters! Please provide them in the following way:")
        print("- the temperature 'T' [K],")
        print("- time 'tau' [ps].")
        sys.exit()

    T = float(sys.argv[1])
    tau = float(sys.argv[2])

    results_D_plus = []
    results_D_minus = []
    results_D_t = []
    results_p_plus = []
    results_p_minus = []
    results_g_av = []

    results_W_0_real = []
    results_W_1_real = []
    results_W_2_real = []
    results_W_3_real = []

    results_W_0_imag = []
    results_W_1_imag = []
    results_W_2_imag = []
    results_W_3_imag = []

    results_phase_W_0_real = []
    results_phase_W_1_real = []
    results_phase_W_2_real = []
    results_phase_W_3_real = []

    results_phase_W_0_imag = []
    results_phase_W_1_imag = []
    results_phase_W_2_imag = []
    results_phase_W_3_imag = []

    results_pure_phase_real = []
    results_pure_phase_imag = []

    os.chdir("./")
    for filename in glob.glob("data_*"):
        if '_T=' + str(T) + '_' in filename and '_tau=' + str(tau) + '_' in filename:
            with open(filename, 'rb') as handle:
                data = pickle.load(handle)
                results_D_plus += data['results_D_plus']
                results_D_minus += data['results_D_minus']
                results_D_t += data['results_D_t']
                results_p_plus += data['results_p_plus']
                results_p_minus += data['results_p_minus']
                results_g_av += data['results_g_av']

                results_W_0_real += data['results_W_0_real']
                results_W_1_real += data['results_W_1_real']
                results_W_2_real += data['results_W_2_real']
                results_W_3_real += data['results_W_3_real']

                results_W_0_imag += data['results_W_0_imag']
                results_W_1_imag += data['results_W_1_imag']
                results_W_2_imag += data['results_W_2_imag']
                results_W_3_imag += data['results_W_3_imag']

                results_phase_W_0_real += data['results_phase_W_0_real']
                results_phase_W_1_real += data['results_phase_W_1_real']
                results_phase_W_2_real += data['results_phase_W_2_real']
                results_phase_W_3_real += data['results_phase_W_3_real']

                results_phase_W_0_imag += data['results_phase_W_0_imag']
                results_phase_W_1_imag += data['results_phase_W_1_imag']
                results_phase_W_2_imag += data['results_phase_W_2_imag']
                results_phase_W_3_imag += data['results_phase_W_3_imag']

                results_pure_phase_real += data['results_pure_phase_real']
                results_pure_phase_imag += data['results_pure_phase_imag']

    ############################################################################
    # Sort all of the lists
    ############################################################################
    results_D_plus = sorted(results_D_plus, key=lambda x: (x[0], x[1]))
    results_D_minus = sorted(results_D_minus, key=lambda x: (x[0], x[1]))
    results_D_t = sorted(results_D_t, key=lambda x: (x[0], x[1]))
    results_p_plus = sorted(results_p_plus, key=lambda x: (x[0], x[1]))
    results_p_minus = sorted(results_p_minus, key=lambda x: (x[0], x[1]))
    results_g_av = sorted(results_g_av, key=lambda x: (x[0], x[1]))

    results_W_0_real = sorted(results_W_0_real, key=lambda x: (x[0], x[1]))
    results_W_1_real = sorted(results_W_1_real, key=lambda x: (x[0], x[1]))
    results_W_2_real = sorted(results_W_2_real, key=lambda x: (x[0], x[1]))
    results_W_3_real = sorted(results_W_3_real, key=lambda x: (x[0], x[1]))

    results_W_0_imag = sorted(results_W_0_imag, key=lambda x: (x[0], x[1]))
    results_W_1_imag = sorted(results_W_1_imag, key=lambda x: (x[0], x[1]))
    results_W_2_imag = sorted(results_W_2_imag, key=lambda x: (x[0], x[1]))
    results_W_3_imag = sorted(results_W_3_imag, key=lambda x: (x[0], x[1]))

    results_phase_W_0_real = sorted(results_phase_W_0_real, key=lambda x: (x[0], x[1]))
    results_phase_W_1_real = sorted(results_phase_W_1_real, key=lambda x: (x[0], x[1]))
    results_phase_W_2_real = sorted(results_phase_W_2_real, key=lambda x: (x[0], x[1]))
    results_phase_W_3_real = sorted(results_phase_W_3_real, key=lambda x: (x[0], x[1]))

    results_phase_W_0_imag = sorted(results_phase_W_0_imag, key=lambda x: (x[0], x[1]))
    results_phase_W_1_imag = sorted(results_phase_W_1_imag, key=lambda x: (x[0], x[1]))
    results_phase_W_2_imag = sorted(results_phase_W_2_imag, key=lambda x: (x[0], x[1]))
    results_phase_W_3_imag = sorted(results_phase_W_3_imag, key=lambda x: (x[0], x[1]))

    results_pure_phase_real = sorted(results_pure_phase_real, key=lambda x: (x[0], x[1]))
    results_pure_phase_imag = sorted(results_pure_phase_imag, key=lambda x: (x[0], x[1]))

    ############################################################################
    # Get all ts
    ############################################################################
    ts = [entry[1] for entry in results_D_plus]

    ############################################################################
    # Cast all data to an numpy array
    ############################################################################
    results_D_plus = np.array(results_D_plus)
    results_D_minus = np.array(results_D_minus)
    results_D_t = np.array(results_D_t)
    results_p_plus = np.array(results_p_plus)
    results_p_minus = np.array(results_p_minus)
    results_g_av = np.array(results_g_av)

    results_W_0_real = np.array(results_W_0_real)
    results_W_1_real = np.array(results_W_1_real)
    results_W_2_real = np.array(results_W_2_real)
    results_W_3_real = np.array(results_W_3_real)

    results_W_0_imag = np.array(results_W_0_imag)
    results_W_1_imag = np.array(results_W_1_imag)
    results_W_2_imag = np.array(results_W_2_imag)
    results_W_3_imag = np.array(results_W_3_imag)

    results_phase_W_0_real = np.array(results_phase_W_0_real)
    results_phase_W_1_real = np.array(results_phase_W_1_real)
    results_phase_W_2_real = np.array(results_phase_W_2_real)
    results_phase_W_3_real = np.array(results_phase_W_3_real)

    results_phase_W_0_imag = np.array(results_phase_W_0_imag)
    results_phase_W_1_imag = np.array(results_phase_W_1_imag)
    results_phase_W_2_imag = np.array(results_phase_W_2_imag)
    results_phase_W_3_imag = np.array(results_phase_W_3_imag)

    results_pure_phase_real = np.array(results_pure_phase_real)
    results_pure_phase_imag = np.array(results_pure_phase_imag)

    ############################################################################
    # Generate plots
    ############################################################################

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_D_plus[:,-1], "-", label=r'$D_+$')
    plt.plot(ts, results_D_minus[:,-1], "-", label=r'$D_-$')
    plt.plot(ts, results_D_t[:,-1], "-", label=r'$D_t$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '$')
    plt.ylabel(r'$D$')
    plt.xlabel(r'$t$')
    plt.legend()
    plt.grid()
    filename = 'Ds_T=' + str(T) + '_tau=' + str(tau) +'.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_p_plus[:,-1], "-", label=r'$p_+$')
    plt.plot(ts, results_p_minus[:,-1], "-", label=r'$p_-$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '$')
    plt.ylabel(r'$p$')
    plt.xlabel(r'$t$')
    plt.legend()
    plt.grid()
    filename = 'probabilities_T=' + str(T) + '_tau=' + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_g_av[:,-1], "-")
    plt.title(r'$\tau\ =\ ' + str(tau) + '$')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t$')
    plt.grid()
    filename = 'g_av_T=' + str(T) + '_tau=' + str(tau) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_phase_W_0_real[:,-1], "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_phase_W_1_real[:,-1], "-", label=r'$e^{-iE\tau/\hbar} \langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_phase_W_2_real[:,-1], "-", label=r'$e^{iE\tau/\hbar} \langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_phase_W_3_real[:,-1], "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_pure_phase_real[:,-1], "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (real)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t$')
    plt.grid()
    plt.legend()
    filename = 'components_real_T=' + str(T) + '_tau=' + str(tau) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_phase_W_0_imag[:,-1], "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_phase_W_1_imag[:,-1], "-", label=r'$e^{-iE\tau/\hbar} \langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_phase_W_2_imag[:,-1], "-", label=r'$e^{iE\tau/\hbar} \langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_phase_W_3_imag[:,-1], "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_pure_phase_imag[:,-1], "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (imag)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t$')
    plt.grid()
    plt.legend()
    filename = 'components_imag_T=' + str(T) + '_tau=' + str(tau) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_W_0_real[:,-1], "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_W_1_real[:,-1], "-", label=r'$\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_W_2_real[:,-1], "-", label=r'$\langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_W_3_real[:,-1], "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_pure_phase_real[:,-1], "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (real)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t$')
    plt.grid()
    plt.legend()
    filename = 'components_no_phase_real_T=' + str(T) + '_tau=' + str(tau) + '.pdf'
    plt.savefig(filename)
    # plt.show()

    # Clear the plot
    plt.clf()

    # plt.figure(figsize=(25,4))
    plt.plot(ts, results_W_0_imag[:,-1], "-", label=r'$\langle W(t-\tau)\rangle$')
    plt.plot(ts, results_W_1_imag[:,-1], "-", label=r'$\langle W(t-\tau) W^\dagger(\tau)\rangle$')
    plt.plot(ts, results_W_2_imag[:,-1], "-", label=r'$\langle W(\tau)W(t-\tau)\rangle$')
    plt.plot(ts, results_W_3_imag[:,-1], "-", label=r'$\langle W(\tau) W(t-\tau) W^\dagger(\tau) \rangle$')

    plt.plot(ts, results_pure_phase_imag[:,-1], "-", label=r'$e^{iE\tau/\hbar}$')
    plt.title(r'$\tau\ =\ ' + str(tau) + '\ (imag)$')
    # plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$t$')
    plt.grid()
    plt.legend()
    filename = 'components_no_phase_imag_T=' + str(T) + '_tau=' + str(tau) + '.pdf'
    plt.savefig(filename)
    # plt.show()

if __name__ == "__main__":
    main()
