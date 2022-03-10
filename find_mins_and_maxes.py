from matplotlib import pyplot as plt
import numpy as np
import sys, os, glob

PERIOD = 0.00413566781

def main():
    os.chdir("./")
    input_file_prefix = "g_av_vs_tau_*"
    input_file_postfix = ".dat"
    for filename in glob.glob(input_file_prefix):
        if input_file_postfix in filename:
            print("filename: ", filename)

            taus = []
            g_avs = []

            taus_min = []
            g_avs_min = []
            taus_max = []
            g_avs_max = []

            prev_prev_g_av = None
            prev_g_av = None
            prev_tau = None

            input_file = open(filename, "r")
            content = input_file.readlines()
            for line in content:
                splitted_line = line.split()
                tau = float(splitted_line[0])
                g_av = float(splitted_line[1])

                taus.append(tau)
                g_avs.append(g_av)

                if prev_g_av is None:
                    prev_g_av = g_av
                    prev_tau = tau
                else:
                    if prev_prev_g_av is not None:
                        if prev_g_av - prev_prev_g_av >= 0: # g_av was ascending one timestep before
                            if g_av - prev_g_av <= 0: # Change from ascend to descend
                                taus_max.append(prev_tau)
                                g_avs_max.append(prev_g_av)
                        if prev_g_av - prev_prev_g_av <= 0: # g_av was descending one timestep before
                            if g_av - prev_g_av >= 0: # Change from descend to ascend
                                taus_min.append(prev_tau)
                                g_avs_min.append(prev_g_av)

                    prev_prev_g_av = prev_g_av
                    prev_g_av = g_av
                    prev_tau = tau

            ####################################################################
            # MIN
            ####################################################################

            # Sieve through the results to get rid of values coming from numerical inaccuracies
            tau_min_batches = []
            g_av_min_batches = []
            tau_min_batch = []
            g_av_min_batch = []
            for i in range(len(taus_min)):
                if tau_min_batch == []:
                    tau_min_batch.append(taus_min[i])
                    g_av_min_batch.append(g_avs_min[i])
                else:
                    if taus_min[i] - taus_min[i-1] <= PERIOD/6:
                        tau_min_batch.append(taus_min[i])
                        g_av_min_batch.append(g_avs_min[i])
                    else:
                        tau_min_batches.append(tau_min_batch)
                        g_av_min_batches.append(g_av_min_batch)
                        tau_min_batch = [taus_min[i]]
                        g_av_min_batch = [g_avs_min[i]]
            tau_min_batches.append(tau_min_batch)
            g_av_min_batches.append(g_av_min_batch)

            for i in range(len(tau_min_batches)):
                tau_min_batch = tau_min_batches[i]
                g_av_min_batch = g_av_min_batches[i]

                min_g_av = np.inf
                min_g_av_index = None
                for j in range(len(g_av_min_batch)):
                    if g_av_min_batch[j] < min_g_av:
                        min_g_av = g_av_min_batch[j]
                        min_g_av_index = j

                tau_min_batches[i] = [tau_min_batch[min_g_av_index]]
                g_av_min_batches[i] = [g_av_min_batch[min_g_av_index]]

            ####################################################################
            # MAX
            ####################################################################

            # Sieve through the results to get rid of values coming from numerical inaccuracies
            tau_max_batches = []
            g_av_max_batches = []
            tau_max_batch = []
            g_av_max_batch = []
            for i in range(len(taus_max)):
                # print("tau: ", taus_max[i])
                if tau_max_batch == []:
                    tau_max_batch.append(taus_max[i])
                    g_av_max_batch.append(g_avs_max[i])
                else:
                    if taus_max[i] - taus_max[i-1] <= PERIOD/6:
                        tau_max_batch.append(taus_max[i])
                        g_av_max_batch.append(g_avs_max[i])
                    else:
                        tau_max_batches.append(tau_max_batch)
                        g_av_max_batches.append(g_av_max_batch)
                        tau_max_batch = [taus_max[i]]
                        g_av_max_batch = [g_avs_max[i]]
            tau_max_batches.append(tau_max_batch)
            g_av_max_batches.append(g_av_max_batch)

            for i in range(len(tau_max_batches)):
                tau_max_batch = tau_max_batches[i]
                g_av_max_batch = g_av_max_batches[i]

                max_g_av = -np.inf
                max_g_av_index = None
                for j in range(len(g_av_max_batch)):
                    if g_av_max_batch[j] > max_g_av:
                        max_g_av = g_av_max_batch[j]
                        max_g_av_index = j

                tau_max_batches[i] = [tau_max_batch[max_g_av_index]]
                g_av_max_batches[i] = [g_av_max_batch[max_g_av_index]]

            taus_min = [val for sublist in tau_min_batches for val in sublist]
            g_avs_min = [val for sublist in g_av_min_batches for val in sublist]

            taus_max = [val for sublist in tau_max_batches for val in sublist]
            g_avs_max = [val for sublist in g_av_max_batches for val in sublist]


            ####################################################################
            # Final sieve to get rid of minima lying close to maxima. This should
            # be only necessary for small values of tau, where the g_av function
            # is kind of flat at the top). For larger taus we get a function similar
            # to sine, where this shouldn't be needed.
            ####################################################################

            taus_min_temp = taus_min.copy()
            g_avs_min_temp = g_avs_min.copy()
            for i in range(len(taus_min_temp)):
                tau_min = taus_min_temp[i]
                g_av_min = g_avs_min_temp[i]
                for tau_max in taus_max:
                    if abs(tau_min - tau_max) < PERIOD/8:
                        taus_min.remove(tau_min)
                        g_avs_min.remove(g_av_min)

            # print("taus_min: ", taus_min)
            # print("taus_max: ", taus_max)

            ####################################################################
            # Find global minimum
            ####################################################################

            global_min_g_av = min(g_avs_min)
            global_min_g_av_index = g_avs_min.index(global_min_g_av)
            print("g_av global minimum: ", global_min_g_av, " (",  global_min_g_av_index, ")")
            print("Tau of g_av's global minimum: ", taus_min[global_min_g_av_index])
            for i in range(len(taus_max)):
                tau = taus_max[i]
                if tau > taus_min[global_min_g_av_index]:
                    break
            print("Tau of g_av's local maximum close to tau of global minimum: ", taus_max[i])

            ####################################################################
            # Save data to files
            ####################################################################

            plt.plot(taus, g_avs, "-", label=r"$g_{av}$")
            plt.plot(taus_min, g_avs_min, ".", label=r"$g_{av}^{MIN}$")
            plt.plot(taus_max, g_avs_max, ".", label=r"$g_{av}^{MAX}$")

            # Mark values used in the generation of g_av vs. t plots
            plt.plot(taus_min[-1], g_avs_min[-1], 'bx', label=r"$picked\ g_{av}^{MIN}$")
            plt.plot(taus_max[-1], g_avs_max[-1], 'rx', label=r"$picked\ g_{av}^{MAX}$")

            plt.legend()
            plt.grid()
            # Save to file
            filename_plot = filename.replace(".dat", "_mins_and_maxes.pdf")
            # plt.savefig(filename_plot)
            plt.clf()
            # plt.show()

            """
            # Save obtained mins and maxes to appropriate files
            filename_mins = filename.replace(".dat", "_mins.txt")
            output_file = open(filename_mins, "w")
            for tau in taus_min:
                output_file.write(str(tau) + "\n")
            output_file.close()

            filename_maxes = filename.replace(".dat", "_maxes.txt")
            output_file = open(filename_maxes, "w")
            for tau in taus_max:
                output_file.write(str(tau) + "\n")
            output_file.close()
            """
            # sys.exit()

if __name__ == '__main__':
    main()
