import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

filename = sys.argv[1]

file = open(filename, "r")
lines = file.readlines()

xs = []
ys = []

x = 0
for line in lines:
        # print("line: ", line)
        line_splitted = line.replace("\n", "").split(" | ")
        # print(line_splitted)
        k = line_splitted[0]
        gk_squared = line_splitted[1]
        # print("k: ", k)
        # print("gk_squared: ", gk_squared)
        # sys.exit()
        xs.append(float(k))
        ys.append(float(gk_squared))

# print("xs: ", xs)
# print("ys: ", ys)
# print("sum of ys: ", np.sum(ys))

gk_squared_max = max(ys)
index = ys.index(gk_squared_max)
k_max = xs[index]

print("gk_squared_max: ", gk_squared_max)
print("k_max: ", k_max)

period = k_max / 2

ks_final = np.arange(period, 2.5, period)
print("ks_final: ", ks_final)

selected_ks_filtered = [] # This list will store real values obtained during launching fortran caluclations
selected_gk_squared_filtered = []
for k in ks_final:
    for i in range(len(xs)):
        _k = xs[i]
        if abs(k - _k) < 1e-6:
            selected_ks_filtered.append(_k)
            selected_gk_squared_filtered.append(ys[i])
            break

# print("selected_ks_filtered: ", selected_ks_filtered)
# print("selected_gk_squared_filtered: ", selected_gk_squared_filtered)

# fig = plt.figure(figsize=[12, 12])
font_size = 20 # Changes the size of all fonts in the plot
tick_size = 20 # Changes the size of all labels on axes in the plot
plt.rc('font', size=font_size)
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

ax = plt.subplot()

# Add rectangles representing the sum approximating integral
for i in range(len(selected_ks_filtered)):
    k = selected_ks_filtered[i]
    gk_squared = selected_gk_squared_filtered[i]
    xy = (k - period/2, 0)# Anchor point
    width = period
    height = gk_squared
    rectangle = Rectangle(xy, width, height, facecolor='0.9', edgecolor='0.5')
    ax.add_patch(rectangle)

# ax.text(selected_ks_filtered[9] - (period/2) - 0.015, selected_gk_squared_filtered[9] + 0.005, r'$\}\alpha$', fontsize=25, rotation=90)

plt.plot(xs,ys, "-")
plt.plot(selected_ks_filtered, selected_gk_squared_filtered, '.')

plt.ylabel(r'$ G(k) \ [nm\ 10^{-2}]$', fontsize=font_size)
plt.xlabel(r'$ k\ [ nm^{-1} ]$', fontsize=font_size)
plt.ylim(0, 0.12)
plt.xlim(xs[0], xs[-1])
plt.tick_params(axis='both', labelsize=tick_size)

plt.xticks([0, 0.15, 0.28, 0.4, 1, 2, 2.4], [0, r'$k_0$', r'$k_1$', r'    $k_2 \cdots$', 1, 2, r'$\cdots k_{18}$'])
xticks = plt.xticks()[-1]
xticks[1].set_fontsize(15)
ax.xaxis.get_major_ticks()[1].tick1line.set_markersize(0)
xticks[2].set_fontsize(15)
ax.xaxis.get_major_ticks()[2].tick1line.set_markersize(0)
xticks[3].set_fontsize(15)
ax.xaxis.get_major_ticks()[3].tick1line.set_markersize(0)
xticks[-1].set_fontsize(15)
ax.xaxis.get_major_ticks()[-1].tick1line.set_markersize(0)

# print(type(xticks[-1]))
plt.yticks([0.05, 0.1], [5, 10])

plt.vlines(selected_ks_filtered, list(np.zeros(len(selected_ks_filtered))), selected_gk_squared_filtered, colors='black', linestyles='--')
# plt.grid()
plt.tight_layout()
# plt.subplots_adjust(left=0.15,
#                 bottom=0.18,
#                 right=0.96,
#                 top=0.9,
#                 wspace=0.02)
plt.savefig("gk_squared.pdf")
# plt.show()
