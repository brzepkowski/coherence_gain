import sys
import numpy as np
import matplotlib.pyplot as plt

filename = sys.argv[1]

file = open(filename, "r")
lines = file.readlines()

xs = []
ys = []

x = 0
for line in lines:
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



plt.plot(xs,ys, ".")
plt.grid()
plt.show()
