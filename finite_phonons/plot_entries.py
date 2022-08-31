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
    if "=== W_one_part ===" in line:
        break
    else:
        xs.append(x)
        x += 1
        ys.append(float(line.replace("\n", "")))

# print("xs: ", xs)
# print("ys: ", ys)
print("sum of ys: ", np.sum(ys))

plt.plot(xs,ys, ".")
plt.show()
