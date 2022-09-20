import numpy as np
import sys

def envelopes(xs, ys, epsilon_min=1e-5, epsilon_max=1e-8, diff=0.015):
    previous_data_point = None
    tendency = None
    min_xs = []
    min_ys = []
    max_xs = []
    max_ys = []

    for i in range(len(ys)):
        data_point = ys[i]
        if previous_data_point is not None:
            if data_point - previous_data_point > 0:
                if tendency is None:
                    tendency = "INCREASE"
                elif tendency == "DECREASE" and abs(data_point - previous_data_point) > epsilon_min:
                    if min_ys == [] or abs(data_point - min_ys[-1]) < diff:
                        tendency = "INCREASE"
                        # print("### CHANGE IN TENDENCY ###")
                        min_xs.append(xs[i-1])
                        min_ys.append(ys[i-1])
            elif data_point - previous_data_point < 0:
                if tendency is None:
                    tendency = "DECREASE"
                elif tendency == "INCREASE" and abs(data_point - previous_data_point) > epsilon_max:
                    if max_ys == [] or abs(data_point - max_ys[-1]) < diff:
                        tendency = "DECREASE"
                        # print("### CHANGE IN TENDENCY ###")
                        max_xs.append(xs[i-1])
                        max_ys.append(ys[i-1])
        previous_data_point = data_point
    return min_xs, min_ys, max_xs, max_ys

def combine_two_envelopes(xs_0, ys_0, xs_1, ys_1):
    xs = []
    ys = []
    while xs_0 != [] and xs_1 != []:
        if xs_0[0] < xs_1[0]:
            xs.append(xs_0.pop(0))
            ys.append(ys_0.pop(0))
        else:
            xs.append(xs_1.pop(0))
            ys.append(ys_1.pop(0))

    if xs_0 != []:
        for i in range(len(xs_0)):
            xs.append(xs_0[i])
            ys.append(ys_0[i])
    else:
        for i in range(len(xs_1)):
            xs.append(xs_1[i])
            ys.append(ys_1[i])

    return xs, ys
