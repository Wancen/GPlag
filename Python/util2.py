import pandas as pd
import time
import matplotlib.pyplot as plt
import numpy as np

def plotts(houseprice_process):
    # Assuming t, x1, and x2 are lists or arrays of the same length
    t = houseprice_process['period_end']
    x1 = houseprice_process['inventory_inverse_scaled']
    x2 = houseprice_process['price_scaled']

    plt.figure(figsize=(10, 6))

    # Plot x1 against t
    plt.plot(t, x1, label='inventory_inverse')

    # Plot x2 against t
    plt.plot(t, x2, label='median_price')

    # Adding a legend
    plt.legend()

    # Adding title and labels
    plt.title("Plot of inventory and price against t")
    plt.xlabel("t")
    plt.ylabel("Values")

    plt.show()

def ccf(x1, x2, maxlag):
    res = []
    n = len(x1)
    c_0 = np.sqrt(np.sum((x1 - np.mean(x1)) ** 2 / n) * np.sum((x2 - np.mean(x2)) ** 2 / n))
    for t in range(-maxlag, maxlag):
        if t <= 0:
            c_t = np.mean(np.roll(x1, n+t) * np.roll(x2, -t-1))
        else:
            c_t = np.mean(np.roll(x1, t) * np.roll(x2, -n+t))
        r_t = c_t / c_0
        res.append(r_t)
    return(res)