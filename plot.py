import matplotlib.pyplot as plt
import csv
import numpy as np

s5 = {"price": [], "err": [], "BSM": [], "tree": [], "err_bsm": []}
s10 = {"price": [], "err": [], "BSM": [], "tree": [], "err_bsm": []}
s20 = {"price": [], "err": [], "BSM": [], "tree": [], "err_bsm": []}
s30 = {"price": [], "err": [], "BSM": [], "tree": [], "err_bsm": []}
s40 = {"price": [], "err": [], "BSM": [], "tree": [], "err_bsm": []}

with open('mc_vol.csv', newline='') as f:
    r = csv.reader(f, delimiter=' ')
    next(r)
    for row in r:
        if row[0] == "0.05":
            s5["price"].append(float(row[2]))
            s5["err"].append(float(row[3]))
            s5["BSM"].append(float(row[4]))
            s5["tree"].append(float(row[5]))
            s5["err_bsm"].append(abs((float(row[4])-float(row[2]))/float(row[4]))*100)
        if row[0] == "0.1":
            s10["price"].append(float(row[2]))
            s10["err"].append(float(row[3]))
            s10["BSM"].append(float(row[4]))
            s10["tree"].append(float(row[5]))
            s10["err_bsm"].append(abs((float(row[4])-float(row[2]))/float(row[4]))*100)
        if row[0] == "0.2":
            s20["price"].append(float(row[2]))
            s20["err"].append(float(row[3]))
            s20["BSM"].append(float(row[4]))
            s20["tree"].append(float(row[5]))
            s20["err_bsm"].append(abs((float(row[4])-float(row[2]))/float(row[4]))*100)
        if row[0] == "0.3":
            s30["price"].append(float(row[2]))
            s30["err"].append(float(row[3]))
            s30["BSM"].append(float(row[4]))
            s30["tree"].append(float(row[5]))
            s30["err_bsm"].append(abs((float(row[4])-float(row[2]))/float(row[4]))*100)
        if row[0] == "0.4":
            s40["price"].append(float(row[2]))
            s40["err"].append(float(row[3]))
            s40["BSM"].append(float(row[4]))
            s40["tree"].append(float(row[5]))
            s40["err_bsm"].append(abs((float(row[4])-float(row[2]))/float(row[4]))*100)
x = [10, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7]
fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(121)
ax1.plot(x, s5["price"], label="σ = 5%")
ax1.fill_between(x, np.array(s5["price"]) - np.array(s5["err"]), np.array(s5["price"]) + np.array(s5["err"]), alpha=0.5)
ax1.plot(x, s10["price"], label="σ = 10%")
ax1.fill_between(x, np.array(s10["price"]) - np.array(s10["err"]), np.array(s10["price"]) + np.array(s10["err"]), alpha=0.5)
ax1.plot(x, s20["price"], label="σ = 20%")
ax1.fill_between(x, np.array(s20["price"]) - np.array(s20["err"]), np.array(s20["price"]) + np.array(s20["err"]), alpha=0.5)
ax1.plot(x, s30["price"], label="σ = 30%")
ax1.fill_between(x, np.array(s30["price"]) - np.array(s30["err"]), np.array(s30["price"]) + np.array(s30["err"]), alpha=0.5)
ax1.plot(x, s40["price"], label="σ = 40%")
ax1.fill_between(x, np.array(s40["price"]) - np.array(s40["err"]), np.array(s40["price"]) + np.array(s40["err"]), alpha=0.5)
ax1.set_xlabel("Number of trials")
ax1.set_ylabel("Put option price in €")
ax1.set_title("Convergence in MC simulation in put option value with different volatilities", fontsize=10)
ax1.set_xscale('log')
# ax1.legend()
ax = fig.add_subplot(122)
ax.plot(x, s5["err_bsm"], label="σ = 5%")
ax.plot(x, s10["err_bsm"], label="σ = 10%")
ax.plot(x, s20["err_bsm"], label="σ = 20%")
ax.plot(x, s30["err_bsm"], label="σ = 30%")
ax.plot(x, s40["err_bsm"], label="σ = 40%")
ax.set_ylabel("Relative error in %")
ax.set_xlabel("Number of trials")
ax.set_title("The relative error of a MC simulation compared to Black-Scholes", fontsize=10)
ax.set_xscale('log')
ax.legend()
fig.tight_layout(pad=1.0)
fig.savefig(f"volatility_test.jpg")

k5 = {"price": [], "err": [], "BSM": [], "tree": [], "err_bsm": []}
k10 = {"price": [], "err": [], "BSM": [], "tree": [], "err_bsm": []}
k20 = {"price": [], "err": [], "BSM": [], "tree": [], "err_bsm": []}
k30 = {"price": [], "err": [], "BSM": [], "tree": [], "err_bsm": []}

with open('mc_strike.csv', newline='') as f:
    r = csv.reader(f, delimiter=' ')
    next(r)
    for row in r:
        if row[0] == "79":
            k5["price"].append(float(row[2]))
            k5["err"].append(float(row[3]))
            k5["BSM"].append(float(row[4]))
            k5["tree"].append(float(row[5]))
            k5["err_bsm"].append(abs((float(row[4])-float(row[2]))/float(row[4]))*100)
        if row[0] == "99":
            k10["price"].append(float(row[2]))
            k10["err"].append(float(row[3]))
            k10["BSM"].append(float(row[4]))
            k10["tree"].append(float(row[5]))
            k10["err_bsm"].append(abs((float(row[4])-float(row[2]))/float(row[4]))*100)
        if row[0] == "119":
            k20["price"].append(float(row[2]))
            k20["err"].append(float(row[3]))
            k20["BSM"].append(float(row[4]))
            k20["tree"].append(float(row[5]))
            k20["err_bsm"].append(abs((float(row[4])-float(row[2]))/float(row[4]))*100)
        if row[0] == "139":
            k30["price"].append(float(row[2]))
            k30["err"].append(float(row[3]))
            k30["BSM"].append(float(row[4]))
            k30["tree"].append(float(row[5]))
            k30["err_bsm"].append(abs((float(row[4])-float(row[2]))/float(row[4]))*100)
x = [10, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7]
fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(121)
ax1.plot(x, k5["price"], label="K = €79")
ax1.fill_between(x, np.array(k5["price"]) - np.array(k5["err"]), np.array(k5["price"]) + np.array(k5["err"]), alpha=0.5)
ax1.plot(x, k10["price"], label="K = €99")
ax1.fill_between(x, np.array(k10["price"]) - np.array(k10["err"]), np.array(k10["price"]) + np.array(k10["err"]), alpha=0.5)
ax1.plot(x, k20["price"], label="K = €119")
ax1.fill_between(x, np.array(k20["price"]) - np.array(k20["err"]), np.array(k20["price"]) + np.array(k20["err"]), alpha=0.5)
ax1.plot(x, k30["price"], label="K = €139")
ax1.fill_between(x, np.array(k30["price"]) - np.array(k30["err"]), np.array(k30["price"]) + np.array(k30["err"]), alpha=0.5)
ax1.set_xlabel("Number of trials")
ax1.set_ylabel("Put option price in €")
ax1.set_title("Convergence in MC simulation in put option value with different strike prices", fontsize=10)
ax1.set_xscale('log')
ax = fig.add_subplot(122)
ax.plot(x, k5["err_bsm"], label="K = €79")
ax.plot(x, k10["err_bsm"], label="K = €99")
ax.plot(x, k20["err_bsm"], label="K = €119")
ax.plot(x, k30["err_bsm"], label="K = €139")
ax.set_ylabel("Relative error in %")
ax.set_xlabel("Number of trials")
ax.set_title("The relative error of a MC simulation compared to Black-Scholes", fontsize=10)
ax.set_xscale('log')
ax.legend()
fig.tight_layout(pad=1.0)
fig.savefig(f"strike_test.jpg")
