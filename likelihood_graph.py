import matplotlib.pyplot as plt
import csv
import numpy as np

# analytical delta as calculated by black-scholes
analytical_delta_dig = 0.018206369779490493

# read and sort the data
data = []

with open('likelihood.csv') as f:
    csv_reader = csv.reader(f, delimiter=' ')
    row_nr = 0
    for row in csv_reader:
        if row_nr > 0:
            data.append((float(row[0]), float(row[1]), float(row[2])))
        row_nr += 1

#plot the data
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)

trials = [trials[0] for trials in data]
means = [trials[1] for trials in data]
std_err = [trials[2] for trials in data]
ax.plot(trials, means)
ax.fill_between(trials, np.array(means) - np.array(std_err), np.array(means) + np.array(std_err), alpha=0.3)
ax.set_xlabel("N sample paths")
ax.set_ylabel("Digital option delta")
ax.set_xscale('log')
ax.set_title("Convergence of MC simulation in determining digital option delta")
plt.axhline(y=0.018206369779490493, color='black', ls='--', lw=0.7)
plt.show()


# plot the data again, but in relative error from the analytical value
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
trials = [trials[0] for trials in data]
means = [100*abs(trials[1]-analytical_delta_dig)/analytical_delta_dig for trials in data]
std_err = [trials[2] for trials in data]
ax.plot(trials, means)
ax.set_xlabel("N sample paths")
ax.set_ylabel("Digital option delta relative error (%)")
ax.set_xscale('log')
ax.set_title("Relative error between MC simulation and analytical value")
plt.axhline(y=0, color='black', ls='--', lw=0.7)
plt.show()

