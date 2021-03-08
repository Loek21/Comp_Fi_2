import matplotlib.pyplot as plt
import csv
import numpy as np

# analytical delta as calculated by black-scholes
analytical_delta_dig = 0.018206369779490493

# read and sort the data
smoothing = {'1':[], '2':[], '3':[], '4':[], '5':[]}

with open('pathwise.csv') as f:
    csv_reader = csv.reader(f, delimiter=' ')
    for row in csv_reader:
        if row[0] == '1':
            smoothing['1'].append((float(row[1]), float(row[2]), float(row[3])))

        elif row[0] == '2':
            smoothing['2'].append((float(row[1]), float(row[2]), float(row[3])))

        elif row[0] == '3':
            smoothing['3'].append((float(row[1]), float(row[2]), float(row[3])))

        elif row[0] == '4':
            smoothing['4'].append((float(row[1]), float(row[2]), float(row[3])))

        elif row[0] == '5':
            smoothing['5'].append((float(row[1]), float(row[2]), float(row[3])))

#plot the data
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
for i in range(1, 6):
    trials = [trials[0] for trials in smoothing[str(i)]]
    means = [trials[1] for trials in smoothing[str(i)]]
    std_err = [trials[2] for trials in smoothing[str(i)]]
    ax.plot(trials, means, label=f"Smoothing {i}")
    ax.fill_between(trials, np.array(means) - np.array(std_err), np.array(means) + np.array(std_err), alpha=0.3)
ax.set_xlabel("N sample paths")
ax.set_ylabel("Digital option delta")
ax.set_xscale('log')
ax.set_title("Convergence of MC simulation in determining digital option delta")
plt.axhline(y=0.018206369779490493, color='black', ls='--', lw=0.7)
plt.legend()
plt.show()


# plot the data again, but in relative error from the analytical value
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
for i in range(1, 6):
    trials = [trials[0] for trials in smoothing[str(i)]]
    means = [100*abs(trials[1]-analytical_delta_dig)/analytical_delta_dig for trials in smoothing[str(i)]]
    std_err = [trials[2] for trials in smoothing[str(i)]]
    ax.plot(trials, means, label=f"Smoothing {i}")
ax.set_xlabel("N sample paths")
ax.set_ylabel("Digital option delta relative error (%)")
ax.set_xscale('log')
ax.set_title("Relative error between MC simulation and analytical value")
plt.axhline(y=0, color='black', ls='--', lw=0.7)
plt.legend()
plt.show()

