import matplotlib.pyplot as plt
import csv
import numpy as np

# analytical deltas as calculated by black-scholes
analytical_delta_put = -0.3262644882651039
analytical_delta_dig = 0.018206369779490493

# choose put or digital option
option_type = 'digital'
#option_type = 'put'

# read and sort the data
put = {'0.01':[], '0.05':[], '0.1':[], '0.5':[]}
digital = {'0.01':[], '0.05':[], '0.1':[], '0.5':[]}

# depending on whether to use the same seed or different seed data, choose between sameseed/difseed file
with open('bump_reval_difseed.csv') as f:
    csv_reader = csv.reader(f, delimiter=' ')
    for row in csv_reader:
        if row[0] == 'put':
            if row[1] == '0.01':
                put['0.01'].append((float(row[2]), float(row[3]), float(row[4])))

            elif row[1] == '0.05':
                put['0.05'].append((float(row[2]), float(row[3]), float(row[4])))

            elif row[1] == '0.1':
                put['0.1'].append((float(row[2]), float(row[3]), float(row[4])))

            elif row[1] == '0.5':
                put['0.5'].append((float(row[2]), float(row[3]), float(row[4])))

        elif row[0] == 'digital':
            if row[1] == '0.01':
                digital['0.01'].append((float(row[2]), float(row[3]), float(row[4])))

            elif row[1] == '0.05':
                digital['0.05'].append((float(row[2]), float(row[3]), float(row[4])))

            elif row[1] == '0.1':
                digital['0.1'].append((float(row[2]), float(row[3]), float(row[4])))

            elif row[1] == '0.5':
                digital['0.5'].append((float(row[2]), float(row[3]), float(row[4])))


# plot the data
if option_type == 'put':
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    for i in [0.01, 0.05, 0.1, 0.5]:
        trials = [trials[0] for trials in put[str(i)]]
        means = [trials[1] for trials in put[str(i)]]
        std_err = [trials[2] for trials in put[str(i)]]
        ax.plot(trials, means, label=f"epsilon {i}")
        ax.fill_between(trials, np.array(means) - np.array(std_err), np.array(means) + np.array(std_err), alpha=0.3)
    ax.set_xlabel("N sample paths")
    ax.set_ylabel("Put option delta")
    ax.set_xscale('log')
    ax.set_title("Convergence of MC simulation in determining put option delta")
    plt.axhline(y=-0.3262644882651039, color='black', ls='--', lw=0.7)
    plt.legend()
    #fig.savefig("mc_convergence_test.jpg")
    plt.show()


    # plot the data again, but in relative error from the analytical value
    # don't forget to change analytical_delta_put to to analytical_value_dig for digital option
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    for i in [0.01, 0.05, 0.1, 0.5]:
        trials = [trials[0] for trials in put[str(i)]]
        means = [100*abs(trials[1]-analytical_delta_put)/abs(analytical_delta_put) for trials in put[str(i)]]
        std_err = [trials[2] for trials in put[str(i)]]
        ax.plot(trials, means, label=f"epsilon {i}")
    ax.set_xlabel("N sample paths")
    ax.set_ylabel("Put option delta relative error (%)")
    ax.set_xscale('log')
    ax.set_title("Relative error between MC simulation and analytical value")
    plt.axhline(y=0, color='black', ls='--', lw=0.7)
    plt.legend()
    #fig.savefig("mc_convergence_test.jpg")
    plt.show()

else:
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    for i in [0.01, 0.05, 0.1, 0.5]:
        trials = [trials[0] for trials in digital[str(i)]]
        means = [trials[1] for trials in digital[str(i)]]
        std_err = [trials[2] for trials in digital[str(i)]]
        ax.plot(trials, means, label=f"epsilon {i}")
        ax.fill_between(trials, np.array(means) - np.array(std_err), np.array(means) + np.array(std_err), alpha=0.3)
    ax.set_xlabel("N sample paths")
    ax.set_ylabel("Digital option delta")
    ax.set_xscale('log')
    ax.set_title("Convergence of MC simulation in determining digital option delta")
    plt.axhline(y=0.018206369779490493, color='black', ls='--', lw=0.7)
    plt.legend()
    #fig.savefig("mc_convergence_test.jpg")
    plt.show()


    # plot the data again, but in relative error from the analytical value
    # don't forget to change analytical_delta_put to to analytical_value_dig for digital option
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    for i in [0.01, 0.05, 0.1, 0.5]:
        trials = [trials[0] for trials in digital[str(i)]]
        means = [100*abs(trials[1]-analytical_delta_dig)/abs(analytical_delta_dig) for trials in digital[str(i)]]
        std_err = [trials[2] for trials in digital[str(i)]]
        ax.plot(trials, means, label=f"epsilon {i}")
    ax.set_xlabel("N sample paths")
    ax.set_ylabel("Digital option delta relative error (%)")
    ax.set_xscale('log')
    ax.set_title("Relative error between MC simulation and analytical value")
    plt.axhline(y=0, color='black', ls='--', lw=0.7)
    plt.legend()
    #fig.savefig("mc_convergence_test.jpg")
    plt.show()

