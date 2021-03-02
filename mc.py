import concurrent.futures
import numpy as np
import matplotlib.pyplot as plt

T = 1
K = 99
S = 100
r = 0.06
sigma = 0.2

def mc_trial(n_trials):
    payoff_list = []
    for n in range(n_trials):
        Z= np.random.normal(0,1)
        payoff = K - S*np.exp((r-0.5*(sigma)**2)*T+sigma*np.sqrt(T)*Z)
        if payoff < 0:
            payoff = 0
        payoff_list.append(payoff)
    return payoff_list

if __name__ == "__main__":
    n_trials = [10, 10**2, 10**3, 10**4, 10**5, 10**6]
    V_mean = []
    V_sterr = []
    for n in n_trials:
        value_list = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
            values = [executor.submit(mc_trial, n) for _ in range(20)]
            for f in concurrent.futures.as_completed(values):
                value_list.append(np.exp(-r*T)*np.mean(f.result()))

        V_mean.append(np.mean(value_list))
        V_sterr.append((np.std(value_list))/np.sqrt(len(value_list)))



    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(n_trials, V_mean)
    ax.fill_between(n_trials, np.array(V_mean) - np.array(V_sterr), np.array(V_mean) + np.array(V_sterr), alpha=0.5)
    ax.set_xlabel("N sample paths")
    ax.set_ylabel("Put option price in €")
    ax.set_xscale('log')
    ax.set_title("Convergence of MC simulation in determining the put option price")
    fig.savefig("mc_convergence_test.jpg")
    plt.show()

    file_aux  = open('mc.csv','a')
    file_aux.write("sigma k trials mean_price sterr_price")
    for i in range(len(V_mean)):
        file_aux.write("\n"+str()+" "+str()+" "+str(n_trials[i])+" "+str(V_mean[i])+" "+str(V_sterr[i]))
    file_aux.close()
