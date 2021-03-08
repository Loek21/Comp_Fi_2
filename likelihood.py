import concurrent.futures
import numpy as np
import matplotlib.pyplot as plt

T = 1
K = 99
S_0 = 100
r = 0.06
sigma = 0.2

def delta_montecarlo_likelihood():
    # draw random number from the normal distribution for the stock price
    # and likelihood ration method
    Z = np.random.normal()

    # simulate the stock price at expiry
    S_T = S_0 * np.exp(((r-0.5*sigma**2)*T + sigma*np.sqrt(T)*Z))

    # get the indicator function for the option payoff
    indicator = 1 if S_T > K else 0

    # calculate the delta for likelihood ratio method
    delta = np.exp(-r*T)*indicator*(Z/(sigma*S_0*np.sqrt(T)))
    #print(delta)

    return delta

# do monte carlo sim for the delta parameter
deltas_avg = []
for x in range(10):
    delta_list = []
    for i in range(10000):
        delta = delta_montecarlo_likelihood()
        delta_list.append(delta)

    delta_param = np.mean(delta_list)
    #print(np.mean(delta_list), np.std(delta_list))

    deltas_avg.append(delta_param)


def run_sim_likelihood(n_trials):
# do monte carlo sim for the delta parameter
    
    delta_list = []
    for _ in range(n_trials):
        delta_list.append(delta_montecarlo_likelihood())

    delta_param = np.mean(delta_list)

    return delta_param

if __name__=="__main__":

    file_aux  = open(f'likelihood.csv','a')
    file_aux.write("n_trials delta_mean delta_sterr")
    trials_list = [10, 100, 1000, 10000, 100000, 1000000]
    for trials in trials_list:
        delta_list = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
            values = [executor.submit(run_sim_likelihood, trials) for _ in range(100)]
            for f in concurrent.futures.as_completed(values):
                delta_list.append(np.mean(f.result()))

        delta_mean = np.mean(delta_list)
        delta_sterr = np.std(delta_list)/10
        file_aux  = open(f'likelihood.csv','a')
        file_aux.write("\n"+str(trials)+" "+str(delta_mean)+" "+str(delta_sterr))
    file_aux.close()


    final_delta = np.mean(deltas_avg)
    analytical_de_dig = 0.018206369779490493
    print("DELTA: ",final_delta)
    relative_err = abs(final_delta-analytical_de_dig)/analytical_de_dig
    print(relative_err*100)
