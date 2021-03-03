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

final_delta = np.mean(deltas_avg)
analytical_de_dig = 0.018206369779490493
print("DELTA: ",final_delta)
relative_err = abs(final_delta-analytical_de_dig)/analytical_de_dig
print(relative_err*100)
