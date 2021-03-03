import concurrent.futures
import numpy as np
import matplotlib.pyplot as plt

T = 1
K = 99
S_0 = 100
r = 0.06
sigma = 0.2

def delta_montecarlo(a):
    # draw random number from the normal distribution
    Z = np.random.normal()

    # simulate the stock price at expiry
    S_T = S_0 * np.exp(((r-0.5*sigma**2)*T + sigma*np.sqrt(T)*Z))

    # we apply smoothing for a digital option and thus need
    # the derivative of the chosen sigmoid function wrt S_T
    dCdST = (a*np.exp(-a*(S_T-K)))/((np.exp(-a*(S_T-K))+1)**2)

    # calculate the delta for pathwise method
    delta = np.exp(-r*T)*dCdST*(S_T/S_0)

    return delta

# do monte carlo sim for the delta parameter
deltas_avg = []
for x in range(10):
    delta_list = []
    for i in range(100000):
        delta_list.append(delta_montecarlo(0.9))

    delta_param = np.mean(delta_list)
    #print(np.mean(delta_list), np.std(delta_list))

    deltas_avg.append(delta_param)

final_delta = np.mean(deltas_avg)
analytical_de_dig = 0.018206369779490493
relative_err = abs(final_delta-analytical_de_dig)/analytical_de_dig
print(relative_err*100)
