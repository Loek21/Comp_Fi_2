import concurrent.futures
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from tree import binomial_tree, put_price_tree

T = 1
# K = 99
S = 100
r = 0.06
sigma = 0.2

def binom_tree(S, K, T, r, sigma):
    """
    Determines the put option price using the Binomial tree method.
    :param S: stock price at time 0
    :param K: strike price option
    :param T: time of expiry
    :param r: risk-free interest rate
    :param sigma: volatility
    :return put option price
    """
    tree = binomial_tree(S, K, T, r, sigma, 100)
    price = put_price_tree(tree, K, T, r, sigma, 100)
    return price[0,0]

def black_scholes(S, K, T, sigma, r):
    """
    Determines the put option price using Black-Scholes model.
    :param S: stock price at time 0
    :param K: strike price option
    :param T: time of expiry
    :param r: risk-free interest rate
    :param sigma: volatility
    :return put option price
    """
    d_1 = (np.log(S/K) + (r+((sigma)**2)/2)*(T))/(sigma*np.sqrt(T))
    d_2 = d_1 - sigma*np.sqrt(T)
    price =  np.exp(-r*(T))*K*norm.cdf(-d_2) - S*norm.cdf(-d_1)
    return price

def mc_trial(n_trials, sigma, K):
    """
    Determines the put option payoff for each trial in MC simulation.
    :param n_trials: number of trials in MC simulation
    :param sigma: volatility
    :param K: strike price option
    :return list of payoffs
    """
    payoff_list = []
    for n in range(n_trials):
        Z= np.random.normal(0,1)
        payoff = K - S*np.exp((r-0.5*(sigma)**2)*T+sigma*np.sqrt(T)*Z)
        if payoff < 0:
            payoff = 0
        payoff_list.append(payoff)
    return payoff_list

if __name__ == "__main__":
    # volatility = [0.05, 0.1 ,0.2, 0.3, 0.4]
    strike = [79, 99, 119, 139]
    n_trials = [10, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7]
    file_aux  = open('mc_strike.csv','a')
    file_aux.write("sigma trials mean_price sterr_price black_scholes_price tree_price")
    # for vol in volatility:
    for K in strike:
        V_mean = []
        V_sterr = []
        for n in n_trials:
            value_list = []
            with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
                values = [executor.submit(mc_trial, n, sigma, K) for _ in range(25)]
                for f in concurrent.futures.as_completed(values):
                    value_list.append(np.exp(-r*T)*np.mean(f.result()))

            V_mean.append(np.mean(value_list))
            V_sterr.append((np.std(value_list))/np.sqrt(len(value_list)))

        BS_price = black_scholes(S, K, T, sigma, r)
        tree_price = binom_tree(S, K, T, r, sigma)
        for i in range(len(V_mean)):
            file_aux.write("\n"+str(K)+" "+str(n_trials[i])+" "+str(V_mean[i])+" "+str(V_sterr[i])+" "+str(BS_price)+" "+str(tree_price))
    file_aux.close()
