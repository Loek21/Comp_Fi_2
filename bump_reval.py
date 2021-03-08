import numpy as np
import random
import time
from scipy.stats import norm
import concurrent.futures

def payoff_bumped(S_0, K, T, r, sigma, epsilon, option_type):
    """
    Returns a call option bumped and unbumped payoff based on spot price, 
    strike price, time till expiry, risk-free interest rate, and volatility.
    """

    # draw two random numbers from the normal distribution
    # for bumped and unbumped prices
    Z_1 = np.random.normal()

    # simulate the stock price at expiry with a bumped initial value
    S_T1 = (S_0+epsilon) * np.exp(((r-0.5*sigma**2)*T + sigma*np.sqrt(T)*Z_1))

    if option_type == 'digital':
        if S_T1 > K:
            payoff_bumped = np.exp(-r)
        else:
            payoff_bumped = 0

    else:  
        payoff_bumped = np.exp(-r)*max(0, K-(S_T1))

    return payoff_bumped

def payoff_unbumped(S_0, K, T, r, sigma, epsilon, option_type):
    """
    Returns a call option bumped and unbumped payoff based on spot price, 
    strike price, time till expiry, risk-free interest rate, and volatility.
    """

    # draw two random numbers from the normal distribution
    # for bumped and unbumped prices
    Z_1 = np.random.normal()

    # simulate the stock price at expiry
    S_T1 = S_0 * np.exp(((r-0.5*sigma**2)*T + sigma*np.sqrt(T)*Z_1))

    if option_type == 'digital':
        if S_T1 > K:
            payoff_unbumped = np.exp(-r)
        else:
            payoff_unbumped = 0

    else:  
        payoff_unbumped = np.exp(-r)*max(0, K-S_T1)

    return payoff_unbumped

def delta_bump_reval(payoff_bumped, payoff_unbumped, epsilon):
    """
    Returns the delta hedge parameter according to the bump 
    and revaluate method.
    """

    delta = (payoff_bumped-payoff_unbumped)/epsilon

    return delta

def run_sim(option_type, epsilon, same_seed, trials):
    bumped_list = []
    unbumped_list = []
    epsilon = epsilon
    same_seed = same_seed
    option_type = option_type

    # if different seeds, can run in the same loop
    if same_seed == False:
        for _ in range(trials):
            
            bumped = payoff_bumped(100, 99, 1, 0.06, 0.2, epsilon, option_type)
            unbumped = payoff_unbumped(100, 99, 1, 0.06, 0.2, epsilon, option_type)
            
            bumped_list.append(bumped)
            unbumped_list.append(unbumped)

    # if same seeds, run both loops with same seeds
    else:
        a = int(np.random.uniform()*1000000)
        np.random.seed(a)
        for _ in range(trials):
            bumped = payoff_bumped(100, 99, 1, 0.06, 0.2, epsilon, option_type)
            bumped_list.append(bumped)

        # if different seeds, can run in the same loop
        np.random.seed(a)
        for _ in range(trials):
            unbumped = payoff_unbumped(100, 99, 1, 0.06, 0.2, epsilon, option_type)
            unbumped_list.append(unbumped)

    delta_param = delta_bump_reval(np.mean(bumped_list), np.mean(unbumped_list), epsilon)

    return delta_param

if __name__=="__main__":
    
    # open file and save data there
    file_aux  = open(f'bump_reval_difseed.csv','a')
    file_aux.write("option epsilon n_trials delta_mean delta_sterr")
    epsilon_list = [0.01, 0.05, 0.1, 0.5]
    option_list = ['put', 'digital']
    trials_list = [10, 100, 1000, 10000, 100000, 1000000]
    for option in option_list:
        for epsilon in epsilon_list:
            for trials in trials_list:
                print(option, epsilon, trials)
                delta_list = []

                # parallel computing
                with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
                    values = [executor.submit(run_sim, option, epsilon, False, trials) for _ in range(100)]
                    for f in concurrent.futures.as_completed(values):
                        delta_list.append(np.mean(f.result()))

                delta_mean = np.mean(delta_list)
                delta_sterr = np.std(delta_list)/10
                file_aux  = open(f'bump_reval_difseed.csv','a')
                file_aux.write("\n"+str(option)+" "+str(epsilon)+" "+str(trials)+" "+str(delta_mean)+" "+str(delta_sterr))
    file_aux.close()
        

