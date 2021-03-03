import numpy as np
import random
import time
from scipy.stats import norm

def payoff_bumped(S_0, K, T, r, sigma, diff_seeds, epsilon, option_type):
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
        if S_T1 + epsilon > K:
            payoff_bumped = np.exp(-r)
        else:
            payoff_bumped = 0

    else:  
        payoff_bumped = np.exp(-r)*max(0, K-(S_T1+epsilon))

    return payoff_bumped

def payoff_unbumped(S_0, K, T, r, sigma, diff_seeds, epsilon, option_type):
    """
    Returns a call option bumped and unbumped payoff based on spot price, 
    strike price, time till expiry, risk-free interest rate, and volatility.
    """

    # draw two random numbers from the normal distribution
    # for bumped and unbumped prices
    Z_1 = np.random.normal()

    # simulate the stock price at expiry
    S_T1 = S_0 * np.exp(((r-0.5*sigma**2)*T + sigma*np.sqrt(T)*Z_1))

    #print(payoff_bumped, payoff_unbumped)
    # analytical_d = norm.cdf(-(np.log(S_0/K) + (r+((0.2)**2)/2)*(T)/(0.2*np.sqrt(T))))
    # print("Analytical delta:", analytical_d)
    # d2 = (np.log(S_0/K)+(r-0.5*sigma**2)*T)/(sigma*np.sqrt(T))
    # digital_option_analytical_price = np.exp(-r)*norm.cdf(d2)
    # digital_option_analytical_delta = (np.exp(-r)*norm.pdf(d2))/(S_0*sigma*np.sqrt(T))

    # print(digital_option_analytical_price, digital_option_analytical_delta)

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

delta_means = []
bumped_list = []
unbumped_list = []
epsilon = 0.5
same_seed = True
option_type = 'digital'

# if different seeds, can run in the same loop
if same_seed == False:
    for i in range(1000):
        
        #bumped, unbumped = payoff(100, 99, 1, 0.06, 0.2, True, epsilon)
        bumped = payoff_bumped(100, 99, 1, 0.06, 0.2, True, epsilon, option_type)
        unbumped = payoff_unbumped(100, 99, 1, 0.06, 0.2, True, epsilon, option_type)
        
        bumped_list.append(bumped)
        unbumped_list.append(unbumped)

# if same seeds, run both loops with same seeds
else:
    a = int(np.random.uniform()*1000000)
    np.random.seed(a)
    for i in range(100000):
        bumped = payoff_bumped(100, 99, 1, 0.06, 0.2, True, epsilon, option_type)
        bumped_list.append(bumped)

    # if different seeds, can run in the same loop
    np.random.seed(a)
    for i in range(100000):
        unbumped = payoff_unbumped(100, 99, 1, 0.06, 0.2, True, epsilon, option_type)
        unbumped_list.append(unbumped)


#analytical_d = 0.659115488305189
analytical_d = -0.34088451169481104
analytical_pr_dig = 0.5639320297620053 
analytical_de_dig = 0.018206369779490493



delta_param = delta_bump_reval(np.mean(bumped_list), np.mean(unbumped_list), epsilon)
print(np.mean(unbumped_list))
print(delta_param)
relative_err = abs(delta_param-analytical_d)/analytical_d
relative_err = abs(delta_param-analytical_de_dig)/analytical_de_dig

print(relative_err*100)
#print(mean)

