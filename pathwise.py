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

def run_sim_pathwise(n_trials, smoothing):
# do monte carlo sim for the delta parameter
    
    delta_list = []
    for _ in range(n_trials):
        delta_list.append(delta_montecarlo(smoothing))

    delta_param = np.mean(delta_list)

    return delta_param

if __name__=="__main__":

    file_aux  = open(f'pathwise.csv','a')
    #file_aux.write("smoothing n_trials delta_mean delta_sterr")
    trials_list = [10, 100, 1000, 10000, 100000, 1000000]
    for smoothing in [4,5]:
        for trials in trials_list:
            delta_list = []
            with concurrent.futures.ProcessPoolExecutor(max_workers=12) as executor:
                values = [executor.submit(run_sim_pathwise, trials, smoothing) for _ in range(100)]
                for f in concurrent.futures.as_completed(values):
                    delta_list.append(np.mean(f.result()))

            delta_mean = np.mean(delta_list)
            delta_sterr = np.std(delta_list)/10
            file_aux  = open(f'pathwise.csv','a')
            file_aux.write("\n"+str(smoothing)+" "+str(trials)+" "+str(delta_mean)+" "+str(delta_sterr))
    file_aux.close()


    analytical_de_dig = 0.018206369779490493
    #relative_err = abs(final_delta-analytical_de_dig)/analytical_de_dig
    #print(relative_err*100)
