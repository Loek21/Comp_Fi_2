import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy.stats import norm
from scipy.stats.mstats import gmean

# Monte carlo Asian
S0 = 100
r = 0.06
T = 1
sigma = 0.2


def sim(n, m, K) -> tuple:
    """
    Creates a stock price vector with m dimensions and determines the geometric and arithmetic Asian call option.
    Also calculates the control variate
    :param n: Number of stock paths to take the average from
    :param m: The number of steps in the stock path, or the dimensions of the stock price vector
    :param K: The strike price
    :return: Tuple with fair values (arithmetic, geometric, control variate)
    """
    sum_ari = 0
    sum_geo = 0
    ari_values = []
    geo_values = []
    dT = T / m

    for i in range(n):
        stock_price = [S0]
        S_temp = S0

        for j in range(m):
            phi = np.random.normal(0, 1)
            S_next = S_temp + S_temp * (r * dT + sigma * phi * np.sqrt(dT))
            stock_price.append(S_next)
            S_temp = S_next

        C_ari = max(0, 1 / (m + 1) * sum(stock_price) - K)
        C_geo = max(0, gmean(stock_price) - K)

        ari_values.append(C_ari)
        geo_values.append(C_geo)

        sum_ari += C_ari
        sum_geo += C_geo

    C_arithmetic = np.exp(-r * T) * sum_ari / n
    C_geometric = np.exp(-r * T) * sum_geo / n

    ari_error = np.std(ari_values) / np.sqrt(n)
    geo_error = np.std(geo_values) / np.sqrt(n)

    # Optimal coefficient for Control variate
    covariance_matrix = np.cov(geo_values, ari_values)
    c = covariance_matrix[0, 1] / np.var(geo_values)

    control_variate = C_arithmetic - c * (
        C_geometric - asianGeometricOption(S0, T, m, K)
    )

    return C_arithmetic, C_geometric, control_variate, ari_error, geo_error


# Geometric Asian option
def asianGeometricOption(S0, T, m, K) -> float:
    """
    Determines the analytical value of the geometric Asian call option
    :param S0: The starting price of the underlying
    :param T: The maturity date or time
    :param m: The number of steps in the stock path, or the dimensions of the stock price vector
    :param K: The strike price
    :return: The fair value of the geometric Asian call option
    """

    sigma_c = sigma * np.sqrt((2 * m + 1) / (6 * m + 6))
    r_c = ((r - 0.5 * sigma ** 2) + sigma_c ** 2) / 2

    d1 = (np.log(S0 / K) + (r_c + 0.5 * sigma_c ** 2) * T) / (np.sqrt(T) * sigma_c)
    d2 = (np.log(S0 / K) + (r_c - 0.5 * sigma_c ** 2) * T) / (np.sqrt(T) * sigma_c)

    return np.exp(-r * T) * (S0 * np.exp(r_c * T) * norm.cdf(d1) - K * norm.cdf(d2))


# Compare geometric and arithmetic asian option as function of m
K = 99
m_list = np.arange(10, 101, 10)

n = 10000  # Samples or number of paths


temp_cont = []
temp_ari = []
temp_geo = []

for m in tqdm(m_list):
    arithmetic, geometric, control_variate, ari_error, geo_error = sim(n, m, K)
    temp_cont.append(control_variate)
    temp_ari.append(arithmetic)
    temp_geo.append(geometric)


plt.plot(m_list, temp_cont, label=f"Control K = {K}")
plt.plot(m_list, temp_ari, label=f"Arithmetic K = {K}")
plt.plot(m_list, temp_geo, label=f"Geometric K = {K}")


plt.xlabel("m", fontsize=22)
plt.ylabel("Option price", fontsize=25)
plt.title("Asian call option price", fontsize=25)
plt.legend(fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

K = 99
n_list = np.arange(1000, 106000, 5000)
m = 80  # Samples or number of paths

temp_cont = []
temp_ari = []
temp_geo = []

for n in tqdm(n_list):
    arithmetic, geometric, control_variate, ari_error, geo_error = sim(n, m, K)
    temp_cont.append(control_variate)
    temp_ari.append(arithmetic)
    temp_geo.append(geometric)

plt.plot(n_list, temp_cont, label=f"Control K = {K}")
plt.plot(n_list, temp_ari, label=f"Arithmetic K = {K}")
plt.plot(n_list, temp_geo, label=f"Geometric K = {K}")


plt.xlabel("n", fontsize=25)
plt.ylabel("Option price", fontsize=25)
plt.title("Asian call option price", fontsize=25)
plt.xscale("log")
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend()
plt.show()


# Compare geometric and arithmetic asian option as function of n
K_list = np.arange(79, 140, 20)
n_list = np.arange(1000, 101000, 5000)

difference2 = []
control_variate_list2 = []

m = 80

for K in K_list:
    temp_diff2 = []
    temp_cont2 = []

    for n in tqdm(n_list):

        arithmetic, geometric, control_variate, ari_error, geo_error = sim(n, m, K)
        temp_diff2.append(abs(geometric - arithmetic))
        temp_cont2.append(abs(control_variate - arithmetic))

    difference2.append(temp_diff2)
    control_variate_list2.append(temp_cont2)


for i, diff in enumerate(difference2):
    plt.plot(n_list, diff, label=f"{K_list[i]}")

plt.xlabel("n", fontsize=25)
plt.ylabel("error", fontsize=25)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xscale("log")
plt.legend(fontsize=22)
plt.title("Difference between geometric and arithmetic average", fontsize=25)
plt.show()
# plt.savefig(f"Difference_M{m}.png", bbox_inches="tight")
plt.clf()

for i, control in enumerate(control_variate_list2):
    plt.plot(n_list, control, label=f"{K_list[i]}")

plt.xlabel("n", fontsize=25)
plt.ylabel("error", fontsize=25)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=22)
plt.xscale("log")
plt.title("Difference between control variate and arithmetic average", fontsize=25)
plt.show()
# plt.savefig(f"Control_M{m}.png", bbox_inches="tight")

