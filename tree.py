import numpy as np

def binomial_tree(S_0, K, T, r, sigma, N):
    """
    Builds the Binomial tree for the stock price.
    :param S_0: spot price at t=0
    :param K: strike price at t=T
    :param T: Time till expiry
    :param r: risk-free interest rate
    :param sigma: historic volatility
    :param N: steps in the tree
    :return: Lower-Triangular matrix representation of binomial tree containing stock prices
    """
    dt = T/N
    tree = np.zeros((N+1,N+1))

    u = np.exp(sigma*np.sqrt(dt))
    d = np.exp(-sigma*np.sqrt(dt))
    for i in range(N+1):
        for j in range(i+1):
            tree[i, j] = S_0*((u**j)*(d**(i-j)))
    return tree

def put_price_tree(tree, K, T, r, sigma, N):
    """
    Calculates the put option price in binomial tree using backwards induction.
    :param tree: Binomial tree containing stock prices
    :param K: strike price at t=T
    :param T: Time till expiry
    :param r: risk-free interest rate
    :param sigma: historic volatility
    :param N: steps in the tree
    :returns: Lower-Triangular matrix representation of binomial tree containing put option prices
    """
    dt = T/N

    # up and down
    u = np.exp(sigma*np.sqrt(dt))
    d = np.exp(-sigma*np.sqrt(dt))

    # risk free probability
    q = (np.exp(r*dt)-d)/(u-d)

    # Binomial tree dimensions
    col = tree.shape[1]
    row = tree.shape[0]

    # Determines put price in the last row of the binomial tree
    for i in range(col):
        S_T = tree[row-1, i]
        tree[row-1, i] = max(0, K-S_T)

    # Backwards induction
    for i in range(row-1)[::-1]:
        for j in range(i+1):
            price_down = tree[i+1, j]
            price_up = tree[i+1, j+1]
            tree[i,j] = np.exp(-r*dt)*(q*price_up+(1-q)*price_down)

    return tree
