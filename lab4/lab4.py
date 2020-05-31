import numpy as np
import math

def invert_mat(B: np.array, A: np.array, x: np.array, i: int) -> np.array:
    A_cap = A.copy()
    A_cap[:, i-1] = x

    l = B @ x
    if (l[i-1] == 0):
        return None
    else:
        l_wave = l.copy()
        l_wave[i-1] = -1
        l_caret = l_wave/(-l[i-1])

        Q = np.identity(len(A_cap))
        Q[:, i-1] = l_caret[:]

        A_cap_inv = Q @ B
        return A_cap_inv

def dual_simplex_method(A, b, c, j):
    jb = j.copy()
    Ab = A[:, jb - 1].copy()
    
    Ab_inv = np.linalg.inv(Ab)
    cb = c[jb - 1]
    y = cb.dot(Ab_inv)

    while True:
        print("- - - - - - - - - - - - - - - - - - -")
        print("Псевдоплан у =", y)

        xb = Ab_inv.dot(b)
        print(xb)
        x = np.zeros(n)
        x[jb - 1] = xb

        neg_i = [i + 1 for i in range(len(x)) if x[i] < 0]
        if len(neg_i) == 0:
            print('Оптимальный план найден.\n', x)
            return

        print("x = {0}".format(x))

        neg_i = neg_i[0]
        delta_y = Ab_inv[list(jb).index(neg_i)]
        print('x[{}] < 0,\n delta(y`) = {}'.format(neg_i, delta_y))
        
        not_b_i = [i for i in range(len(c)-1) if i not in set(jb - 1)]
        print('Для кажного небазисного j:')
        mu = dict(zip(not_b_i,[np.array(np.matmul(delta_y, A[:,i])) for i in not_b_i]))
        print("mu: ", *mu.values(), sep=' ')
        
        if all([mu_i >= 0 for mu_i in mu.values()]):
            print('Unbounded')
            return None

        not_b_i_and_neg_mu = [i for i in not_b_i if mu[i] < 0]
        sigma = [(c[i] - np.array(np.matmul(y, A[:,i]))) / mu[i] for i in not_b_i_and_neg_mu]
        print("Для каждого небазисного j при mu[j] < 0:\n sigma =", sigma)
        sigma_i = dict(zip(sigma, not_b_i_and_neg_mu))
        sigma0 = min(sigma)
        j0 = sigma_i[sigma0] + 1

        jb[list(jb).index(neg_i)] = j0
        print("jb: ", jb)
        y = y + sigma0 * delta_y

        Ab_inv = invert_mat(Ab_inv, Ab, A[:,j0 - 1], list(jb).index(j0) + 1)
        if Ab_inv is None:             
            print("Unbounded")
            return None
        Ab[:, list(jb).index(j0) - 1] = A[:, j0 - 1]
    return None


if __name__ == "__main__":
    m, n = tuple(map(int, input().split()))
    A = np.zeros((m,n))

    for i in range(m):
        A[i] = np.array([*map(float, input().split())], dtype=float)
    B = np.array([*map(float, input().split())], dtype=float)
    c = np.array([*map(float, input().split())], dtype=float)
    j = np.array([*map(int, input().split())], dtype=int)

    dual_simplex_method(A, B, c, j)