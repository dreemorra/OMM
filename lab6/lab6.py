import numpy as np
np.set_printoptions(precision = 3)

def quad_method(c, D, A, jb, jb_star, x):
    n = len(c)
    while True:
        m = len(A)
        k = len(jb_star)

        if len(np.setdiff1d(jb, jb_star)) > 0:
            print("Unbounded")
            return

        c_x = c + D @ x
        A_b_inv = np.linalg.inv(A[:, jb - 1])
        u_x = -c_x[jb - 1] @ A_b_inv
        delta_x = u_x @ A + c_x
        
        H = np.zeros((k + m, k + m))
        H[:k, :k] = D[jb_star - 1][:, jb_star - 1]
        H[k:, :k] = A[:, jb_star - 1]
        H[:k, k:] = A[:, jb_star - 1].T
        
        if np.linalg.det(H) == 0:
            print("Unbounded")
            return

        negs = np.where(delta_x < 0)[0]
        if len(negs) == 0:
            print("Bounded")
            for i in x: print(f"{i}", end=' ')
            return
            
        j0 = negs[0] + 1

        b_star = np.zeros(k + m)
        b_star[:k] = D[jb_star - 1, j0 - 1]
        b_star[k:] = A[:, j0 - 1]
        t_x = -np.linalg.inv(H) @ b_star

        l = np.zeros(n)
        l[jb_star - 1] = t_x[:k]
        l[j0 - 1] = 1
        d = l @ D @ l

        theta_j0 = np.abs(delta_x[j0 - 1]) / d if d > 0 else np.inf
        theta = np.zeros(k)
        theta[:] = np.inf
        
        negs = np.where(l[jb_star - 1] < 0)[0]
        theta[negs] = -x[jb_star - 1][negs] / l[jb_star - 1][negs]
        
        if min(theta_j0, theta.min()) == np.inf:
            print("Unbounded")
            return
        
        j_star = j0 if theta_j0 < theta.min() else jb_star[np.where(theta == theta.min())[0][0]]
        x = x + min(theta_j0, theta.min()) * l

        if j0 == j_star:
            jb_star = np.hstack((jb_star, j_star))
        elif len(jb[jb == j_star]) == 0:
            jb_star = jb_star[jb_star != j_star]
        else:
            s = np.where(jb == j_star)[0][0]
            jp = np.nan
            for j in np.setdiff1d(jb_star, jb):
                if (A_b_inv @ A[:, j - 1])[0] != 0:
                    jp = j
                    break
            if not np.isnan(jp):
                jb[jb == j_star] = jp
                jb_star = jb_star[jb_star != j_star]
            else:
                jb[jb == j_star] = j0
                jb_star[jb_star == j_star] = j0

if __name__ == "__main__":
    m, n = map(int, input().split())
    A = np.zeros((m, n))
    D = np.zeros((n, n))
    for i in range(m):
        A[i] = np.array([*map(float, input().split())], dtype=float)
    B = np.array([*map(float, input().split())], dtype=float)
    c = np.array([*map(float, input().split())], dtype=float)
    for i in range(n):
        D[i] = np.array([*map(float, input().split())], dtype=float)
    x = np.array([*map(float, input().split())], dtype=float).T
    jb = np.array([*map(int, input().split())], dtype=int)
    jb_star = np.array([*map(int, input().split())], dtype=int)
    
    quad_method(c, D, A, jb, jb_star, x)
