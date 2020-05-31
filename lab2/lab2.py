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

def main_phase(c: np.array, A: np.array, Ab_inv: np.array, x: np.array, j: np.array):
    Ab = np.zeros((m, len(j)))
    # формируем базисную матрицу Ab 
    for k in range(len(j)):
        Ab[:, k] = A[:, j[k]-1]

    # находим обратную матрицу
    if Ab_inv is None:
        Ab_inv = np.linalg.inv(Ab)
    print('\nAb_inv:\n', Ab_inv)
        
    while True:
        # из с выбираем элементы с соотв. номерами в j
        cb = np.array([c[i - 1] for i in j])
        print("\ncb:\n", cb)
        # вектор потенциалов
        u = np.array(cb @ Ab_inv)
        print("\nu:\n", u)
        # вектор оценок
        delta = np.array(u @ A) - c
        print("\ndelta:\n", delta)

        # проверяем, является ли план оптимальным(хотя бы одна небазисная компонента отрицательна); иначе ищем небазисную компоненту с наименьшим индексом
        j0 = -1    #здесь будет храниться её индекс
        delta_min = math.inf
        for i in range(len(c)):
            if i + 1 not in j:
                if delta[i] < 0 and delta[i] < delta_min:
                    j0 = i
                    delta_min = delta[i]
        if j0 == -1:
            print("\nBounded")
            print("Optimal x:", x)
            return 
        print("\nj0:\n", j0)

        z = Ab_inv @ A[:,j0]
        print("\nz:\n", z)
        # проверяем, ограничены ли оптимальные планы (хотя бы одна theta не inf)
        thetas = [float(x[j[i] - 1] / z[i]) if z[i] > 0 else math.inf for i in range(len(z))]
        print("\nthetas:\n", thetas)

        theta0 = min(thetas)
        print("\nmin theta:\n", theta0)

        if theta0 == math.inf:
            print("\nUnbounded")
            return

        #берем индекс theta0 и обновляем в этом месте базиса элемент, подставляя j0
        j_min = thetas.index(theta0)
        j[j_min] = j0 + 1
        
        #обновляем x для следующей итерации
        x_new = np.zeros(len(x), dtype=float)
        x_new[j - 1] = x[j - 1] - theta0 * z
        x_new[j0] = theta0
        x = x_new
        print("\nCurrent x:\n", x_new)
        print("--------------------------------------------")

        #вычисляем обратную матрицу для следующей итерации
        Ab_inv = invert_mat(Ab_inv, Ab, A[:,j0], j_min + 1)
        if Ab_inv is None:              #если она не обратима, то мы не сможем никак дальше найти оптимальный план
            print("\nUnbounded")
            return


if __name__ == "__main__":
    m, n = tuple(map(int, input().split()))
    A = np.zeros((m,n))

    for i in range(m):
        A[i] = np.array([*map(float, input().split())], dtype=float)
    B = np.array([*map(float, input().split())], dtype=float)
    c = np.array([*map(float, input().split())], dtype=float)
    x = np.array([*map(float, input().split())], dtype=float)
    j = np.array([*map(int, input().split())])

    main_phase(c, A, None, x, j)