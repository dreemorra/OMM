import numpy as np
import math

def find_graph_cycle(m, n, u_b):
    '''Находит цикл в графе'''
    u_b = u_b.copy()
    rows, cols = np.zeros(m), np.zeros(n)
    for i, j in u_b:
        rows[i] += 1
        cols[j] += 1
    while True:
        is_ex = True
        for k in range(m):
            if rows[k] == 1:
                is_ex = False
                for i, j in u_b:
                    if i == k:
                        cols[j] -= 1
                        rows[i] = 0
                        u_b.remove((i, j))
                        break
        for k in range(n):
            if cols[k] == 1:
                is_ex = False
                for i, j in u_b:
                    if j == k:
                        rows[i] -= 1
                        cols[j] = 0
                        u_b.remove((i, j))
                        break
        if is_ex:
            return u_b
        if len(u_b) < 4:
            return None
    
def calculate_potentials(C, u_b):
    '''Решает систему линейных уравнений и находит потенциалы'''
    u = [math.inf]*len(C)
    v = [math.inf]*len(C[0])
    u[0] = 0
    for _ in range(len(u_b)):
        for i, j in u_b:
            if v[j] == math.inf and u[i] != math.inf:
                v[j] = C[i][j] - u[i]
                break
            elif u[i] == math.inf and v[j] != math.inf:
                u[i] = C[i][j] - v[j]
                break
    return u, v

def split_graph_cycle(i_0, j_0, cycle):
    neg, pos = set(), set()
    pos.add((i_0, j_0))
    for _ in range(len(cycle) >> 1):
        for i, j in cycle:
            if i == i_0 and j != j_0:
                neg_i, neg_j = i, j
                break
        neg.add((neg_i, neg_j))
        for i, j in cycle:
            if j == neg_j and i != neg_i:
                i_0, j_0 = i, j
                break
        pos.add((i_0, j_0))
    return neg, pos

def initial_phase(a, b):
    '''начальная фаза метода, x - базисный план перевозок, u - подмножество базисных позиций'''
    m, n = len(a), len(b)
    x = np.zeros([m, n])
    u_b = set()

    i = j = 0
    while i < m and j < n:
        u_b.add((i, j))
        if a[i] < b[j]:
            x[i][j] = a[i]
            b[j] -= a[i]
            i += 1
        else:
            x[i][j] = b[j]
            a[i] -= b[j]
            j += 1

    not_u = set((i, j) for j in range(n)
                   for i in range(m) if (i, j) not in u_b)
    if len(u_b) != m + n - 1:
        for i, j in not_u:
            u_b.add((i, j))
            if find_graph_cycle(m, n, u_b) is None:
                if len(u_b) == m + n - 1:
                    break
            else:
                u_b.remove((i, j))
        not_u.difference_update(u_b)
    return x, u_b, not_u

def potential_method(a, b, c):
    #проверяем условие баланса и балансируем, если надо
    m, n = len(a), len(b)
    d = sum(a) - sum(b)

    if d < 0:
        a.append(abs(d))
        m += 1
        c.append([0 for _ in range(n)])
    elif d > 0:
        b.append(d)
        n += 1
        for r in c:
            r.append(0)

    #начальная фаза
    X, u_b, not_ub = initial_phase(a, b)
    print("Базисный план:\n", X)
    if len(not_ub) == 0:
        return X

    print("v, u - потенциалы для любого i, j прин. u_b")
    print("u_b: ", u_b)

    #основная фаза
    while True:
        u, v = calculate_potentials(c, u_b)
        print("u: ", u)
        print("v: ", v)
        delta = math.inf
        i0, j0 = min(not_ub, key=lambda x: c[x[0]][x[1]]-u[x[0]]-v[x[1]])
        delta = c[i0][j0]-u[i0]-v[j0]
        if delta >= 0:
            print("Оптимальный план найден.")
            return X
            
        u_b.add((i0, j0))
        print("Зафиксированная (теперь базисная) поцизия: ", i0, j0)
        not_ub.remove((i0, j0))
        cycle = find_graph_cycle(m, n, u_b)
        print("Цикл в графе: ", cycle)
        neg, pos = split_graph_cycle(i0, j0, cycle)
        i_star, j_star = min(neg, key=lambda elem: X[elem[0]][elem[1]])
        theta = X[i_star][j_star]
        for elem in pos:
            X[elem[0]][elem[1]] += theta
        for elem in neg:
            X[elem[0]][elem[1]] -= theta
        print("Измененный план: \n", X)
        u_b.remove((i_star, j_star))
        not_ub.add((i_star, j_star))
        print("Текущий базис: ", u_b)
        print("- - - - - - - - - - - - - - - - -")

if __name__ == "__main__":
    m, n = map(int, input().split())
    c = [list(map(int, input().split())) for _ in range(m)]
    a = list(map(int, input().split()))
    b = list(map(int, input().split()))

    print("\nВектор предложений: ", a)
    print("Вектор заявок: ", b)
    print("стоимости перевозок из a в b:\n", c)

    result = potential_method(a, b, c).astype(int)

    print("\nОптимальный план:")
    for r in result[:m]:
        print(*r[:n])
