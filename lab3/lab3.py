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

    for k in range(len(j)):
        Ab[:, k] = A[:, j[k]-1]

    if Ab_inv is None:
        Ab_inv = np.linalg.inv(Ab)
        
    while True:
        cb = np.array([c[i - 1] for i in j])
        u = np.array(cb @ Ab_inv)
        delta = np.array(u @ A) - c

        j0 = -1    
        delta_min = math.inf
        for i in range(len(c)):
            if i + 1 not in j:
                if delta[i] < 0 and delta[i] < delta_min:
                    j0 = i
                    delta_min = delta[i]
        if j0 == -1:
            print("Bounded")
            return x, j

        z = Ab_inv @ A[:,j0]
        
        thetas = [float(x[j[i] - 1] / z[i]) if z[i] > 0 else math.inf for i in range(len(z))]
        theta0 = min(thetas)
        if theta0 == math.inf:
            print("Unbounded")
            return None, None

        j_min = thetas.index(theta0)
        j[j_min] = j0 + 1
        
        x_new = np.zeros(len(x), dtype=float)
        x_new[j - 1] = x[j - 1] - theta0 * z
        x_new[j0] = theta0
        x = x_new

        Ab_inv = invert_mat(Ab_inv, Ab, A[:,j0], j_min + 1)
        if Ab_inv is None:             
            print("Unbounded")
            return None, None

def initial_phase(A: np.array, b: np.array, c: np.array) -> np.array:
    # приводим к неотрицательному виду все b
    for index, b_item in enumerate(b):
        if b_item < 0:
            b[index] = b_item * (-1)
            A[:, index] *= -1
    
    x_temp = np.concatenate((np.zeros(len(c)), b))
    A_temp = np.concatenate((A, np.eye(len(A))), axis=1)
    c_temp = np.concatenate((np.zeros(len(c)), np.array([-1 for _ in range(len(A))])))
    Jb_temp = np.array([i for i in range(len(c) + 1, len(c) + len(A) + 1)])

    # применяем основную фазу
    x_new, Jb_new = main_phase(c_temp, A_temp, None, x_temp, Jb_temp)

    # базисный допустимый план для основной задачи 
    x_result = x_new[:len(c)]
    
    # проверяем, совместна ли задача
    if any(elem != 0 for elem in x_new[len(c):]):
        print('No solution')
        return
    
    # ко всем индексам из Jb, которые не являются индексами родных переменных, применяем корректирующий алгоритм
    l = []
    j_donot_fit = [(index, j) for index, j in enumerate(Jb_new) if j not in range(len(c))]
    for k in j_donot_fit:
        #для каждого небазисного индекса j родной переменной вычисляем вектор l
        i = k[1] - len(c) -1
        j_not_in_b = [j for j in range(len(c)) if j not in Jb_new]
        Ab_inv = np.linalg.inv(A_temp.copy()[:, Jb_new-1])
        
        for j_not_in_b_item in j_not_in_b:
            l.append(int(np.matmul(Ab_inv, A_temp[:, j_not_in_b_item])[i]))
    
        # уравнение i избыточное в A, A_temp и b, а элемент k - в Jb
        if all(l_item == 0 for l_item in l):
            A = np.delete(A, (i), axis=0)
            A_temp = np.delete(A_temp, (i), axis=0)
            b = np.delete(b, i)
            Jb_new = np.delete(Jb_new, k[0])
        #иначе находим индекс ненулевого j и заменяем jk на j
        else:
            target_i = next((l_item for l_item in l if l_item != 0), None)
            Jb_new[k[0]] = target_i
    
    print("x:", x_result)
    print("Jb:", Jb_new)
    print("A:\n", A)
    print("b:", b)

if __name__ == "__main__":
    m, n = tuple(map(int, input().split()))
    A = np.zeros((m,n))

    for i in range(m):
        A[i] = np.array([*map(float, input().split())], dtype=float)
    B = np.array([*map(float, input().split())], dtype=float)
    c = np.array([*map(float, input().split())], dtype=float)

    initial_phase(A, B, c)