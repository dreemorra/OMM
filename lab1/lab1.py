import numpy as np
import re

np.set_printoptions(precision=3)

def read_matrix(fname: str) -> np.array:
    """Чтение матрицы из файла."""
    matrix = []
    with open(fname, 'r') as file_mat:
        for line in file_mat.readlines():
            matrix.append(re.split('\s+', line.strip()))
    matrix = np.array(matrix, dtype=float)
    return matrix

def apply_mat(A: np.array, x: np.array):
    """Перемножение матрицы и вектора."""
    return np.array([np.sum(x*A[i]) for i in range(len(A))])

def is_invertible(A_inv: np.array, A_cap: np.array, x: np.array, i: int) -> bool:
    """Проверяем, является ли матрица А с крышечкой обратимой; 
    если обратима, то вычисляем ее обратную матрицу"""
    A_inv = A_inv.copy()
    A_cap = A_cap.copy()
    x = x.copy()

    l = A_inv @ x
    if (l[i, 0] == 0):
        return False
    else:
        l_wave = l.copy()
        l_wave[i, 0] = -1
        l_caret = l_wave/(-l[i])

        Q = np.identity(len(A_cap))
        Q[:, i] = l_caret[:, 0]

        A_cap_inv = Q @ A_inv
        print(f"inversed A_cap matrix:\n{A_cap_inv}")
        return True


if __name__ == "__main__":
    A = read_matrix("./matrix.txt")
    A_inv = read_matrix("./matrix_A.txt")
    A_cap = read_matrix("./matrix_A_cap.txt")
    x = read_matrix("./x.txt")
    print(f"A:\n{A}")
    print(f"A inversed:\n{A_inv}")
    print(f"A with cap:\n{A_cap}")
    print(f"x vector:\n{x}")
    i = 1
    is_invertible(A_inv, A_cap, x, i)