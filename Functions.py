import numpy as np
from sympy import *


M, S = symbols('M S')  # Объявляем символьные переменные


def matrix_e(stiff):
    k_2z = stiff['GFy'].values[0]
    k_3y = stiff['EIz'].values[0]
    k_3z = stiff['EIy'].values[0]
    k_t = stiff['GIt'].values[0]  # Запись жесткостей

    E_matrix = np.array([[k_2z, -k_2z * sin(M * S / k_3z), 0, k_2z * cos(M * S / k_3z), 0],
                         [-k_2z * sin(M * S / k_3z), k_2z * (sin(M * S / k_3z)) ** 2, 0,
                          -k_2z * sin(M * S / k_3z) * cos(M * S / k_3z), -1 / 2 * M],
                         [0, 0, k_t * (cos(M * S / k_3z)) ** 2 + k_3y * (sin(M * S / k_3z)) ** 2, 1 / 2 * M,
                          (k_t - k_3y) * sin(M * S / k_3z) * cos(M * S / k_3z)],
                         [k_2z * cos(M * S / k_3z), -k_2z * sin(M * S / k_3z) * cos(M * S / k_3z), 1 / 2 * M,
                          k_2z * (cos(M * S / k_3z)) ** 2, 0],
                         [0, -1 / 2 * M, (k_t - k_3y) * sin(M * S / k_3z) * cos(M * S / k_3z), 0,
                          k_t * (sin(M * S / k_3z)) ** 2 + k_3y * (cos(M * S / k_3z)) ** 2]])

    return E_matrix


def k_g_matrix(e_mat, l_e):
    N_matrix = np.array([[-1/l_e, 0, 0, 1/l_e, 0, 0],
                         [0, 1-(S/l_e), 0, 0, S/l_e, 0],
                         [0, -1/l_e, 0, 0, 1/l_e, 0],
                         [0, 0, 1-(S/l_e), 0, 0, S/l_e],
                         [0, 0, -1/l_e, 0, 0, 1/l_e]])
    #print(N_matrix.T)
    #print(e_mat)
    #print(N_matrix)
    pre_int_stiff = np.dot(np.dot(N_matrix.T, e_mat), N_matrix)
    #print(pre_int_stiff)

    K_el = [[], [], [], [], [], []]
    G_el = [[], [], [], [], [], []]

    for i in range(len(pre_int_stiff)):
        for j in range(len(pre_int_stiff[i])):
            k = integrate(pre_int_stiff[i][j], (S, 0, l_e))
            g = diff(k, M)
            K_el[i].append(k)
            G_el[i].append(g)

    return K_el, G_el


def k_add(stiff, l_beam):
    k_3z = stiff['EIy'].values[0]
    k_add_matrix = (-sin(M*l_beam/k_3z))/(2*(1-cos(M*l_beam/k_3z)))
    return k_add_matrix, diff(k_add_matrix, M)
