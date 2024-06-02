import pandas as pd
from Functions import *
from scipy import sparse
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve
from numpy import linalg as la

beam_l = 5  # Общая длина балки
elem_n = 50  # Количество конечных элементов
elem_l = beam_l / elem_n  # Длина одного конечного элемента
node_coordinates = [i for i in np.arange(0, beam_l + elem_l, elem_l)]  # Генерация координат узлов балки

K_frame = pd.read_csv('I-beam (B) GOST 26020-83.csv')  # Импорт сортамента из .csv файла
b_type = 'Yana'  # Тип жесткости стержня
gost_row = K_frame.loc[K_frame['Type'] == b_type]  # Создание строки с жесткостями стержня из сортамента

E_matrix = matrix_e(gost_row)  # Расчет матрицы Е

K_G_matrix = k_g_matrix(E_matrix, elem_l)  # Расчет матриц К и G

k_last = k_add(gost_row, beam_l)  # Расчет добавочного слагаемого

M_i = 2 * 0.7 * np.pi * sqrt(gost_row['GIt'].values[0] * gost_row['EIz'].values[0]) / beam_l
# M_i = 30
# Вычисление первоначального М
M_ex = 2 * np.pi / (beam_l * sqrt(((1 / gost_row['EIy'].values[0]) - (1 / gost_row['GIt'].values[0])) *
                                  ((1 / gost_row['EIy'].values[0]) - (1 / gost_row['EIz'].values[0]))))
# Точное значение М

M_old = [M_i]   # Массив моментов по итерациям

while True:
    print(M_old[-1])
    for i in range(len(K_G_matrix)):
        for j in range(len(K_G_matrix[i])):
            for k in range(len(K_G_matrix[i][j])):
                K_G_matrix[i][j][k] = K_G_matrix[i][j][k].subs(M, M_i)  # Подставляем М в матрицы К и G для элементов
    #print(K_G_matrix)

    k_global = sparse.dok_matrix((3 * (elem_n + 1), (3 * (elem_n + 1))))  # Матрица жесткости пустая
    g_global = sparse.dok_matrix((3 * (elem_n + 1), (3 * (elem_n + 1))))  # Матрица геометрической жесткости пустая

    k_elem_sym = 0.5 * (np.asarray(K_G_matrix[0]) + np.asarray(K_G_matrix[0]).T)  # Симметризация матрицы К для элемента
    g_elem_sym = 0.5 * (np.asarray(K_G_matrix[1]) + np.asarray(K_G_matrix[1]).T)  # Симметризация матрицы G для элемента

    for i in range(elem_n):
        k_global[int(i) * 3:int(i + 1) * 3 + 3, int(i) * 3:int(i + 1) * 3 + 3] += np.asarray(
            K_G_matrix[0])  # Сборка К глобальной
        g_global[int(i) * 3:int(i + 1) * 3 + 3, int(i) * 3:int(i + 1) * 3 + 3] += np.asarray(
            K_G_matrix[1])  # Сборка G глобальной

    k_global = k_global[3:, 3:]
    k_global[k_global.shape[0] - 1, k_global.shape[0] - 1] += k_last[0].subs(M, M_i)
    k_global[k_global.shape[0] - 2, k_global.shape[0] - 2] += k_last[0].subs(M, M_i)
    # Удаление первых трех столбцов и строк, добаление добавочной матрицы К_М

    g_global = g_global[3:, 3:]
    g_global[g_global.shape[0] - 1, g_global.shape[0] - 1] += k_last[1].subs(M, M_i)
    g_global[g_global.shape[0] - 2, g_global.shape[0] - 2] += k_last[1].subs(M, M_i)
    # Удаление первых трех столбцов и строк, добаление добавочной матрицы G_М

    k_global = k_global.tocsc()
    g_global = g_global.tocsc()
    # Переброс матриц в другой sparse тип

    x_0 = sparse.dok_matrix((k_global.shape[0], 1))
    x_0[1, 0] = 1
    x_0 = x_0.tocsc()
    # Запись столбца начального приближения

    b_matrix = sparse.csc_matrix(g_global.dot(x_0))
    # Вычисление матрицы В
    y_i = sparse.csc_matrix(spsolve(k_global, b_matrix).reshape((k_global.shape[0], 1)))
    # Вычисление y_i
    x_i = (y_i / sparse.linalg.norm(y_i)).tocsc()
    # Вычисление i-ого приближения
    lam_old = [0]
    # Значения лямбды по итерации; используется для сравнения результатов

    while True:
        # Цикл для вычисления лямбды
        b_matrix = (g_global.dot(x_i)).tocsc()
        y_i = sparse.csc_matrix(spsolve(k_global, b_matrix).reshape((k_global.shape[0], 1)))
        lam_i = (sparse.csc_matrix.transpose(x_i).dot(y_i))[0, 0]
        x_i = (y_i / sparse.linalg.norm(y_i)).tocsc()

        if abs((lam_i - lam_old[-1]) / lam_i) <= 0.000001:
            #ins = (inv(k_global).dot(g_global)).toarray()
            #print(la.eigvals(ins))
            #insa = (inv(g_global).dot(k_global)).toarray()
            #print(la.eigvals(insa))
            #print('\n')
            #print(k_global.toarray())
            #print(g_global.toarray())
           # print(lam_i)
            lam_old.append(lam_i)   # Добавление лямбды в общий список по итерации
            break
        else:
            lam_old.append(lam_i)   # Добавление лямбды в общий список по итерации

    M_i = M_i + (-1/lam_i)  # Вычисление нового момента

    if abs((M_i - M_old[-1]) / M_i) <= 0.001:
        M_old.append(M_i)  # Добавление момента в общий список по i-ой итерации
        break
    else:
        M_old.append(M_i)  # Добавление момента в общий список по i-ой итерации

print('Полученный момент: ' + str(M_i))
print('Аналитика: ' + str(M_ex))
print('Количество итераций: ' + str(len(M_old)))
print(M_old)

    # ins = (inv(k_global).dot(g_global)).toarray()
    # print(la.eigvals(ins))
