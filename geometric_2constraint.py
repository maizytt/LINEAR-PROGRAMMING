from linear_programming_class import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from fractions import Fraction
import sympy as sp
import warnings
warnings.filterwarnings("ignore")

def convert_into_canonical_form(P):
    # Ham muc tieu
    if not P.isMin:
        for i in range(P.numberVariable):
            P.z[i] = -np.array(P.z[i])
    
    # Rang buoc ve dau
    for i in range(P.numberVariable - P.numberNewVariable):
        if P.signVariableConstraint[i] == -1:  # xi <= 0
            # Ham muc tieu
            P.z[i] = -np.array(P.z[i])

            # Rang buoc dang thuc
            for j in range(P.numberConstraint):
                P.A[j][i] = -P.A[j][i]
            P.signVariableConstraint[i] = 1

        elif P.signVariableConstraint[i] == 0:  # xi tuy y => Add column xi-
            P.numberVariable += 1
            P.numberNewVariable += 1
            # Ham muc tieu
            P.z = np.append(P.z, -P.z[i])
            P.signVariableConstraint[i] = 1
            P.signVariableConstraint.append(1)
            # Rang buoc dang thuc
            P.A = np.hstack((P.A, -P.A[:, i:i+1]))

    # Rang buoc dang thuc
    for i in range(P.numberConstraint):
        
        if P.signInequalityConstraint[i] == -1:  # <= bi
            
            P.numberVariable += 1
            P.numberNewVariable += 1
            P.signVariableConstraint[i] = 1
            P.signVariableConstraint.append(1)
            new_column = np.zeros((P.A.shape[0], 1))
            new_column[i, 0] = 1
            P.A = np.hstack((P.A, new_column))
            P.signInequalityConstraint[i] = 0  # Gan lai dau cho RB
            P_z_list = P.z.tolist()
            P_z_list.append(0)
            P.z = np.array(P_z_list)

        elif P.signInequalityConstraint[i] == 1:  # >= bi 
            
            P.numberVariable += 1
            P.numberNewVariable += 1
            P.signVariableConstraint[i] = 1
            P.signVariableConstraint.append(1)
            new_column = np.zeros((P.A.shape[0], 1))
            new_column[i, 0] = -1
            P.A = np.hstack((P.A, new_column))
            P.signInequalityConstraint[i] = 0  # Gan lai dau cho RB
            P_z_list = P.z.tolist()
            P_z_list.append(0)
            P.z = np.array(P_z_list)
            

def plot_geometric_2constraint(P,indices_inside_hull, X):
    x_values = [point[0] for point in X[indices_inside_hull]]
    y_values = [point[1] for point in X[indices_inside_hull]]

    plt.fill(x_values + [x_values[0]], y_values + [y_values[0]], 'skyblue')
    plt.axvline(x=P.b[0], color='b', linestyle='--', label= f'x = ({P.b[0]}, z)')
    for i, (x, y) in enumerate(zip(x_values, y_values), start=0):
        plt.scatter(x, y, color='black', s=30)
        plt.annotate(f'A{indices_inside_hull[i]+1}', (x, y), xytext=(5, 5), textcoords='offset points')

    # Ghi nhãn trục x, y
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Geometric metod with 2 constraints')

    # Hiển thị đồ thị
    plt.legend()
    plt.show()

def geometric_2constraint(P):
    
    convert_into_canonical_form(P)

    x1, x2 = sp.symbols('x1 x2')
    z = P.z[0]*x1 + P.z[1]*x2

        # Tính các hệ số áp dụng tính chất tập lồi
    sum_constraint = np.sum(P.A, axis=0) / np.sum(P.b)

    X = np.zeros((P.numberVariable,2))
    for i in range(P.numberVariable):
        X[i][0] = P.A[0][i]/sum_constraint[i] # Lấy 1 trong 2 ràng buộc
        X[i][1] = P.z[i]/sum_constraint[i] # Lấy z
    y = [P.b[0], z]


    # Tìm bao lồi của tập điểm
    hull = ConvexHull(X)

    # Lấy các chỉ số của các điểm nằm trong bao lồi
    indices_inside_hull = hull.vertices

    # In các điểm nằm trong bao lồi
    points_inside_hull = X[indices_inside_hull]

    z_values = []
    index = []

    for i in indices_inside_hull:
        for j in indices_inside_hull:
            if i != j :
                if X[i][0] > X[j][0] and X[j][0] <= P.b[0] and P.b[0] <= X[i][0]:
                    coef = np.linalg.solve([[X[i][0],1], [X[j][0], 1]], [X[i][1],X[j][1]])
                    z_values.append(P.b[0]*coef[0] + coef[1])
                    index.append([i, j])
                elif X[i][0] < X[j][0] and X[j][0] >= P.b[0] and P.b[0] >= X[i][0]:
                    coef = np.linalg.solve([[X[i][0],1], [X[j][0], 1]], [X[i][1],X[j][1]])
                    z_values.append(P.b[0]*coef[0] + coef[1])
                    index.append([i, j])


    # Tìm các nghiệm của bài toán
    solution = np.zeros(P.numberVariable)
    if P.isMin:
        optimize_value = np.min(z_values)
        x1, x2 = index[np.argmin(z_values)]
        temp = np.linalg.solve([[1,1], [X[x1][0], X[x2][0]]], [1, P.b[0]])
        solution[x1] = Fraction(temp[0]/sum_constraint[x1]).limit_denominator()
        solution[x2] = Fraction(temp[1]/sum_constraint[x2]).limit_denominator()
    else:
        optimize_value = -np.min(z_values)
        x1, x2 = index[np.argmin(z_values)]
        temp = np.linalg.solve([[1,1], [X[x1][0], X[x2][0]]], [1, P.b[0]])
        solution[x1] = Fraction(temp[0]/sum_constraint[x1]).limit_denominator()
        solution[x2] = Fraction(temp[1]/sum_constraint[x2]).limit_denominator()


    print("Các nghiệm của bài toán QHTT (P) :")
    for i in range(len(solution)):
        print(f"x{i+1} =", solution[i])

    print(f"\nGiá trị tối ưu {'Min' if P.isMin else 'Max'} là ", optimize_value)


    plot_geometric_2constraint(P,indices_inside_hull,X)


