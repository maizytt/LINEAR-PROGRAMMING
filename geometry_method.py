import numpy as np
from linear_programming_class import *
import matplotlib.pyplot as plt
import sympy as sp

def convert_into_normal_form(P):
    # Ham muc tieu
    if not P.isMin:
        for i in range(P.numberVariable):
            P.z[i] = -np.array(P.z[i])
    

    # Rang buoc ve dau
    for i in range(P.numberVariable - P.numberNewVariable):
        if P.signVariableConstraint[i] == -1:  # xi <= 0
            # Ham muc tieu
            P.z[i] = -P.z[i]

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
        if P.signInequalityConstraint[i] == 1:  # >= bi
            P.A[i] = -P.A[i]
            P.b[i] = -P.b[i]
            P.signInequalityConstraint[i] = -1  # Gan lai dau cho RB
        elif P.signInequalityConstraint[i] == 0:  # = bi => Add new constraint
            P.numberConstraint += 1
            new_constraint = np.array([-P.A[i]])
            P.A = np.vstack((P.A, new_constraint))
            P.b = np.append(P.b, -P.b[i])

            P.signInequalityConstraint[i] = -1
            P.signInequalityConstraint = np.append(P.signInequalityConstraint, -1)



def  plot_geometric(P,interset_points):     
    x1_values = np.linspace(-10, 10, 1000)
    x2_values = np.linspace(-10, 10, 1000)
    X1, X2 = np.meshgrid(x1_values, x2_values)

    condition = []
    for i in range(P.numberConstraint):
        if P.signInequalityConstraint[i] == 1:  # >= bi
            condition.append(P.A[i][0]*X1 + P.A[i][1]*X2 >= P.b[i])
            
        elif P.signInequalityConstraint[i] == 0:  # = bi   
            condition.append(P.A[i][0]*X1 + P.A[i][1]*X2 == P.b[i])
        else:
            condition.append(P.A[i][0]*X1 + P.A[i][1]*X2 <= P.b[i])

    for i  in range(P.numberVariable):
        if P.signVariableConstraint[i] == 1:
            if i == 0:
                condition.append(X1 >= 0)
            else:
                condition.append(X2 >= 0)
        elif P.signVariableConstraint[i] == -1:
            if i == 0:
                condition.append(X1 <= 0)
            else:
                condition.append(X2 <= 0)

    intersection = np.all(condition, axis=0)

    # Vẽ miền xác định
    plt.figure(figsize=(8,6))
    for i in range(len(condition)):
        plt.imshow(condition[i], extent=(-10, 10, -10, 10), origin='lower', alpha=0.3)

    plt.imshow(intersection, extent=(-10, 10, -10, 10), origin='lower',cmap = 'plasma', alpha=0.2)
    for point in interset_points:
        plt.plot(point[0], point[1], 'ro', color = 'red')

    #Trượt hàm mục tiêu
    z_values = P.z[0]* X1 - P.z[1]* X2
    contour_levels = np.linspace(-10,10, 5)
    plt.contour(X1, X2, z_values, levels=contour_levels, colors='blue')
    
    # Đặt tiêu đề và nhãn cho trục x và trục y
    plt.title('Acceptable domain')
    plt.xlabel('x1')
    plt.ylabel('x2')

    # Hiển thị đồ thị
    plt.grid(True)
    plt.show()


def get_contraint(P, x1, x2):
    condition = []
    for i in range(P.numberConstraint):
        if P.signInequalityConstraint[i] == 1:  # >= bi
            condition.append(P.A[i][0]*x1 + P.A[i][1]*x2 - P.b[i] >= 0)
        elif P.signInequalityConstraint[i] == 0:  # = bi   
            condition.append(P.A[i][0]*x1 + P.A[i][1]*x2 - P.b[i] == 0)
        else:
            condition.append(P.A[i][0]*x1 + P.A[i][1]*x2 - P.b[i] <= 0)

    # Thêm điều kiện từ các ràng buộc về biên của biến
    for i in range(P.numberVariable):
        if P.signVariableConstraint[i] == 1:
            if i == 0:
                condition.append(x1 >= 0)
            else:
                condition.append(x2 >= 0)
        elif P.signVariableConstraint[i] == -1:
            if i == 0:
                condition.append(x1 <= 0)
            else:
                condition.append(x2 <= 0)  
    return condition

def find_intersection_point(P, constraints):
    x1, x2 = sp.symbols('x1 x2')
    interset_points = [[0,0]]
    for i in range(P.numberConstraint):   
        interset_points.append( [0,P.b[i]/P.A[i][1]]) # giao với x1 = 0
        interset_points.append([P.b[i]/P.A[i][0],0]) # giao với x2 = 0
        for j in range(1,P.numberConstraint): # giao các ràng buộc với nhau
            for j in range(1,P.numberConstraint): # giao các ràng buộc với nhau
                if i != j:
                    if P.A[i][0]==P.A[j][0] and P.A[i][1]==P.A[j][1]:
                        continue
                    else:
                        interset_points.append(np.linalg.solve([P.A[i],P.A[j]], [P.b[i],P.b[j]]))
    constraint_funcs = [sp.lambdify([x1, x2], constraint, 'numpy') for constraint in  constraints]

    intersec = []
    for point in interset_points:
        point_np = np.array(point)
        if all(constraint_func(*point_np) == True for constraint_func in constraint_funcs):
            intersec.append(point)
    return intersec

def polygon_area(points):
    n = len(points)
    area = 0.0

    for i in range(n):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % n]
        area += x1 * y2
        area -= y1 * x2

    area = abs(area) / 2.0
    return area

def get_points_in_intersection(P, constraints, xmin, xmax, ymin, ymax,  num_samples = 100000):
    x1, x2 = sp.symbols('x1 x2')
    constraint_funcs = [sp.lambdify([x1, x2], constraint, 'numpy') for constraint in constraints]


    x1_samples = np.linspace(xmin, xmax, num_samples)
    x2_samples = np.linspace(ymin, ymax, num_samples)

    points_in_intersection = []
    for i in range(num_samples):
        point = np.array([x1_samples[i], x2_samples[i]])
        if all(constraint_func(*point) == True for constraint_func in constraint_funcs):
            points_in_intersection.append(point)

    return points_in_intersection

def find_line_equation(point1, point2):
    
    x1, y1 = point1
    x2, y2 = point2
    
    if x2 - x1 != 0:
        
        return [(y1-y2), (x2 - x1), x1*(y1-y2) + y1*(x2-x1)]
    elif x1 == x2:
        return [1,0,x1]  
    elif y1 == y2:
        return[0,1,y1]
    else:
        return None
def geometric_method(P):
    x1, x2 = sp.symbols('x1 x2')
    z = P.z[0]*x1 + P.z[1]*x2

    # Lay ra cac rang buoc
    condition = get_contraint(P,x1,x2)

    # Mien chap nhan duoc
    intersection = sp.And(*condition)

    # Cac dinh cua mien chap nhan duoc
    interset_points = find_intersection_point(P, condition)

    area_intersection  = polygon_area(interset_points)
    # Tim gia tri toi uu 
    z_values = []
    for point in interset_points:
        z_values.append(P.z[0]*point[0] + P.z[1]*point[1])

    if len(interset_points) != 0:
        if P.isMin:
            optimize_value = np.min(z_values)
            index = np.where(np.array(z_values) == optimize_value)[0]
        else:
            optimize_value = np.max(z_values)
            index = np.where(np.array(z_values) == optimize_value)[0]


    # Giai bai toan QHTT
    if len(interset_points) == 0:
        print("Bai toan vo nghiem")
    elif area_intersection == np.inf:
        print("Mien chap nhan duoc khong gioi noi.") 
        # Truot z tren mien giao cac rang buoc
        f = sp.solve(z.diff(x1), x1, domain=intersection)
        if len(f) != 0:
            print(f"Gia tri toi uu la {f}")
        else:
            if P.isMin:
                if P.z[0] >= 0 and P.z[1] >= 0:
                    point_check = get_points_in_intersection(P, condition, xmin = -10**6 , xmax = interset_points[index[0]][0], ymin = -10**6 , ymax = interset_points[index[0]][1],  num_samples = 100000)
                    for point in point_check:
                        if optimize_value > P.z[0]*point[0] + P.z[1]*point[1]:
                            print(f"Gia tri toi uu min cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {-np.inf}")
                            break
                        else:
                            print(f"Nghiem cua bai toan QHTT (P) la: x1 = {interset_points[index[0]][0]} va x2 = {interset_points[index[0]][1]}")
                            print(f"Gia tri toi uu min cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {optimize_value}")
                elif P.z[0] <= 0 and P.z[1] <= 0:
                    point_check = get_points_in_intersection(P, condition, xmin = interset_points[index[0]][0] , xmax = 10**6, ymin = interset_points[index[0]][1] , ymax = 10**6,  num_samples = 100000)
                    for point in point_check:
                        if optimize_value > P.z[0]*point[0] + P.z[1]*point[1]:
                            print(f"Gia tri toi uu min cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {-np.inf}")
                            break
                        else:
                            print(f"Nghiem cua bai toan QHTT (P) la: x1 = {interset_points[index[0]][0]} va x2 = {interset_points[index[0]][1]}")
                            print(f"Gia tri toi uu min cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {optimize_value}")
                elif P.z[0] >= 0 and P.z[1] <= 0:
                    point_check = get_points_in_intersection(P, condition, xmin = -10**6 , xmax = interset_points[index[0]][0], ymin = interset_points[index[0]][1] , ymax = 10**6,  num_samples = 100000)
                    for point in point_check:
                        if optimize_value > P.z[0]*point[0] + P.z[1]*point[1]:
                            print(f"Gia tri toi uu min cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {-np.inf}")
                            break
                        else:
                            print(f"Nghiem cua bai toan QHTT (P) la: x1 = {interset_points[index[0]][0]} va x2 = {interset_points[index[0]][1]}")
                            print(f"Gia tri toi uu min cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {optimize_value}")
                elif P.z[0] <= 0 and P.z[1] >= 0:
                    point_check = get_points_in_intersection(P, condition, xmin = interset_points[index[0]][0] , xmax = 10**6, ymin = -10**6 , ymax = interset_points[index[0]][1],  num_samples = 100000)
                    for point in point_check:
                        if optimize_value > P.z[0]*point[0] + P.z[1]*point[1]:
                            print(f"Gia tri toi uu min cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {-np.inf}")
                            break
                        else:
                            print(f"Nghiem cua bai toan QHTT (P) la: x1 = {interset_points[index[0]][0]} va x2 = {interset_points[index[0]][1]}")
                            print(f"Gia tri toi uu min cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {optimize_value}")
                
                    
            else:
                if P.z[0] >= 0 and P.z[1] >= 0:
                    point_check = get_points_in_intersection(P, condition, xmin = interset_points[index[0]][0] , xmax = 10**6, ymin = interset_points[index[0]][1] , ymax = 10**6,  num_samples = 100000)
                    for point in point_check:
                        if optimize_value < P.z[0]*point[0] + P.z[1]*point[1]:
                            print(f"Gia tri toi uu max cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {np.inf}")
                            break
                        else:
                            print(f"Nghiem cua bai toan QHTT (P) la: x1 = {interset_points[index[0]][0]} va x2 = {interset_points[index[0]][1]}")
                            print(f"Gia tri toi uu max cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {optimize_value}")
                elif P.z[0] <= 0 and P.z[1] <= 0:
                    point_check = get_points_in_intersection(P, condition, xmin = -10**6 , xmax = interset_points[index[0]][0], ymin = -10**6 , ymax = interset_points[index[0]][1],  num_samples = 100000)
                    for point in point_check:
                        if optimize_value < P.z[0]*point[0] + P.z[1]*point[1]:
                            print(f"Gia tri toi uu max cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {np.inf}")
                            break
                        else:
                            print(f"Nghiem cua bai toan QHTT (P) la: x1 = {interset_points[index[0]][0]} va x2 = {interset_points[index[0]][1]}")
                            print(f"Gia tri toi uu max cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {optimize_value}")
                elif P.z[0] >= 0 and P.z[1] <= 0:
                    point_check = get_points_in_intersection(P, condition, xmin = interset_points[index[0]][0] , xmax = 10**6, ymin = -10**6 , ymax = interset_points[index[0]][1],  num_samples = 100000)
                    for point in point_check:
                        if optimize_value < P.z[0]*point[0] + P.z[1]*point[1]:
                            print(f"Gia tri toi uu max cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {np.inf}")
                            break
                        else:
                            print(f"Nghiem cua bai toan QHTT (P) la: x1 = {interset_points[index[0]][0]} va x2 = {interset_points[index[0]][1]}")
                            print(f"Gia tri toi uu max cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {optimize_value}")
                elif P.z[0] <= 0 and P.z[1] >= 0:
                    point_check = get_points_in_intersection(P, condition, xmin = -10**6 , xmax = interset_points[index[0]][0], ymin = interset_points[index[0]][1] , ymax = 10**6,  num_samples = 100000)
                    for point in point_check:
                        if optimize_value < P.z[0]*point[0] + P.z[1]*point[1]:
                            print(f"Gia tri toi uu max cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {np.inf}")
                            break
                        else:
                            print(f"Nghiem cua bai toan QHTT (P) la: x1 = {interset_points[index[0]][0]} va x2 = {interset_points[index[0]][1]}")
                            print(f"Gia tri toi uu max cua {P.z[0]}*x1 + ({P.z[1]})*x2 la {optimize_value}")


    else:
        if len(index) == 1:
            print(f"Nghiem cua bai toan QHTT (P) la: x1 = {interset_points[index[0]][0]} va x2 = {interset_points[index[0]][1]}")
            print(f"Gia tri toi uu la {optimize_value}")
        else:
            print("Bai toan QHTT (P) co vo so nghiem.")
            a, b, c = find_line_equation((interset_points[index[0]][0], interset_points[index[0]][1]), (interset_points[index[1]][0], interset_points[index[1]][1]))
            print(f"Bai toan QHTT (P) co vo so nghiem thuoc duong thang {a}*x1 + ({b})*x2 = {c}")
    plot_geometric(P, interset_points)

