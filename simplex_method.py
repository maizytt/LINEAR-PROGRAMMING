from linear_programming_class import *
import numpy as np


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
            
def choose_algorithm(PStand):
    bland = False
    for i in range(PStand.numberConstraint):
        if PStand.b[i] < 0:  # 2 Pha
            return 2

        if PStand.b[i] == 0:  # Bland
            bland = True
    return 1 if bland else 0  # Don hinh
def convert_into_table_form(P, Table):
    # Đưa hàm mục tiêu vào bảng
    for j in range(P.numberVariable):
        Table[0][j] = P.z[j]

    # Đưa các ràng buộc dạng thức vào bảng
    for i in range(1, 1 + P.numberConstraint):
        # Đối với các biến ban đầu (Biến không cơ sở) Xj
        for j in range(P.numberVariable):
            Table[i][j] = P.A[i - 1][j]
        # Đối với các biến cơ sở Wj
        for j in range(P.numberVariable, P.numberVariable + P.numberConstraint):
            if i - 1 == j - P.numberVariable:
                Table[i][j] = 1
        # Đặt giá trị b
        Table[i][P.numberVariable + P.numberConstraint] = P.b[i - 1]
def convert_into_table_form(p, table):
    table[0, :p.numberVariable] = p.z
    for i in range(1, p.numberConstraint + 1):
        table[i, :p.numberVariable] = p.A[i - 1]
        table[i, p.numberVariable + i - 1] = 1
        table[i, -1] = p.b[i - 1]
def find_arg_min_ratio(table, num_rows, num_cols, y_pivot, phase1):
    x_pivot = -1
    min_ratio = -1
    ratio = 0

    for i in range(1, num_rows):
        if table[i][y_pivot] > 0:
            min_ratio = table[i][num_cols - 1] / table[i][y_pivot]
            x_pivot = i
            break

    if x_pivot == -1:  # Co bien vao nhung khong co bien ra -> bai toan khong gioi noi
        return -1

    for i in range(i + 1, num_rows):
        if table[i][y_pivot] > 0:
            ratio = table[i][num_cols - 1] / table[i][y_pivot]
            if ratio < min_ratio:
                min_ratio = ratio
                x_pivot = i
            if phase1:
                if ratio == min_ratio and table[i][num_cols - 2] == 1:
                    x_pivot = i

    return x_pivot
def choose_pivot_dantzig(table, num_rows, num_cols, phase1):
    min_c = 0
    y_pivot = -1

    for i in range(num_cols - 1):
        if table[0][i] < 0 and table[0][i] < min_c:
            min_c = table[0][i]
            y_pivot = i
    
    if y_pivot == -1:  # Ham muc tieu da thoa dieu kien dung
        return 0, -1, -1

    # Bien vao la bien o cot y_pivot
    x_pivot = find_arg_min_ratio(table, num_rows, num_cols, y_pivot, phase1)

    if x_pivot == -1:
        return -1, -1, -1  # Co bien vao nhung khong co bien ra -> infity

    return 1, x_pivot, y_pivot

def choose_pivot_bland(table, num_rows, num_cols):
    y_pivot = -1

    for i in range(num_cols - 1):
        if table[0][i] < 0:
            y_pivot = i
            break

    if y_pivot == -1:  # Da thoa dieu kien dung
        return 0, -1, -1

    x_pivot = find_arg_min_ratio(table, num_rows, num_cols, y_pivot, False)

    if x_pivot == -1:
        return -1, -1, -1  # Co bien vao nhung khong co bien ra -> infity

    return 1, x_pivot, y_pivot
def rotate_pivot(Table, num_rows, num_cols, xPivot, yPivot):
    for i in range(num_rows):
        if( i!= xPivot):
            coef = -Table[i][yPivot] / Table[xPivot][yPivot]
            for j in range(num_cols):
                Table[i][j] += coef * Table[xPivot][j]
        else:
            coef = Table[xPivot][yPivot]
            for j in range(num_cols):
                Table[xPivot][j] /= coef


def find_pivot_of_column(table, row_start, row_end, col_start, col_end):
    for j in range(col_start, col_end):
        pivot_found = False
        for i in range(row_start, row_end):
            if table[i][j] != 0:
                pivot_found = True
                break
        if pivot_found:
            return j
    return -1
def find_pivot_of_column_2_phase(Table, col):
    xPivot = -1
    flag = False
    for i in range(1, Table.shape[0]):
        if Table[i, col] == 0:
            continue

        if Table[i, col] == 1:
            if flag == False:
                xPivot = i
                flag = True
            else:
                return -1
        else:
            return -1
        
    return xPivot
def dantzig(table, num_rows, num_cols, phase1 = 0):
    x_pivot, y_pivot = -1, -1
    count = 0
    while True:
        check, x_pivot, y_pivot = choose_pivot_dantzig(table, num_rows, num_cols, phase1)
        count+=1
        if check == 1:
            rotate_pivot(table, num_rows, num_cols, x_pivot, y_pivot)
    
        else:  # check = -1: co bien vao khong co bien ra, 0: dung
            return -check, count  # 1: khong gioi noi, 0: dung
    check = 0 
    return check, count

def bland(table, num_rows, num_cols):
    x_pivot, y_pivot = -1, -1
    count = 0
    while True:
        check, x_pivot, y_pivot = choose_pivot_bland(table, num_rows, num_cols)
        count+=1
        if check != 1:  # check = -1: co bien vao khong co bien ra, 0: dung
            return -check, count  # 1: khong gioi noi, 0: dung
        else:
            rotate_pivot(table, num_rows, num_cols, x_pivot, y_pivot)
    return 0, count
def two_phase(table, num_rows, num_cols):
    # Lap bang moi
    table_p1 = [[0.0] * (num_cols + 1) for _ in range(num_rows)]
    table_p1[0][num_cols - 1] = 1.0
    table_p1[0][num_cols] = table[0][num_cols - 1]
    for i in range(1, num_rows):
        for j in range(num_cols - 1):
            table_p1[i][j] = table[i][j]
        table_p1[i][num_cols - 1] = -1.0
        table_p1[i][num_cols] = table[i][num_cols - 1]

    x_pivot, y_pivot = -1, num_cols - 1
    min_b = 0
    for i in range(1, num_rows):
        if table[i][y_pivot] < min_b:
            min_b = table[i][y_pivot]
            x_pivot = i

    rotate_pivot(table_p1, num_rows, num_cols + 1, x_pivot, y_pivot)
    check, phase1_loop = dantzig(table_p1, num_rows, num_cols + 1, 1)
    print('So lan lap o pha 1 la: ',phase1_loop)

    # Da dat "tu vung" toi uu
    for j in range(num_cols - 1):
        if table_p1[0][j] != 0:
            return -1  # Bai toan vo nghiem

    # Chuyen sang Pha 2    
    for i in range(1, num_rows):
        for j in range(num_cols - 1):
            table[i][j] = table_p1[i][j]
        table[i][num_cols - 1] = table_p1[i][num_cols]

    for j in range(num_cols):
        x_pivot = find_pivot_of_column_2_phase(table, j)
        if x_pivot == -1:
            continue
        rotate_pivot(table, num_rows, num_cols, x_pivot, j)
    
    
    # Da thu duoc tu vung moi
    # Ap dung Dantzig de giai
    check, phase2_loop = dantzig(table, num_rows, num_cols,0)
    print('So lan lap o pha 2 la: ',phase2_loop)
    return check
def check_one_root(problem, table, num_cols, pivots):
    for i in range(num_cols - 1):
        if i >= problem.numberVariable - problem.numberNewVariable and i < problem.numberVariable:
            continue
        if pivots[i] == -1 and abs(table[0][i]) < 1e-4 and problem.signVariableConstraint[i] != 0:
            return False
    return True

def find_name_variable(problem, table, num_cols, index):
    name = ""
    if index < problem.numberVariable - problem.numberNewVariable:
        name = "x" + str(index + 1)
        return 1, name
    elif index + 1 > problem.numberVariable and index + 1 < num_cols:
        name = "w" + str(index + 1 - problem.numberVariable)
        return 1, name
    return 0, name

def output_of_table(problem, table, num_rows, num_cols, result):
    print("\n", "-" * 15, " KET QUA BAI TOAN ", "-" * 15)

    if result == 1:  # Khong gioi noi
        if problem.isMin:
            print("   Bai toan Khong gioi noi. \n   Tuc la Min z = -∞")
        else:
            print("   Bai toan Khong gioi noi. \n   Tuc la Max z = +∞")
        return

    elif result == 0:  # Bai toan dat dieu kien dung: Nghiem duy nhat hoac VSN
        # In GTTU
        if problem.isMin:
            print(f"   Min z = {-table[0][num_cols - 1]}\n")
        else:
            print(f"   Max z = {table[0][num_cols - 1]}\n")

        # In nghiem toi uu
        pivots = [find_pivot_of_column(table, 1, num_rows, j, j + 1) for j in range(num_cols - 1)]

        if check_one_root(problem, table, num_cols, pivots):  # Nghiem duy nhat
            print("   Nghiem toi uu la: \n")
            for j in range(problem.numberVariable - problem.numberNewVariable):
                if table[0][j] != 0:
                    print(f"   x{j + 1} = {0:.4f}")
                    continue
                count = 0
                index = 0
                for i in range(1, num_rows):
                    if table[i][j] != 0:
                        count += 1
                        index = i
                temp = round(table[index][num_cols - 1] / table[index][j] * 1e4) / 1e4
                if problem.signVariableConstraint[j] == -1:
                    print(f"   x{j + 1} = {-temp:.4f}")
                else:
                    print(f"   x{j + 1} = {temp:.4f}")
            
        else:  # VSN
            print("   Bai toan co vo so nghiem.")
            print("   Nghiem toi uu co dang: \n")

            signs = [1] * (num_cols - 1)
            for i in range(problem.numberVariable - problem.numberNewVariable):
                if problem.signVariableConstraint[i] < 0:
                    signs[i] = -1

            for i in range(problem.numberVariable - problem.numberNewVariable):
                if pivots[i] == -1:
                    if abs(table[0][i]) > 1e-4:
                        print(f"   x{i + 1} = 0")
                    else:
                        if problem.signVariableConstraint[i] == 0:
                            print(f"   x{i + 1} tu do")
                        else:
                            print(f"   x{i + 1} {problem.convert_int_to_inequality_sign(problem.signVariableConstraint[i])} 0 ")
                else:
                    print(f"   x{i + 1} = {signs[i] * table[pivots[i]][num_cols - 1]} ", end="")
                    for j in range(num_cols - 1):
                        if abs(table[0][j]) > 1e-4 or pivots[j] != -1 or j == i:
                            continue
                        status, name_var = find_name_variable(problem, table, num_cols, j)
                        if status == 0:  # Bien moi do bien tu do sinh ra thi khong can xuat
                            continue
                        else:
                            if -signs[i] * signs[j] * table[pivots[i]][j] == 0:
                                continue
                            if -signs[i] * signs[j] * table[pivots[i]][j] > 0:
                                print(" + ", end="")
                            print(f"{-signs[i] * signs[j] * table[pivots[i]][j]}{name_var} ", end="")
                    print()

            print("\n   Thoa (cac) dieu kien:\n")
            for i in range(num_cols - 1):
                if i >= problem.numberVariable - problem.numberNewVariable and i < problem.numberVariable:
                    continue
                if i < problem.numberVariable - problem.numberNewVariable and problem.signVariableConstraint[i] == 0:
                    continue
                if pivots[i] == -1:
                    if i < problem.numberVariable - problem.numberNewVariable:
                        continue
                    if abs(table[0][i]) < 1e-4:
                        status, name_var = find_name_variable(problem, table, num_cols, i)
                        if status == 1:
                            print("   ", end="")
                            if i >= problem.numberVariable - problem.numberNewVariable:
                                print(f"{name_var} >= 0 ")
                            else:
                                if problem.signVariableConstraint[i] == 0:
                                    print(f"{name_var} tu do")
                                else:
                                    print(f"{name_var} {problem.convert_int_to_inequality_sign(problem.signVariableConstraint[i])} 0 ")
                else:
                    print(f"   {signs[i] * table[pivots[i]][num_cols - 1]} ", end="")
                    for j in range(num_cols - 1):
                        if abs(table[0][j]) > 1e-4 or pivots[j] != -1:
                            continue  # Neu la bien phu thuoc hoac bien nam tren ham muc tieu (=0) khi ko xuat
                        status, name_var = find_name_variable(problem, table, num_cols, j)
                        if status == 0:  # Bien moi do bien tu do sinh ra thi khong can xuat
                            continue
                        else:
                            if -signs[i] * signs[j] * table[pivots[i]][j] == 0:
                                continue
                            if -signs[i] * signs[j] * table[pivots[i]][j] > 0:
                                print(" + ", end="")
                            print(f"{-signs[i] * signs[j] * table[pivots[i]][j]}{name_var} ", end="")
                    print(" >= 0 ")
    else:
        print("   Bai toan vo nghiem!")

    print("\n", "-" * 48, "-")
def solve_problem(P):
    # Đưa P về dạng chuẩn
    convert_into_normal_form(P)

    # Đưa P ở dạng chuẩn thành bảng
    Table = np.zeros((P.numberConstraint + 1, P.numberVariable + P.numberConstraint + 1))

    number_row = P.numberConstraint + 1
    number_col = P.numberConstraint + P.numberVariable + 1


    convert_into_table_form(P, Table)

    # Chọn thuật toán để xoay
    x_pivot = -1
    y_pivot = -1

    if choose_algorithm(P) == 0:  # Đơn hình
        check, loop = dantzig(Table, number_row, number_col)
        print('So lan lap tu vung la: ',loop)
    elif choose_algorithm(P) == 1:  # Bland
        check, loop = bland(Table, number_row, number_col)
        print('So lan lap tu vung la: ',loop)
    else:  # 2 Pha
        check = two_phase(Table, number_row, number_col)

    output_of_table(P, Table, number_row, number_col, check)