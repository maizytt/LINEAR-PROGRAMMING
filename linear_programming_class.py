import numpy as np
import copy

class Problem:
    def __init__(self):
        self.numberVariable = 0
        self.numberConstraint = 0
        self.z = []
        self.isMin = False
        self.A = []
        self.b = []
        self.signInequalityConstraint = []
        self.signVariableConstraint = []
        self.numberNewVariable = 0

    def convert_inequality_sign_to_int(self, sign):
        if sign == "<=":
            return -1
        elif sign == "=":
            return 0
        else:
            return 1

    def convert_int_to_inequality_sign(self, sign):
        if sign == -1:
            return "<="
        elif sign == 0:
            return "="
        else:
            return ">="

    def input_problem(self):
        print("\n>>>>>>>>>>>>>>>>\t\tNHAP BAI TOAN\t\t<<<<<<<<<<<<<<<<<<<<<\n")
        self.numberNewVariable = 0
        self.numberVariable = int(input("So luong bien: "))
        self.numberConstraint = int(input("So rang buoc (loai tru rang buoc ve dau cua cac bien): "))
        
        
        print("\n\t\t------------\tHAM MUC TIEU\t---------------")
        min_max = input("\nNeu tim Gia tri nho nhat - min, vui long nhap \"Min\" nguoc lai nhap \"Max\" (khong phan biet chu hoa - chu thuong): ").lower()
        self.isMin = (min_max == "min")
        print("\nNhap he so ham muc tieu tuong ung voi (x1, x2, ..., xn), cac he so cach nhau boi <khoang trang - space>:")
        self.z = np.array([float(x) for x in input("Nhap: ").split()])
        
        self.A = np.zeros((self.numberConstraint, self.numberVariable))
        self.b = np.zeros(self.numberConstraint)
        self.signInequalityConstraint = np.zeros(self.numberConstraint, dtype=int)
        
        print("\n\t\t---------------\tRANG BUOC (BAT) DANG THUC\t---------------\n")
        for j in range(self.numberConstraint):
            print(f"\nNhap cac he so tuong ung voi (x1, x2, ..., xn) cua RB{j + 1}, cac he so cach nhau boi <khoang cach>: ", end="")
            self.A[j] = np.array([float(x) for x in input().split()])
            tempS = input("Nhap dau cua rang buoc (>=, =, <=): ")
            self.signInequalityConstraint[j] = self.convert_inequality_sign_to_int(tempS)
            self.b[j] = float(input(f"Nhap b{j + 1}: "))
        
        print("\t\t---------------RANG BUOC VE DAU CUA BIEN---------------\n")
        print("Nhap lan luot cac rang buoc ve dau cua cac bien (x1, x2, ..., xn).")
        print("Nhap 1 neu bien >= 0, -1 neu bien <= 0, 0 neu bien tu do, cac he so cach nhau boi <khoang cach>:")

        self.signVariableConstraint = [int(x) for x in input("Nhap: ").split()]
        print("\n------------------------------------------------------------------------------\n")
            
            
    def print_problem(self):
        print("\n" + ">" * 30 + " BAI TOAN QUY HOACH TUYEN TINH " + "<" * 30 + "\n")
        print(f"{'Min z =' if self.isMin else 'Max z ='} ", end="")
        for i, coef in enumerate(self.z):
            temp = f"{'+ ' if coef >= 1 and i > 0 else ''}{coef}*x{i + 1}"
            print(f"{temp:>9}", end="")
        print("\n")
        
        for j in range(self.numberConstraint):
            for i, coef in enumerate(self.A[j]):
                temp = f"{'+ ' if coef >= 1 and i > 0 else ''}{coef}*x{i + 1}"
                print(f"{temp:>9}", end="")
            print(f" {self.convert_int_to_inequality_sign(self.signInequalityConstraint[j]):>4}{self.b[j]:>6.2f} ({j + 1}),")
        
        for i, sign in enumerate(self.signVariableConstraint):
            if sign == 1:
                print(f"{'':>15} x{i + 1}{'>=     0':>22}")
            elif sign == -1:
                print(f"{'':>15} x{i + 1}{'<=     0':>22}")
        
        print("\n" + "-" * 61 + "\n")

    def convert_into_normal_form(self):
        # Ham muc tieu
        if not self.isMin:
            for i in range(self.numberVariable):
                self.z[i] = -np.array(self.z[i])
        

        # Rang buoc ve dau
        for i in range(self.numberVariable - self.numberNewVariable):
            if self.signVariableConstraint[i] == -1:  # xi <= 0
                # Ham muc tieu
                self.z[i] = -self.z[i]

                # Rang buoc dang thuc
                for j in range(self.numberConstraint):
                    self.A[j][i] = -self.A[j][i]
                self.signVariableConstraint[i] = 1

            elif self.signVariableConstraint[i] == 0:  # xi tuy y => Add column xi-
                self.numberVariable += 1
                self.numberNewVariable += 1
                # Ham muc tieu
                self.z = np.append(self.z, -self.z[i])
                self.signVariableConstraint[i] = 1
                self.signVariableConstraint.append(1)
                # Rang buoc dang thuc
                self.A = np.hstack((self.A, -self.A[:, i:i+1]))

        # Rang buoc dang thuc
        for i in range(self.numberConstraint):
            if self.signInequalityConstraint[i] == 1:  # >= bi
                self.A[i] = -self.A[i]
                self.b[i] = -self.b[i]
                self.signInequalityConstraint[i] = -1  # Gan lai dau cho RB
            elif self.signInequalityConstraint[i] == 0:  # = bi => Add new constraint
                self.numberConstraint += 1
                new_constraint = np.array([-self.A[i]])
                self.A = np.vstack((self.A, new_constraint))
                self.b = np.append(self.b, -self.b[i])

                self.signInequalityConstraint[i] = -1
                self.signInequalityConstraint = np.append(self.signInequalityConstraint, -1)

    # Ham dua ve dang chuan tac
    def convert_into_canonical_form(self):
        # Ham muc tieu
        if not self.isMin:
            for i in range(self.numberVariable):
                self.z[i] = -np.array(self.z[i])
        
        # Rang buoc ve dau
        for i in range(self.numberVariable - self.numberNewVariable):
            if self.signVariableConstraint[i] == -1:  # xi <= 0
                # Ham muc tieu
                self.z[i] = -np.array(self.z[i])

                # Rang buoc dang thuc
                for j in range(self.numberConstraint):
                    self.A[j][i] = -self.A[j][i]
                self.signVariableConstraint[i] = 1

            elif self.signVariableConstraint[i] == 0:  # xi tuy y => Add column xi-
                self.numberVariable += 1
                self.numberNewVariable += 1
                # Ham muc tieu
                self.z = np.append(self.z, -self.z[i])
                self.signVariableConstraint[i] = 1
                self.signVariableConstraint.append(1)
                # Rang buoc dang thuc
                self.A = np.hstack((self.A, -self.A[:, i:i+1]))

        # Rang buoc dang thuc
        for i in range(self.numberConstraint):
            
            if self.signInequalityConstraint[i] == -1:  # <= bi
                
                self.numberVariable += 1
                self.numberNewVariable += 1
                self.signVariableConstraint[i] = 1
                self.signVariableConstraint.append(1)
                new_column = np.zeros((self.A.shape[0], 1))
                new_column[i, 0] = 1
                self.A = np.hstack((self.A, new_column))
                self.signInequalityConstraint[i] = 0  # Gan lai dau cho RB
                (self.z).append(0)

            elif self.signInequalityConstraint[i] == 1:  # >= bi 
                
                self.numberVariable += 1
                self.numberNewVariable += 1
                self.signVariableConstraint[i] = 1
                self.signVariableConstraint.append(1)
                new_column = np.zeros((self.A.shape[0], 1))
                new_column[i, 0] = -1
                self.A = np.hstack((self.A, new_column))
                self.signInequalityConstraint[i] = 0  # Gan lai dau cho RB
                (self.z).append(0)



    def copy(self):
        return copy.deepcopy(self)
    

    