from linear_programming_class import *
import geometric_2constraint as g2cmt
import geometry_method as gmt
import simplex_method as smt

def main():
    print("\n>>>>>>>>>>>>>>>>>>>> BAI TOAN QUY HOACH TUYEN TINH <<<<<<<<<<<<<<<<<<<<<<<<")
    check_run = False
    while True:
        if check_run == False:
            P = Problem()
            P.input_problem()
            P.print_problem()
            P_copy = P.copy()
            sum_constraint = np.sum(P.A, axis=0) / np.sum(P.b)
            result =True
            if all(x >= 0 for x in sum_constraint):
                result = True
            else:
                result = False
        else:
            cnt = input("Ban co muon tiep tuc thuc hien bai toan nay khong (co / khong)? ")
            if cnt.lower() == "khong":
                print("Moi ban nhap bai toan moi: ")
                P = Problem()
                P.input_problem()
                P.print_problem()
            else:
                P = P_copy.copy()
            
        # Menu program
        if P.numberVariable == 2: # Trường hợp 2 biến
            print("Phuong phap cho bai toan nay la: ")
            print("1. Phuong phap hinh hoc.")
            print("2. Phuong phap don hinh.")

            chonbaitap = int(input("Moi ban chon, vui long nhap so tuong ung voi phuong phap ma ban chon: "))
            
        elif P.numberConstraint== 2 and result == True: # nhiều hơn 2 biến có 2 ràng buộc
            print("Phuong phap cho bai toan nay la: ")
            print("1.  Phuong phap hinh hoc voi 2 rang buoc.")
            print("2.  Phuong phap don hinh.")
            chonbaitap = int(input("Moi ban chon, vui long nhap so tuong ung voi phuong phap ma ban chon: ")) 
            if chonbaitap == 1:
                chonbaitap = 3
        else:
            print("Phuong phap cho bai toan nay la phuong phap don hinh.")
            chonbaitap = 2 #defaut là đơn hình 

        if chonbaitap == 1:
            print("Phuong phap hinh hoc:")
            # ham phuong phap hinh hoc
            gmt.geometric_method(P)
            tieptuc = input("Ban co muon tiep tuc khong (co / khong)? ")
            if tieptuc.lower() == "khong":
                print("Ket thuc chuong trinh")
                break
            check_run = True
        elif chonbaitap == 2:
            print("Phuong phap don hinh: ")
            # ham phuong phap don hinh
            smt.solve_problem(P)
            tieptuc = input("Ban co muon tiep tuc khong (co / khong)? ")
            if tieptuc.lower() == "khong":
                print("Ket thuc chuong trinh")
                break
            check_run = True
        elif chonbaitap == 3:
            print("Phuong phap hinh hoc voi 2 rang buoc:")
            # Hinh hoc 2 rang buoc loi
            g2cmt.geometric_2constraint(P)
            tieptuc = input("Ban co muon tiep tuc khong (co / khong)? ")
            if tieptuc.lower() == "khong":
                print("Ket thuc chuong trinh")
                break
            check_run = True
        else:
            print("Ban da chon sai moi ban chon lai!")


if __name__ == "__main__":
    main()
