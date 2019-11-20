def main ():
    
    L = [1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
    dx = 0.01
    for i in (range(len(L))):
        l = L[i]
        offset = int(l/dx) / 5
        print("=====================================")
        print("Cable length %g cm" % (l))
        for i in range(5):
            print("%d" % (i*offset))
        print("=====================================")
    

if __name__ == "__main__":
    main()