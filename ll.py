t=int(input())
for z in range(t):
    n=int(input())
    print("O",end="")
    k=n-1
    i=0
    row=7
    while(i<k):
        
        print(".",end="")
        
        i+=1
        row-=1
        if(row==0):
            print()
            row=8
    k=64-n
    i=0
    while(i<k):
        print("X",end="")
        
        i+=1
        row-=1
        if(row==0):
            print()
            row=8
    print()