

def integration(n,x,func,coor):

    func_coefs = interpol(func,coor,0,0)                    #first interpolant

    for i in range(n):                                      #integrate n times
        func_coefs = int_ana(func_coefs,coor)

    if x<coor[0]:
        return 0
    else:
        for j in range(1,len(coor)):
            if coor[j] >= x:                                     #this to find the actual value corresponding to our x
                output = polynomial(x,func_coefs[j-1])
                break

    return output