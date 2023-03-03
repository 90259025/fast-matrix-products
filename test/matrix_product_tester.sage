import random
mat_dim = 64
function_degree = 4
len_prod = 2^12

# optional argument to reduce answer mod p
# p = 0
p = Primes().next(2^11)

USE_BIG_INT = True

A_in = [_ for _ in range(function_degree+1)]
A_out =  [[0 for _ in range(mat_dim)] for _ in range(mat_dim)]
for _ in range(function_degree+1):
    A_in[_] = list(matrix.identity(mat_dim))


R.<x> = ZZ[]

for i in range(mat_dim):
    for j in range(mat_dim):
        f = sum([random.randint(-10,9)*x^i for i in range(function_degree+1)])
        for ind in range(function_degree+1):
            A_in[ind][i][j] = f(ind)
        A_out[i][j] = f(x)
A_in = [matrix(_) for _ in A_in]
A_out = matrix(A_out)

if p > 0:
    A_out_step = A_out(0) % p
    for i in range(1,len_prod):
        A_out_step = (A_out_step*A_out(i))%p
        
    A_out = A_out_step % p
else:
    A_out = prod([A_out(i) for i in range(len_prod)])

printstr =""
s = "["
for ind in range((function_degree+1)):
    line_str = "["
    for i in range(mat_dim):
        for j in range(mat_dim):
            line_str += str(A_in[ind][i][j])
            if j != mat_dim-1:
                line_str += " "
            elif i != mat_dim-1:
                line_str += "; "
            else:
                line_str += "], "
    s += line_str
s += "]"
s = s.replace("], ]","]]")
if p > 0:
    printstr += "println(matrix_product(" + s + " " + ") == "
else:
    printstr += "println(matrix_product(" + s + " " + ") == "

s = "["
line_str = ""
for i in range(mat_dim):
    for j in range(mat_dim):
        line_str += str(A_out[i][j])
        if j != mat_dim-1:
            line_str += " "
        elif i != mat_dim-1:
            line_str += "; "
        else:
            line_str += " ;;; "
s += line_str
s += "]"
s = s.replace("], ]","]")

printstr += s+")"
printstr = printstr.replace('[[','[').replace(']])','])').replace('], [',' ;;; ').replace(']]','').replace('],','').replace(' ;;; ]',']')

# if p > 0:
#     printstr = printstr.replace(') == ', "]), " + str(len_prod) + '),' + str(p) + ') == ')
# else:
printstr = printstr.replace(') == ', "]), " + str(len_prod) + ') == ')

if not USE_BIG_INT:
    print(printstr)

else:
# optional BigInt conversion

# if you're doing this in a sage notebook this will probably crash for large matrices
    
    if p > 0:
        printstr = printstr.replace('([','(IntModQ.([').replace('== [',' == IntModQ.([')
    else:
        printstr = printstr.replace('([','(BigInt.([').replace('== [',' == BigInt.([')
    printstr += ')'
    print(printstr)

    # with open('text_dump.txt','w') as o:
        # o.write(printstr)
