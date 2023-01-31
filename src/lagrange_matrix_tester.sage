import random
mat_dim = 64
function_degree = 10

USE_BIG_INT = True

A_in = [_ for _ in range(function_degree+1)]
A_out =  [_ for _ in range(3*(function_degree+1))]
for _ in range(function_degree+1):
    A_in[_] = list(matrix.identity(mat_dim))
for _ in range(3*(function_degree+1)):
    A_out[_] = list(matrix.identity(mat_dim))

R.<x> = ZZ[]

for i in range(mat_dim):
    for j in range(mat_dim):
        f = sum([random.randint(-10,9)*x^i for i in range(function_degree+1)])
        for ind in range(function_degree+1):
            A_in[ind][i][j] = f(ind)
        for ind in range(3*(function_degree+1)):
            A_out[ind][i][j] = f(ind+function_degree+1)
A_in = [matrix(_) for _ in A_in]
A_out = [matrix(_) for _ in A_out]

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
printstr += "println(lagrange_matrix(" + s + ") == "

s = "["
for ind in range(3*(function_degree+1)):
    line_str = ""
    for i in range(mat_dim):
        for j in range(mat_dim):
            line_str += str(A_out[ind][i][j])
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

if not USE_BIG_INT:
    print(printstr)

else:
# optional BigInt conversion

# if you're doing this in a sage notebook this will probably crash for large matrices

    printstr = printstr.replace('([','(BigInt.([').replace('== [',') == BigInt.([')
    printstr += ')'
    print(printstr)

    # with open('text_dump.txt','w') as o:
        # o.write(printstr)
