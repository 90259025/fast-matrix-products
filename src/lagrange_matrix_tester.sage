import random
mat_dim = 5
function_degree = 2

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
    line_str = "["
    for i in range(mat_dim):
        for j in range(mat_dim):
            line_str += str(A_out[ind][i][j])
            if j != mat_dim-1:
                line_str += " "
            elif i != mat_dim-1:
                line_str += "; "
            else:
                line_str += "], "
    s += line_str
s += "]"
s = s.replace("], ]","]]")

printstr += s+")"
print(printstr)