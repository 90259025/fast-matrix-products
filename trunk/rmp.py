# from https://hal.inria.fr/inria-00072675/document

# IN: [[a0,...,a_{2n-1}], [x0,...,x_{n-1}]
# OUT: [c_{n-1},...,c_{2n-2}] where c_k = sum_{i+j=k}a_i x_j
# If n = 2^k this takes exactly 3^k ring operations

# Assume n is a power of two for now
def RMP(a,x):
    n = len(a)//2
    if n == 1:
        return [a[0] * x[0]]
    alpha = RMP([a[i] + a[n//2 + i] for i in range(n)], x[n//2:])
    beta = RMP(a[n//2:3*n//2],[x[i] - x[i+n//2] for i in range(n//2)])
    gamma = RMP([a[i] + a[i+n//2] for i in range(n//2,3*n//2)], x[:n//2])
    h = [alpha[i] + beta[i] for i in range(n//2)]
    l = [gamma[i] - beta[i] for i in range(n//2)]
    return h + l #concatenated as lists
    
# can be tested with
# var('t')
# a = [1,2,3,4,8,9,26,7]
# x = [4,6,38,-1]
# print(expand(sum(a[i] * t^i for i in range(len(a)))*sum(x[i] * t^i for i in range(len(x)))))
# RMP(a,x)