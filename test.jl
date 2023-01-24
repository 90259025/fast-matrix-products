include("functions.jl")

# test MP with a bunch of randomly generated inputs

# -t^3 + 38*t^2 + 6*t + 4
A = [4, 6, 38, -1]
# 7*t^7 + 26*t^6 + 9*t^5 + 8*t^4 + 4*t^3 + 3*t^2 + 2*t + 1)
B = [1, 2, 3, 4, 8, 9, 26, 7]
# -7*t^10 + 240*t^9 + 1021*t^8 + 518*t^7 + 458*t^6 + 233*t^5 + 168*t^4 + 109*t^3 + 62*t^2 + 14*t + 4
println(MP(A,B)==[109, 168, 233, 458])

A = [-9, -7, -7]
B = [-5, -3, 2, -8, -9, 9]
println(MP(A,B)==[38, 79, 123])

A = [3, 1, -2, 10, 4, 10, 0, 3, -6, -5, 9, -4, -3, -4, -8, 1, -7, -2, -5, -5, 0, 2, 4, -2, -1, 8, 6, 8, 8, 6, 2] 
B = [-4, 9, 7, -5, -5, 3, 10, 7, -7, -6, 5, -4, -7, 9, -2, 7, 4, 9, 5, 2, 7, -7, -1, 9, -7, 1, -4, 2, -7, -7, 9, -1, 6, 2, 8, 1, -1, 6, -8, 7, -3, 1, 9, 7, 9, 0, 3, -5, -9, -4, -5, -2, -9, 2, -1, 2, -4, 9, 4, -4, -5, 2]
println(MP(A,B)==[138, -321, -82, 10, -188, 248, -199, 168, 12, -48, 329, -64, 391, 73, 343, 290, 89, 328, -74, 141, -248, -251, -100, -399, -100, -340, -97, -151, -53, 62, 107])

A = [1, 2, 10, -9, 7, -9, -2, -9, -5, 3, 5, -4, 9, 0, 3, -7, 5, 6, -7, -4, 1, -9, 4, 10, 10, 2, -5, 1, 6, -1, 6, 7, 10, -9, 10, -9, -1, -9, -5, 1, -5, -4, 0, 9, 10, 3, 2, 10, -6, 10, 4, 7, -5, -7, -4, 5, -3, 3, -4, 7, 1, 9, 3, -6]
B = [7, -9, 10, -4, 2, 9, 2, 1, 1, 10, -1, -3, 10, -4, 9, 3, -8, 5, -4, 3, -5, 9, 9, 0, -9, 8, -7, -4, 4, -3, 2, -6, 3, -8, 8, 7, -6, 9, 9, -2, -1, 1, 1, -7, 9, -3, 2, 6, -3, -7, -6, -5, -7, -9, -8, 6, -9, -6, 9, -7, -5, 6, -4, 2, 2, 7, 7, -5, -8, 9, 2, 5, 5, -3, 8, -7, -1, -9, 8, -7, -4, 5, -5, -8, -6, -7, -2, 8, -2, -6, 5, 4, 8, -6, 0, 2, -2, -5, -7, 0, -5, 5, 7, -4, -5, 6, 10, 2, -8, -7, 8, 0, -3, 5, 0, -8, -2, 4, 1, -7, -2, 9, 7, -7, -7, -3, 7, 6]
println(MP(A,B)==[-243, 556, -149, 334, 32, 273, 200, -85, 249, -13, -334, -16, -269, -172, -622, 394, -454, -12, -18, 130, 65, -290, -136, 155, 18, 156, 601, 9, 113, -74, 125, 296, -163, 78, 4, -273, -355, -414, 38, -192, 103, 217, 54, -67, -486, 151, -381, -76, -538, 181, -391, 93, 411, 104, -117, 112, 222, -14, -83, 258, 179, -119, -78, -44])

