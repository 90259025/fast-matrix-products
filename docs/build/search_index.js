var documenterSearchIndex = {"docs":
[{"location":"","page":"-","title":"-","text":"CurrentModule = Main.OurModuleName\r\nDocTestSetup = quote\r\n    println(\"hi!\")\r\n    using Main.OurModuleName\r\nend","category":"page"},{"location":"","page":"-","title":"-","text":"Modules = [Main.OurModuleName]","category":"page"},{"location":"#Main.OurModuleName.MP-Tuple{Array{<:Union{Main.OurModuleName.IntModQ, Integer, Rational}}, Array{<:Union{Main.OurModuleName.IntModQ, Integer, Rational}}}","page":"-","title":"Main.OurModuleName.MP","text":"MP(y, z)\n\nCompute the middle product of polynomials y and z.\n\nFor vectors y = y_1y_n and z = z_1z_2n-1 of polynomial coefficients, return c_n-1 c_2n-1, where c_k = sum_i + j = k a_i b_j.\n\nExamples\n\njulia> Main.OurModuleName.MP([4, 6, 38, -1], [1, 2, 3, 4, 8, 9, 26, 7])\n4-element Vector{Int64}:\n 109\n 168\n 233\n 458\n\njulia> Main.OurModuleName.MP([-7//10, 1, -2], [9//4, -1//4, -2//5, 9//8, -4//9, -1//2])\n3-element Vector{Rational{Int64}}:\n -447//100\n  -11//16\n  161//72\n\n\n\n\n\n","category":"method"},{"location":"#Main.OurModuleName.lagrange-Union{Tuple{Vector{T1}}, Tuple{T1}} where T1<:Union{Main.OurModuleName.IntModQ, Integer, Rational}","page":"-","title":"Main.OurModuleName.lagrange","text":"lagrange(p)\n\nLagrange interpolate a vector p of polynomial outputs to a new position degree+1.\n\nFor a polynomial P, apply Lagrange interpolation to the values p = P(0) P(1)  P(d) to produce a vector of shifted values P(d+1) P(d + 2)  P(2d + 1).\n\nExample\n\njulia> Main.OurModuleName.lagrange([1, -7, -31])\n3-element Vector{Int64}:\n  -71\n -127\n -199\n\n\n\n\n\n","category":"method"},{"location":"#Main.OurModuleName.lagrange_matrix-Union{Tuple{Array{T1, 3}}, Tuple{T1}} where T1<:Union{Main.OurModuleName.IntModQ, Integer, Rational}","page":"-","title":"Main.OurModuleName.lagrange_matrix","text":"lagrange_matrix(A)\n\nLagrange interpolate all entries of a list of d+1 matrices and return 3(d+1) matrices\n\nApply the following to each entry of the matrices, with a = d+1 2d+1 3d+1: For a polynomial P, apply Lagrange interpolation to the values p = P(0) P(1)  P(d) to produce a vector of shifted values P(a) P(a + 1)  P(a + d).\n\nExample\n\njulia> Main.OurModuleName.lagrange_matrix([0:2 3:5 6:8 ;;; 9:11 12:14 15:17])\n3×3×6 Array{Int64, 3}:\n[:, :, 1] =\n 18  21  24\n 19  22  25\n 20  23  26\n\n[:, :, 2] =\n 27  30  33\n 28  31  34\n 29  32  35\n\n[:, :, 3] =\n 36  39  42\n 37  40  43\n 38  41  44\n\n[:, :, 4] =\n 45  48  51\n 46  49  52\n 47  50  53\n\n[:, :, 5] =\n 54  57  60\n 55  58  61\n 56  59  62\n\n[:, :, 6] =\n 63  66  69\n 64  67  70\n 65  68  71\n\n\n\n\n\n","category":"method"},{"location":"#Main.OurModuleName.matrix_product-Union{Tuple{S}, Tuple{T}, Tuple{Array{T, 3}, S}} where {T<:Union{Main.OurModuleName.IntModQ, Integer, Rational}, S<:Integer}","page":"-","title":"Main.OurModuleName.matrix_product","text":"matrix_product(T, a)\n\nCompute the product T(a-1)T(a-2)...T(a) a matrix T.\n\nExample\n\njulia> Main.OurModuleName.matrix_product(BigInt.([0:2 3:5 6:8 ;;; 9:11 12:14 15:17 ]), 2^2)\n3×3 Matrix{BigInt}:\n  688662   762210   835758\n  907092  1003968  1100844\n 1125522  1245726  1365930\n\n\n\n\n\n","category":"method"}]
}
