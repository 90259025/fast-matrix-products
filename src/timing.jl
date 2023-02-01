# small testing program to compare MP against DSP

using DSP
include("functions.jl")

A = [-6 -4 -8 -9  7]#  0  7  9  1 12]
B = [-7 -3 -9  0 12  9 10 16  5 10]

@time DSP.conv(A,B)
@time MP(A,B)

C = [5 1 -7 -8 3 7 5 -2 -4 8 1 9 3 -10 -17 0 -22 -2 31 11 65 34 -13 -103 3 -112 28 80 25 191 103 24 -293 38 -326 122 161 37 417 228 131 -629 129 -718 316 280 41 773 427 338 -1153 300 -1342 646 443 31 1289 718 675 -1907 575 -2252 1148 656 1 1995 1119 1172 -2933 978 -3502 1858 925 -55 2921 1648 1859 -4273 1533 -5146 2812 1256 -143 4097 2323 2766 -5969 2264 -7238 4046 1655 -269 5553 3162 3923 -8063 3195 -9832 5596 2128 -439 7319 4183 5360 -10597 4350 -12982 7498 2681 -659 9425 5404 7107 -13613 5753 -16742 9788 3320 -935 11901 6843 9194 -17153 7428 -21166 12502 4051 -1273 14777 8518 11651 -21259 9399 -26308 15676 4880 -1679 18083 10447 14508 -25973 11690 -32222 19346]
D = C.+1
C = C[1:length(C)÷2]

println("using conv, length", length(C))

@time DSP.conv(C,D)

println("using MP, length", length(C))
@time MP(C,D)

for i = 1:10
	global C = [C C]
	global D = [D D]
	println("using conv, length", length(C))

	@time DSP.conv(C,D)

	println("using MP, length", length(C))
	@time MP(C,D)
end

# results from last run were
#   2.124152 seconds (4.47 M allocations: 226.390 MiB, 14.01% gc time, 99.98% compilation time)
#   0.870826 seconds (1.32 M allocations: 62.266 MiB, 15.68% gc time, 99.98% compilation time)
# using conv, length72
#   0.622055 seconds (1.73 M allocations: 87.912 MiB, 3.69% gc time, 98.55% compilation time)
# using MP, length72
#   0.274561 seconds (162.55 k allocations: 6.605 MiB, 99.63% compilation time)
# using conv, length144
#   0.017035 seconds (89 allocations: 343.828 KiB)
# using MP, length144
#   0.002264 seconds (37.13 k allocations: 2.599 MiB)
# using conv, length288
#   0.031882 seconds (89 allocations: 701.047 KiB)
# using MP, length288
#   0.017939 seconds (111.41 k allocations: 7.829 MiB, 60.52% gc time)
# using conv, length576
#   0.067901 seconds (93 allocations: 1.364 MiB)
# using MP, length576
#   0.025036 seconds (334.24 k allocations: 23.548 MiB)
# using conv, length1152
#   0.142963 seconds (94 allocations: 2.867 MiB)
# using MP, length1152
#   0.070843 seconds (1.00 M allocations: 70.766 MiB, 20.55% gc time)
# using conv, length2304
#   0.375174 seconds (95 allocations: 5.729 MiB)
# using MP, length2304
#   0.171537 seconds (3.01 M allocations: 212.538 MiB, 11.44% gc time)
# using conv, length4608
#   0.760401 seconds (95 allocations: 11.454 MiB)
# using MP, length4608
#   0.707867 seconds (9.02 M allocations: 638.088 MiB, 12.86% gc time)
# using conv, length9216
#   1.586631 seconds (95 allocations: 22.903 MiB)
# using MP, length9216
#   2.048973 seconds (27.07 M allocations: 1.870 GiB, 11.61% gc time)
# using conv, length18432
#   3.339414 seconds (95 allocations: 45.801 MiB, 0.12% gc time)
# using MP, length18432
#   5.978904 seconds (81.22 M allocations: 5.613 GiB, 18.06% gc time)
# using conv, length36864
#   8.827531 seconds (95 allocations: 91.598 MiB, 0.08% gc time)
# using MP, length36864
#  14.234480 seconds (243.67 M allocations: 16.842 GiB, 11.75% gc time)
# using conv, length73728
#  22.802135 seconds (95 allocations: 183.192 MiB, 0.03% gc time)
# using MP, length73728
#  51.077333 seconds (731.00 M allocations: 50.534 GiB, 9.78% gc time)