include("../src/functions.jl")

using BenchmarkTools


test_MP = false 
test_lagrange = false
test_matrix_product = false
test_matrix_product_slow = false
test_matrix_product_mod_q = false
test_matrix_product_mod_q_really_slow = false
test_matrix_product_really_slow = false

# test_MP = true 
# test_lagrange = true
# test_matrix_product = true
# test_matrix_product_slow = true
test_matrix_product_mod_q = true
test_matrix_product_mod_q_really_slow = true
# test_matrix_product_really_slow = true




if test_MP
    println("Testing middle product...")
    # test MP with a bunch of randomly generated inputs

    # -t^3 + 38*t^2 + 6*t + 4
    A = [4, 6, 38, -1]
    # 7*t^7 + 26*t^6 + 9*t^5 + 8*t^4 + 4*t^3 + 3*t^2 + 2*t + 1)
    B = [1, 2, 3, 4, 8, 9, 26, 7]
    # -7*t^10 + 240*t^9 + 1021*t^8 + 518*t^7 + 458*t^6 + 233*t^5 + 168*t^4 + 109*t^3 + 62*t^2 + 14*t + 4
    println(MP(A, B) == [109, 168, 233, 458])
    println(MP(IntModQ.(A), IntModQ.(B)) == IntModQ.(MP(A, B)))

    A = [-9, -7, -7]
    B = [-5, -3, 2, -8, -9, 9]
    println(MP(A,B)==[38, 79, 123])
    println(MP(IntModQ.(A), IntModQ.(B)) == IntModQ.(MP(A, B)))

    A = [3, 1, -2, 10, 4, 10, 0, 3, -6, -5, 9, -4, -3, -4, -8, 1, -7, -2, -5, -5, 0, 2, 4, -2, -1, 8, 6, 8, 8, 6, 2] 
    B = [-4, 9, 7, -5, -5, 3, 10, 7, -7, -6, 5, -4, -7, 9, -2, 7, 4, 9, 5, 2, 7, -7, -1, 9, -7, 1, -4, 2, -7, -7, 9, -1, 6, 2, 8, 1, -1, 6, -8, 7, -3, 1, 9, 7, 9, 0, 3, -5, -9, -4, -5, -2, -9, 2, -1, 2, -4, 9, 4, -4, -5, 2]
    println(MP(A,B)==[138, -321, -82, 10, -188, 248, -199, 168, 12, -48, 329, -64, 391, 73, 343, 290, 89, 328, -74, 141, -248, -251, -100, -399, -100, -340, -97, -151, -53, 62, 107])
    println(MP(IntModQ.(A), IntModQ.(B)) == IntModQ.(MP(A, B)))

    A = [1, 2, 10, -9, 7, -9, -2, -9, -5, 3, 5, -4, 9, 0, 3, -7, 5, 6, -7, -4, 1, -9, 4, 10, 10, 2, -5, 1, 6, -1, 6, 7, 10, -9, 10, -9, -1, -9, -5, 1, -5, -4, 0, 9, 10, 3, 2, 10, -6, 10, 4, 7, -5, -7, -4, 5, -3, 3, -4, 7, 1, 9, 3, -6]
    B = [7, -9, 10, -4, 2, 9, 2, 1, 1, 10, -1, -3, 10, -4, 9, 3, -8, 5, -4, 3, -5, 9, 9, 0, -9, 8, -7, -4, 4, -3, 2, -6, 3, -8, 8, 7, -6, 9, 9, -2, -1, 1, 1, -7, 9, -3, 2, 6, -3, -7, -6, -5, -7, -9, -8, 6, -9, -6, 9, -7, -5, 6, -4, 2, 2, 7, 7, -5, -8, 9, 2, 5, 5, -3, 8, -7, -1, -9, 8, -7, -4, 5, -5, -8, -6, -7, -2, 8, -2, -6, 5, 4, 8, -6, 0, 2, -2, -5, -7, 0, -5, 5, 7, -4, -5, 6, 10, 2, -8, -7, 8, 0, -3, 5, 0, -8, -2, 4, 1, -7, -2, 9, 7, -7, -7, -3, 7, 6]
    println(MP(A,B)==[-243, 556, -149, 334, 32, 273, 200, -85, 249, -13, -334, -16, -269, -172, -622, 394, -454, -12, -18, 130, 65, -290, -136, 155, 18, 156, 601, 9, 113, -74, 125, 296, -163, 78, 4, -273, -355, -414, 38, -192, 103, 217, 54, -67, -486, 151, -381, -76, -538, 181, -391, 93, 411, 104, -117, 112, 222, -14, -83, 258, 179, -119, -78, -44])
    println(MP(IntModQ.(A), IntModQ.(B)) == IntModQ.(MP(A, B)))

    A = [-7//10, 1, -2] 
    B = [9//4, -1//4, -2//5, 9//8, -4//9, -1//2]
    println(MP(A,B) == [-447//100, -11//16, 161//72])
end



if test_lagrange
    println("Testing Lagrange interpolation...")
    # Test lagrange
    #println(lagrange([1//1,2,3])==[9//1,10//1,11//1])
    println(lagrange(IntModQ.([1,2,3])) == IntModQ.(lagrange([1,2,3])))
    #println(lagrange([0//1, 127, 2186, 16383, 78124, 279935, 823542, 2097151])==[19487170//1,  35831807//1,  62748516//1,  105413503//1,  170859374//1,  268435455//1,  410338672//1,  612220031//1])
    # maybe this should test if it equals ints instead of rationals
    println(lagrange([1,8,27,64])==[125//1, 216//1, 343//1, 512//1])
    println(lagrange(IntModQ.([1,8,27,64])) == IntModQ.(lagrange([1,8,27,64])))

    println(lagrange([1//4, -7//4, -31//4])==[-71//4, -127//4, -199//4])
    println(lagrange([-7//2, -797//360, -179947//90, -8894831//120, -85280963//90, -493110857//72, -1036653841//30, -48973005053//360, -40185312907//90, -51009804373//40]) == [-58738621183//18, -2750188672237//360, -498472099889//30, -12230002801613//360, -5930345300803//90, -2931007133459//24, -19579854846683//90, -134741703303133//360, -6243471811443//10, -364748503947917//360])

end

if test_matrix_product
    println("Testing matrix product")
    println(matrix_product(BigInt.([0:2 3:5 6:8 ;;; 9:11 12:14 15:17 ]), 2)==BigInt[96 123 150; 126 162 198; 156 201 246])
    println(matrix_product(BigInt.([0:2 3:5 6:8 ;;; 9:11 12:14 15:17 ]), 2^6)==BigInt[ 67826000765456369170712628352974773878226653389534849285491327727605255841356333691411401979742066166545865337319507611473143691259091192338466020753449038154654856809893066702848 68184234440547279012250523801548761166222171901823712345710233863673651945626467606514319300106141491497929396470561568507817460835360022057843538621599050528636869817079305338880 68542468115638188853788419250122748454217690414112575405929139999742048049896601521617236620470216816449993455621615525542491230411628851777221056489749062902618882824265543974912; 89339072384741737275296403532737925829442454479634481730185956085392065149689416356410906359023268980777590359650940891214509295389330845795403195495270435763317387315963953676288 89810930726209786025931948251140819441623361894827017093661607554086443454793826242216766826808322006452693139910882059907878364450680297343025545425559643199399956651535742533632 90282789067677834776567492969543713053804269310019552457137259022780821759898236128022627294593375032127795920170823228601247433512029748890647895355848850635482525987107531390976;110852144004027105379880178712501077780658255569734114174880584443178874458022499021410410738304471795009315381982374170955874899519570499252340370237091833371979917822034840649728 111437627011872293039613372700732877717024551887830321841612981244499234963961184877919214353510502521407456883351202551307939268066000572628207552229520235870163043485992179728384 112023110019717480699346566688964677653390848205926529508345378045819595469899870734428017968716533247805598384720030931660003636612430646004074734221948638368346169149949518807040])
    println(mod.(matrix_product(BigInt.([0:2 3:5 6:8 ;;; 9:11 12:14 15:17 ]), 2^8),389)==BigInt[385 384 383;368 106 233;351 217 83])
    println(matrix_product(BigInt.([9 3 -4 0 -5; 0 -6 -6 -6 5; 9 0 3 6 -6; 8 -9 9 -3 -6; -6 7 -6 -9 -2 ;;; 7 1 2 -4 -9; -5 -12 -13 1 6; 5 9 -7 -3 -15; 2 -10 11 -2 -8; -6 3 -14 -16 -12 ]), 8)  == BigInt.([6310130646821 -66269513953 10920862184868 14273655142439 20131255320171; 3349865639168 1000465803395 5369195287318 6469662473094 9624086613414; 1096285155135 -327979329759 2138306477583 2608189330287 3606051420576; 2380346656383 1019582406835 3434682629988 4762133848821 7031906325642; 9406955627743 310732742814 15757365394583 21520161646192 30290168379768]))
    println(mod.(matrix_product(BigInt.([3 -7 7 4 -10; 4 -6 3 9 9; -2 0 1 6 2; 1 -7 5 4 1; 6 8 -5 -4 3 ;;; 12 2 14 3 -17; 11 -5 11 6 11; -5 -7 10 0 2; 3 -12 14 7 -8; 4 -2 -7 2 -7 ]), 8),389)  == BigInt.([114 23 205 54 2; 23 103 307 193 336; 17 134 238 54 301; 180 183 375 245 30; 315 181 202 192 358]))    
    println(mod.(matrix_product(BigInt.([9 9 -4; 6 1 1; 2 2 -8 ;;; 0 4 -20; -6 13 6; -1 -12 9 ;;; -25 3 -48; -22 43 11; 10 -40 42 ]), 2),2053)  == BigInt.([2003 201 1891; 2046 25 1948; 2049 130 1953]))
    println(matrix_product(BigInt.([-7 -8 -7; 0 6 3; -8 -6 -1 ;;; -13 -9 -16; -3 7 19; -6 -10 7 ;;; -27 -10 -23; 4 16 53; -12 -4 15 ]), 2)  == BigInt.([157 77 -89; -36 12 135; 128 40 7]))
    println(matrix_product(BigInt.([1 -3 0; 1 -1 8; -1 -2 4 ;;; 3 -12 9; 8 1 10; -4 -6 -7 ;;; 9 -31 24; 27 -3 14; -17 -22 -38 ]), 16)  == BigInt.([26232097201847637501656274121044930537 -252761813710026155381080448027519775048 -252158182194751717460963166378577530153; 87931715687714391870923057642403710511 -991443661719657974170019250074987959386 -1066541778933255595259716322568147718527; 49033727672125983633121590997052304384 -1118752398295368069898807765066599654123 -1223672762795933290451505791107355310744]))
    println(mod.(matrix_product(BigInt.([4 4 -5 4 1; -9 -10 9 -2 0; -3 -4 -5 5 -3; -9 4 -7 4 -1; -9 -8 0 -9 1 ;;; 3 13 1 1 -8; -8 0 13 -4 7; -22 -1 -17 5 7; -16 18 -22 1 -10; 0 -14 -13 -9 2 ;;; -18 124 31 -156 15; -11 140 39 26 106; -93 126 5 -81 127; -73 150 -79 34 -71; 97 -122 -118 33 -141 ;;; -185 541 103 -833 232; -78 734 99 136 519; -240 635 157 -433 591; -366 676 -238 235 -220; 438 -560 -531 171 -764 ;;; -720 1588 235 -2612 949; -341 2322 205 398 1612; -463 1928 583 -1327 1753; -1225 2064 -583 856 -493; 1251 -1700 -1612 483 -2395 ]), 256),2053)  == BigInt.([30 1832 46 1524 704; 1950 975 1437 1939 1120; 1116 1999 1033 811 773; 1635 1107 642 63 685; 961 525 1309 1212 2006]))
    println(mod.(matrix_product(BigInt.([9 9 -4; 6 1 1; 2 2 -8 ;;; 0 4 -20; -6 13 6; -1 -12 9 ;;; -25 3 -48; -22 43 11; 10 -40 42 ]), 2),2053)  == BigInt.([2003 201 1891; 2046 25 1948; 2049 130 1953]))

    println(matrix_product(BigInt.([-7 -8 -7; 0 6 3; -8 -6 -1 ;;; -13 -9 -16; -3 7 19; -6 -10 7 ;;; -27 -10 -23; 4 16 53; -12 -4 15 ]), 2)  == BigInt.([157 77 -89; -36 12 135; 128 40 7]))
    println(matrix_product(BigInt.([1 -3 0; 1 -1 8; -1 -2 4 ;;; 3 -12 9; 8 1 10; -4 -6 -7 ;;; 9 -31 24; 27 -3 14; -17 -22 -38 ]), 16)  == BigInt.([26232097201847637501656274121044930537 -252761813710026155381080448027519775048 -252158182194751717460963166378577530153; 87931715687714391870923057642403710511 -991443661719657974170019250074987959386 -1066541778933255595259716322568147718527; 49033727672125983633121590997052304384 -1118752398295368069898807765066599654123 -1223672762795933290451505791107355310744]))

    # these work :(
    println(matrix_product(BigInt.([9 -10 2 -10 -4; 5 6 4 6 -8; 5 3 -1 8 7; 7 -3 -2 -10 -1; 4 1 -4 8 8 ;;; 12 -1 3 -23 2; 12 16 7 7 -3; 10 -5 9 22 -5; 16 -12 6 -11 -4; 8 0 -1 1 2 ;;; 23 18 -4 -46 26; 21 40 28 6 -2; 7 -19 35 50 -37; 35 -29 16 -6 -15; 26 1 2 -14 -12 ]), 8)  == BigInt.([2514795400589982 6276325888390667 -5005956512331281 -8660388752240875 6490100197051590; -1397398190457354 -7701891907328895 -163317580696595 5837588421589353 -3770682992743052; -5038168102813528 -13205158126189960 12728539747107478 19319952001675195 -15851797240870838; -768540367870240 -7278170171738767 -3014035567670593 3490621444281859 -1469620985139727; -4794719087470696 -12453255368989930 12816778060172705 18719497506760717 -15642153421738003]))
    println(matrix_product(BigInt.([6 -1 4 5 0; 6 6 -6 9 -2; 2 8 6 6 -10; -5 -4 -10 6 3; 1 -1 3 0 -5 ;;; 0 -4 21 -3 -10; 12 4 -7 8 -2; 6 4 13 24 -11; 3 1 -24 11 1; 6 -8 1 6 -18 ;;; -68 7 74 -43 -28; 38 -44 -18 -45 12; 46 -42 40 104 -32; 33 40 -60 54 -1; -7 -37 -49 52 -23 ;;; -258 32 193 -163 -54; 120 -192 -63 -198 46; 158 -166 117 294 -103; 97 167 -130 159 3; -62 -112 -195 186 4 ]), 8)  == BigInt.([8727017520137673444 -3974849912531237820 -341925176866990514 455623479198359700 -1149755949346523661; 88455489595158670404 -79163411568133604172 -39834637151588873038 9168257115014588682 4581244973676471531; -58180942634656289304 37204280785188420680 9789339334064091724 -32440291570507708868 7249896573211784836; 80550283783751383620 -78868428246146014396 -48707192722843682176 2892054471606742561 10553221402804284184; -80950542361962043968 64845054942616640928 27766642523684449612 -18465554548822424424 849381354778630617])) 

    println(mod.(matrix_product(BigInt.([4 4 -5 4 1; -9 -10 9 -2 0; -3 -4 -5 5 -3; -9 4 -7 4 -1; -9 -8 0 -9 1 ;;; 3 13 1 1 -8; -8 0 13 -4 7; -22 -1 -17 5 7; -16 18 -22 1 -10; 0 -14 -13 -9 2 ;;; -18 124 31 -156 15; -11 140 39 26 106; -93 126 5 -81 127; -73 150 -79 34 -71; 97 -122 -118 33 -141 ;;; -185 541 103 -833 232; -78 734 99 136 519; -240 635 157 -433 591; -366 676 -238 235 -220; 438 -560 -531 171 -764 ;;; -720 1588 235 -2612 949; -341 2322 205 398 1612; -463 1928 583 -1327 1753; -1225 2064 -583 856 -493; 1251 -1700 -1612 483 -2395 ]), 256),2053)  == BigInt.([30 1832 46 1524 704; 1950 975 1437 1939 1120; 1116 1999 1033 811 773; 1635 1107 642 63 685; 961 525 1309 1212 2006]))

    
end

if test_matrix_product_slow
    # println("slow...")
    # print("linear, 5 by 5, product length 2^16: ")
    # @time println(mod.(matrix_product(BigInt.([-6 7 6 1 4; 7 0 7 6 8; 8 -10 -10 -4 -2; -1 -2 -2 -2 -6; 8 6 5 6 6 ;;; -1 7 3 5 -2; 2 -6 5 -1 7; 2 -16 -10 1 -6; 3 4 -11 2 3; 7 2 2 7 0 ]), 65536),65537)  == BigInt.([3838 34641 30854 14399 56148; 21866 7865 55611 59958 5428; 36113 28614 13942 5160 6730; 21496 43923 42481 11067 28051; 719 51431 41049 44630 64748]))
    characteristic_q = 2053
    print("linear, 5 by 5, prduct length 2^21: ")
    @time println(matrix_product(IntModQ.([7 -6 9 7 -10; -4 -4 7 -10 6; -10 1 -9 3 7; -3 -3 -3 2 -8; 7 -2 -4 -3 4 ;;; 11 -5 13 -2 -8; -5 -1 8 -13 11; -11 4 0 10 8; -13 6 -2 9 -15; -2 -2 -14 4 0 ]), 2097152)  == IntModQ.([1612 406 442 54 1983; 1440 1618 1502 731 1199; 218 1429 1137 983 973; 991 550 228 922 79; 1186 327 1014 200 2]))
end

if test_matrix_product_mod_q
    println("Testing mod q")
    characteristic_q = 2053
    # this one needs q = 389
    # println(matrix_product(IntModQ.([9 9 -4; 6 1 1; 2 2 -8 ;;; 0 4 -20; -6 13 6; -1 -12 9 ;;; -25 3 -48; -22 43 11; 10 -40 42 ]), 2)  == IntModQ.(matrix_product([9 9 -4; 6 1 1; 2 2 -8 ;;; 0 4 -20; -6 13 6; -1 -12 9 ;;; -25 3 -48; -22 43 11; 10 -40 42 ], 2)))

    println(matrix_product([IntModQ.([-6 -10 -7 0 -7; -4 -2 -9 2 -5; 1 8 7 6 6; -3 6 -7 -2 -2; 2 1 -7 -7 7]), IntModQ.([-15 -17 -4 -10 2; 0 -3 -16 3 -6; 7 -1 8 10 6; 5 -3 -17 -6 6; 2 -5 -17 -14 12 ])], 16)  == IntModQ.([639 371 881 237 1756; 218 248 705 1397 272; 702 207 1112 1129 1246; 132 1548 164 48 1383; 1974 843 690 1231 793]))
    println(matrix_product([IntModQ.([7 -2 -4; -10 9 -2; 5 -6 9]), IntModQ.([11 -8 -7; -20 9 0; 5 -10 8 ])], 1048576)  == IntModQ.([1671 373 958; 1122 1996 73; 1950 325 461]))
    @time println(matrix_product([IntModQ.([-10 4 -8; 8 -6 1; -6 -1 0]), IntModQ.([-8 0 -6; 15 -15 -9; -15 0 -1])], 16777216)  == IntModQ.([415 2032 607; 579 298 1038; 1391 964 1739]))
    println(matrix_product([IntModQ.([-3 5 -7; 4 -7 6; -2 -3 -10]), IntModQ.([10 3 -10; 1 4 15; 8 -5 -16]), IntModQ.([31 -5 -9; -22 29 28; 22 -7 -14])], 524288)  == IntModQ.([1270 1302 1507; 203 1882 185; 620 210 1459]))
    @time println(matrix_product([IntModQ.([8 2 -4; -9 -9 -4; -1 -3 -4]), IntModQ.([7 7 -4; -13 -20 5; -28 -35 -7]), IntModQ.([-14 -18 -102; 1 -81 86; -235 -213 -62]), IntModQ.([-103 -187 -568; 39 -390 455; -988 -795 -391]), IntModQ.([-332 -710 -1864; 83 -1313 1496; -2869 -2183 -1384])], 524288)  == IntModQ.([208 1496 2050; 1591 817 51; 1908 1956 860]))
    @time println(matrix_product([IntModQ.([8 0 -5 1 -3 2 4 -6 2 -3; -7 3 8 -5 -7 -4 -3 4 6 -5; -5 2 -8 0 0 1 9 -4 5 0; -3 -10 -4 -4 -9 -9 6 5 -7 0; -4 -10 0 1 -1 -7 8 0 -5 6; -8 4 9 -8 3 -1 -6 7 7 -10; -9 1 -7 -3 -9 -1 6 -2 -7 -2; 9 6 1 6 6 9 2 -9 -9 -7; -3 -8 -4 -3 -9 0 0 0 0 9; -1 -10 8 5 4 2 -6 -8 8 2]), IntModQ.([12 -12 0 -2 -9 6 16 -21 24 -16; -17 6 10 -14 -10 -14 -7 -11 0 -29; -7 6 -3 3 11 21 27 -1 -4 -3; -1 -29 8 -3 -23 -17 0 -25 11 -13; -5 1 -1 15 -6 3 24 -4 17 -3; -11 10 -1 -27 -7 14 -20 1 -1 -2; -18 18 8 -8 -10 10 4 10 8 -20; 15 -5 9 -2 8 4 6 -17 -15 -6; 4 -14 -17 1 6 6 -15 -21 -5 19; 14 -10 -8 -6 -8 16 -30 -4 19 1]), IntModQ.([-30 -60 25 -9 -33 78 198 -60 94 -51; -17 151 84 -113 -89 -86 -19 -54 14 -217; 69 116 92 86 148 131 43 -58 -7 40; 47 -154 54 -80 -135 69 12 -223 121 -60; 28 26 -78 119 -23 75 208 -118 199 -138; -2 72 17 -90 -141 181 -160 -7 -97 0; -65 95 109 7 -79 77 -26 -6 123 -124; 101 -146 123 106 52 -87 -64 -99 -59 -93; 95 -54 -146 -5 157 94 -138 -162 -38 155; 161 -38 -144 -119 -90 124 -172 52 174 -14]), IntModQ.([-322 -198 160 4 -27 374 886 -39 242 -78; 113 774 422 -554 -424 -364 -33 -155 84 -935; 511 596 571 387 699 535 -45 -373 62 177; 249 -631 146 -457 -561 495 84 -919 497 -159; 155 41 -435 487 -22 347 896 -666 859 -747; 73 310 183 -221 -741 854 -720 -23 -419 -70; -222 316 458 24 -438 320 -60 -158 548 -404; 489 -765 577 660 234 -498 -496 -435 -195 -472; 552 -230 -631 -57 738 468 -567 -591 -195 717; 758 -136 -616 -598 -374 500 -624 316 743 -175]), IntModQ.([-1236 -504 567 85 129 1146 2608 222 498 -19; 613 2403 1360 -1781 -1315 -1064 -43 -368 246 -2789; 1823 1878 1944 1116 2144 1581 -435 -1312 317 456; 761 -1874 272 -1500 -1661 1651 282 -2635 1409 -304; 460 -26 -1420 1389 75 1053 2616 -2188 2507 -2418; 316 892 713 -444 -2389 2603 -2162 -77 -1177 -350; -585 789 1313 -23 -1477 931 -26 -626 1613 -974; 1545 -2450 1749 2206 722 -1631 -1794 -1301 -477 -1467; 1873 -740 -1856 -239 2235 1452 -1620 -1548 -644 2197; 2339 -370 -1760 -1875 -1040 1414 -1674 1040 2140 -734 ])], 524288)  == IntModQ.([79 1101 1153 1326 192 841 1315 956 789 1810; 49 1892 1441 996 85 1543 226 935 1045 167; 1293 1627 726 1924 923 665 1192 1532 24 81; 556 547 1490 887 1400 1474 1860 1856 90 43; 329 1364 1128 1780 1827 80 1331 1622 667 287; 1872 1192 202 1714 79 1329 1781 2027 1385 386; 94 1075 1596 1445 335 963 2005 1264 903 609; 296 1465 23 156 1435 114 1030 1827 639 146; 402 1245 1157 1799 1892 589 565 762 1791 1872; 1471 213 867 951 2012 1814 5 1566 631 1001]))
    
end

if test_matrix_product_mod_q_really_slow
    include("../test/stress_test/really_big_matrix_product_mod_.jl")
end

if test_matrix_product_really_slow
    print("degree 4, 64 by 64, product length 2^8: ")
    include("../test/stress_test/really_big_matrix_product_test.jl")
end
