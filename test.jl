function factorial(n::{<:Integer})
    if n == 0
        return 1
    end
    
    return factorial(n - 1)
end

println(factorial(5))