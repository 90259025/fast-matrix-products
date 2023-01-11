# Returns a brightness for the complex number z in the Mandelbrot set. Optionally allows c to be set to compute the brightness for a Julia set instead.
function mandelbrotbrightness(z; c = z, maxiterations = 100)
    brightness = 0
    
    current_z = z
    
    for iteration = 1:maxiterations
        current_z = current_z^2 + c
        
        if abs(current_z) >= 2
            brightness = 0
            break
        end
    end
    
    return brightness
end