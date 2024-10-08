# ==============================< PACKAGES >============================================================
using Printf

include("tInterval.jl")

const verbose = false

# ==============================< Declarations >========================================================

sinc(x::Float64)::Float64 = return (x == 0) ? (1) : (sin(x) / x)
sinc_d_1st(x::Float64)::Float64 = return (x == 0) ? (0) : ((x * cos(x) - sin(x)) / (x^2))
sinc_d_2nd(x::Float64)::Float64 = return (x == 0) ? (-0.33333333315479) : (((x^2 - 2) * -x * sin(x) - 2 * x^2 * cos(x)) / (x^4))

# ==============================< Newton >==============================================================

N(x::Float64)::Float64 = return (x - (sinc_d_1st(x) / sinc_d_2nd(x))) 

function newton_method(x::Float64, threshold::Float64 = 0.0000000000001)::Float64
    tmp ::Union{Float64, Nothing}   = x-1       # xk-1

    while !(tmp <= x <= tmp + threshold)
        tmp = x
        x = N(x)
    end

    return x
end

# ==============================< Newton Interval >=====================================================

N(f::Function, df::Function, x::tInterval{Float64}, c::Float64 = m(x)) = return c - f(c)/df(x)

θ(x::tInterval{Float64}, δ::Float64 = 1.1, λ::Float64 = 10^-12) = m(x) + δ * (x - m(x)) + λ * tInterval(-1, 1)

function NewtonInterval(x::tInterval, max_iter::Int64 = 1000, τ::Float64 = 10^-8)
    i::Int64 = 0
    while (i < max_iter) && (w(x) >= τ)
        x = N(sinc_d_1st, sinc_d_2nd, θ(x))
        i += 1
    end

    return x
end

# ==============================< FUNCTIONS >===========================================================
function approx_sinc_Newton(a::Float64, b::Float64)

    res_min::Float64 = -1.0
    res_max::Float64 = -1.0

    threshold::Float64 = 0.0000000000001 # accuracy
    
    # Case 1: a > b ------------------------------------------------------------------------------------
    if a > b

        verbose && println("a > b")
        return res_min, res_max # -1.0: sentinel value for indicate an impossible case

    # Case 2: a = b ------------------------------------------------------------------------------------
    elseif a==b
        if (a==0) && (b==0)
            verbose && println("a == b == 0")
            return 1., 1.
        end

        verbose && println("a == b != 0")
    
        res = (a == 0) ? (1) : (sin(a)/a) # res == min == max 
        return res, res
        
    # Case 3: 0 ∈ [a,b] --------------------------------------------------------------------------------
    elseif (a <= 0 <= b)
        verbose && println("a <= 0 <= b")

        res_max = 1
        
        min_sinc = newton_method(3π/2, threshold) # lowest point of sinc (first minimum ~ 3π/2 -> use newton method starting at x = 3π/2)

        if min_sinc < abs(a) || min_sinc < b # is a or b further than the first minimum (lowest point of the function)

            verbose && println("0 <= y <= b or 0 <= y <= abs(a) | sinc(y) -> lowes point of sinc")

            return sin(min_sinc)/min_sinc, res_max
        else
            fa = (a == 0) ? 1 : sin(a)/a
            fb = (b == 0) ? 1 : sin(b)/b

            res_min = min(fa, fb)
            
            return res_min, res_max
        end
    else
        # Case 4: a ≤ 0 et b ≤ 0 -----------------------------------------------------------------------
        if ( a <= 0 && b <= 0)
            verbose && println("BEFORE NEG: a = ", a, " et b = ", b)

            a, b = abs(b), abs(a) # transformation in positif number

            verbose && println("AFTER NEG: a = ", a, " et b = ", b)
        end
        
        # Case 5: a ≥ 0 et b ≥ 0 -----------------------------------------------------------------------
        verbose && println("POS: a = ", a, " et b = ", b)

        (b - a > 2π) && (b = a + 2π)

        ka::Float64 = ceil((a-π)/π)
        kb::Float64 = ceil((b-π)/π)

        if ka == kb # a and b are both in the range of 1 extremum
            extremum::Float64 = newton_method(((1+2ka)π)/2, threshold)

            fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
            fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

            if a <= extremum <= b # a and b are from both sides of the extremum 
                verbose && println("ka == kb; a <= extremum <= b")

                res_min = (ka%2 == 0) ? min(fa, fb) : sin(extremum)/extremum # if k even then extremum is a maximum, hence min = min(f(a), f(b)). Otherwise min = f(extremum).
                res_max = (ka%2 == 0) ? sin(extremum)/extremum : max(fa, fb) # if k even then extremum is a maximum, hence max = f(extremum). Otherwise max = max(f(a), f(b)). 

                return res_min, res_max
            else
                verbose && println("ka == kb; a, b <= extremum OR extremum <= a, b")

                res_min = min(fb, fa)
                res_max = max(fb, fa)

                return res_min, res_max
            end
        else
            extremum_ka::Float64 = newton_method(((1+2ka)π)/2, threshold)
            extremum_kb::Float64 = newton_method(((1+2kb)π)/2, threshold)

            if (a <= extremum_ka <= b <= extremum_kb) && (ka + 1 == kb)

                verbose && println("(a <= extremum_ka <= b <= extremum_kb) && (ka + 1 == kb)")

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
                fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

                res_min = (ka%2 == 0) ? fb : sin(extremum_ka)/extremum_ka
                res_max = (ka%2 == 0) ? ((ka == 0) ? fa : sin(extremum_ka)/extremum_ka) : fb

                return res_min, res_max

            elseif (a <= extremum_ka) && (((extremum_kb <= b) && (ka + 1 == kb)) || (ka + 2 == kb))

                verbose && println("(a <= extremum_ka) && (((extremum_kb <= b) && (ka + 1 == kb)) || (ka + 2 == kb))")

                (ka + 2 == kb) && (kb -= 1; extremum_kb = newton_method((((1+2kb)π)/2), threshold))

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)

                res_min = (ka%2 == 0) ? sin(extremum_kb)/extremum_kb : sin(extremum_ka)/extremum_ka
                res_max = (ka%2 == 0) ? ((ka == 0) ? fa : sin(extremum_ka)/extremum_ka) : sin(extremum_kb)/extremum_kb

                return res_min, res_max

            elseif (extremum_ka <= a <= b <= extremum_kb ) && (ka + 1 == kb) # a and b from different range but no extremum in between

                verbose && println("(extremum_ka <= a <= b <= extremum_kb ) && (ka + 1 == kb)")

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
                fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

                res_min = min(fb, fa)
                res_max = max(fb, fa)

                return res_min, res_max

            elseif (extremum_ka <= a <= extremum_kb <= b) && (ka + 1 == kb)

                verbose && println("(extremum_ka <= a <= extremum_kb <= b) && (ka + 1 == kb)")

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
                fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

                res_min = (kb%2 == 0) ? fa : sin(extremum_kb)/extremum_kb
                res_max = (kb%2 == 0) ? sin(extremum_kb)/extremum_kb : fa

                return res_min, res_max

            elseif (extremum_ka <= a <= b <= extremum_kb) && (ka + 2 == kb)

                verbose && println("(extremum_ka <= a <= b <= extremum_kb) && (ka + 2 == kb)")

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
                fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

                kc = ka + 1 # range between ka and kb
                extremum_kc = newton_method(((1+2kc)π)/2, threshold)

                res_min = (kc%2 == 0) ? min(fa, fb) : sin(extremum_kc)/extremum_kc
                res_max = (kc%2 == 0) ? sin(extremum_kc)/extremum_kc : max(fa, fb)

                return res_min, res_max

            elseif (extremum_ka <= a <= extremum_kb <= b) && (ka + 2 == kb)

                verbose && println("(extremum_ka <= a <= extremum_kb <= b) && (ka + 2 == kb)")

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
                fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

                kc = ka + 1 # range between ka and kb
                extremum_kc = newton_method(((1+2kc)π)/2, threshold)

                res_min = (kc%2 == 0) ? min(fa, sin(extremum_kb)/extremum_kb) : sin(extremum_kc)/extremum_kc
                res_max = (kc%2 == 0) ? sin(extremum_kc)/extremum_kc : max(fa, sin(extremum_kb)/extremum_kb)

                return res_min, res_max
                
            end
        end
    end
end

function approx_sinc_NewtonInterval(a::Float64, b::Float64)::Union{Nothing, tInterval{Float64}}

    res_min::Float64 = -1.0
    res_max::Float64 = -1.0
    
    # Case 1: a > b ------------------------------------------------------------------------------------
    if a > b

        verbose && println("a > b")
        return tInterval{typeof(a)}(nothing, nothing) # situation not handled ∅

    # Case 2: a = b ------------------------------------------------------------------------------------
    elseif a == b
        if (a==0) && (b==0)
            verbose && println("a == b == 0")
            return tInterval(1., 1.) # 1. == min == max
        else
            verbose && println("a == b != 0")
            tmp = (a == 0) ? (1) : (sin(a)/a) # tmp == min == max
            return tInterval(tmp, tmp)
        end
    # Case 3: 0 ∈ [a,b] --------------------------------------------------------------------------------
    elseif (a <= 0 <= b)
        verbose && println("a <= 0 <= b")

        res_max = 1
        
        # use Newton interval method to determine the lowest point of sinc. 
        sinc_minimum::tInterval{Float64} = NewtonInterval(tInterval(4., 3π/2)) # lowest point of sinc at least belong in [4, 3π/2] → graphic opservation
    
        if sinc_minimum.l ≤ abs(a) || sinc_minimum.l ≤ b
            # set the extremum as lowest image of [a, b]

            verbose && println("0 < y ≤ b or 0 < y ≤ |a| | sinc(y) -> lowest point of sinc. \nThus y belong in [a, b] and is the lowest image of [a, b]")

            res_min = min(sinc(sinc_minimum.l), sinc(sinc_minimum.u))

            return tInterval(res_min, res_max)
        else
            # extremum don't belong in [a, b]
            
            verbose && println("0 ≤ b < y and 0 ≤ |a| < y | sinc(y) -> lowest point of sinc. \nThus y don't belong in [a, b] and therfore the lowest image of [a, b] is min(sinc(a), sinc(b))")

            res_min = min(sinc(a), sinc(b))
            
            return tInterval(res_min, res_max)
        end
    else
        # Case 4: a < b ≤ 0 ----------------------------------------------------------------------------
        if ( a < b ≤ 0)
            verbose && println("a < b ≤ 0: \nReturning to a strictly positive situation → a = $a et b = $b become a = $(-b) et b = $(-a)")

            a, b = abs(b), abs(a) # return to a strictly positive situation
        end
        
        # Case 5: a ≥ 0 et b ≥ 0 -----------------------------------------------------------------------
        if (b - a > 2π + 0.5)
            # println("w([a, b]) > 2π + 0.5: \nBecause the variation of sinc only decrease as x∈[a, b] grow larger and the maximum distance between two extremum of sinc is smaller than 2π+0.5, the upper bound b could be set to b := a + 2π + 0.5")
            b = a + 2π + 0.5
        end

        ka::Float64 = ceil((a-π)/π)
        kb::Float64 = ceil((b-π)/π)

        if ka == kb # a and b are both in the range of 1 extremum
            extremum::tInterval{Float64} = NewtonInterval(tInterval(((1+2ka)π)/2 - 0.3, ((1+2ka)π)/2))

            fa = sinc(a) # compute f(a)
            fb = sinc(b) # compute f(b)

            if a ≤ extremum.l ≤ extremum.u ≤ b # a and b are from both sides of the extremum 
                verbose && println("ka == kb; a <= extremum <= b")
                tmp = (fa, fb, sinc(extremum.l), sinc(extremum.u))

                return tInterval(minimum(tmp), maximum(tmp))
            else
                verbose && println("ka == kb; a, b <= extremum OR extremum <= a, b")

                return tInterval(min(fb, fa), max(fb, fa))
            end
        else
            extremum_ka::tInterval{Float64} = (ka == 0) ? tInterval(sinc(a), sinc(a)) : NewtonInterval(tInterval(((1+2ka)π)/2 - 0.3, ((1+2ka)π)/2))
            extremum_kb::tInterval{Float64} = (kb == 0) ? tInterval(sinc(b), sinc(b)) : NewtonInterval(tInterval(((1+2kb)π)/2 - 0.3, ((1+2kb)π)/2))

            verbose && println("ka = $(extremum_ka), kb = $(extremum_kb), a = $a, b = $b, ka = $ka, kb = $kb")

            if (ka + 1 == kb)
                if (a ≤ extremum_ka.l && b ≤ extremum_kb.l)  # ane  ^   | b ^   |   ^ 
                    verbose && println("(a ≤ extremum_ka.l ≤ b ≤ extremum_kb.l) && (ka + 1 == kb)")

                    tmp = (sinc(a), sinc(b), sinc(extremum_ka.l), sinc(extremum_ka.u))

                    return tInterval(minimum(tmp), maximum(tmp))

                elseif (a ≤ extremum_ka.l ≤ extremum_kb.u ≤ b) # a ^   |   ^ b |   ^ 
                    verbose && println("(a <= extremum_ka) && (((extremum_kb <= b) && (ka + 1 == kb)) || (ka + 2 == kb))")

                    tmp = (sinc(a), sinc(b), sinc(extremum_ka.l), sinc(extremum_ka.u), sinc(extremum_kb.l), sinc(extremum_kb.u))

                    return tInterval(minimum(tmp), maximum(tmp))

                elseif (extremum_ka.u ≤ a ≤ b ≤ extremum_kb.l ) #   ^ a | b ^   |   ^ 
                    verbose && println("(extremum_ka <= a <= b <= extremum_kb ) && (ka + 1 == kb)")

                    tmp = (sinc(a), sinc(b))

                    return tInterval(minimum(tmp), maximum(tmp))

                elseif (extremum_ka.u ≤ a ≤ extremum_kb.u ≤ b) #   ^ a |   ^ b |   ^ 
                    verbose && println("(extremum_ka <= a <= extremum_kb <= b) && (ka + 1 == kb)")

                    tmp = (sinc(a), sinc(b), sinc(extremum_kb.l), sinc(extremum_kb.u))

                    return tInterval(minimum(tmp), maximum(tmp))
                end

            else # ka + 2 = kb
                kc = ka + 1 # range between ka and kb
                extremum_kc::tInterval{Float64} = NewtonInterval(tInterval(((1+2kc)π)/2 - 0.3, ((1+2kc)π)/2))

                if (a ≤ extremum_ka.l ≤ b ≤ extremum_kb.l) # a ^   |   ^   | b ^ 
                    verbose && println("(a ≤ extremum_ka.l ≤ b ≤ extremum_kb.l) && (ka + 2 == kb)")

                    tmp = (sinc(a), sinc(b), sinc(extremum_ka.l), sinc(extremum_ka.u), sinc(extremum_kc.l), sinc(extremum_kc.u))

                    return tInterval(minimum(tmp), maximum(tmp))

                elseif (a ≤ extremum_ka.l ≤ b ≤ extremum_kb.l) # a ^   |   ^   |   ^ b
                    verbose && println("(a ≤ extremum_ka.l ≤ b ≤ extremum_kb.l) && (ka + 2 == kb)")

                    tmp = (sinc(a), sinc(b), sinc(extremum_ka.l), sinc(extremum_ka.u), sinc(extremum_kc.l), sinc(extremum_kc.u))

                    return tInterval(minimum(tmp), maximum(tmp))

                elseif (extremum_ka.u ≤ a ≤ b ≤ extremum_kb.l) #   ^ a |   ^   | b ^ 
                    verbose && println("(extremum_ka.u ≤ a ≤ b ≤ extremum_kb.l) && (ka + 2 == kb)")

                    tmp = (sinc(a), sinc(b), sinc(extremum_kc.l), sinc(extremum_kc.u))

                    return tInterval(minimum(tmp), maximum(tmp))

                elseif (extremum_ka.u ≤ a ≤ extremum_kb.u ≤ b) #   ^ a |   ^   |   ^ b
                    verbose && println("(extremum_ka.u ≤ a ≤ extremum_kb.u ≤ b) && (ka + 2 == kb)")

                    tmp = (sinc(a), sinc(b), sinc(extremum_kb.l), sinc(extremum_kb.u), sinc(extremum_kc.l), sinc(extremum_kc.u))

                    return tInterval(minimum(tmp), maximum(tmp))
                    
                end
            end
        end
    end
end



# ==============================< MAIN >================================================================

function main()

    println("\nLaunch of the code...")

    # Assignment of values ​​to the x terminal of the interval
    a::Float64 = 1.2
    b::Float64 = 4.001

    borneInf, borneSup = approx_sinc_Newton(a,b)

    println("\nNewton method:")
    if (borneInf == -1) && (borneSup == -1)
        println("Not possible")
    else
        println("sinc([", a, ",", b, "]) = [",borneInf, ",", borneSup,"]" )
    end

    println("\nInterval Newton method:")

    res = approx_sinc_NewtonInterval(a, b)

    if res.l == nothing || res.u == nothing
        println("Not possible")
    else
        println("sinc([$a, $b]) = [$(round(res.l, digits=4)), $(round(res.u, digits=4))]")
    end

    println("\n...End of the code")

    return nothing
end
