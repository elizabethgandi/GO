# ========================================== PACKAGES ===============================================================
using Printf

# sinc(x::Float64)::Float64 = return (x == 0) ? (1) : (sin(x) / x)
sinc_derive_premiere(x::Float64)::Float64 = return (x == 0) ? (0) : ((x * cos(x) - sin(x)) / (x^2))
sinc_derive_seconde(x::Float64)::Float64 = return (x == 0) ? (-0.33333333315479) : (((x^2 - 2) * -x * sin(x) - 2 * x^2 * cos(x)) / (x^4))

N(x::Float64)::Float64 = return (x - (sinc_derive_premiere(x) / sinc_derive_seconde(x))) 

function newton_method(x::Float64, threshold::Float64 = 0.0000000000001)::Float64
    tmp ::Union{Float64, Nothing}   = x-1       # xk-1

    while !(tmp <= x <= tmp + threshold)
        tmp = x
        x = N(x)
    end

    return x
end

# ========================================== FONCTIONS ==============================================================
function approximation_fct_sinus_cardinal(a::Float64, b::Float64)

    res_min::Float64 = -1.0
    res_max::Float64 = -1.0

    threshold::Float64 = 0.0000000000001 # accuracy
    
    # Cas 1: a > b ------------------------------------------------------------------------------------
    if a > b

        println("a > b")
        return res_min, res_max # -1.0 pour l'instant est une valeur sentinelle pour indiquer impossible 

    # Cas 2: a = b ------------------------------------------------------------------------------------
    elseif a==b
        if (a==0) && (b==0)
            println("a == b == 0")
            return 1., 1.
        end

        println("a == b != 0")
    
        res = (a == 0) ? (1) : (sin(a)/a) # res == min == max é
        return res, res
        
    # Cas 3: 0 ∈ [a,b] ------------------------------------------------------------------------------------
    elseif (a <= 0 <= b)
        println("a <= 0 <= b")

        res_max = 1

        # chercher le min
        
        min_sinc = newton_method(3π/2, threshold) # lowest point of sinc (first minimum ~ 3π/2 -> use newton method starting at x = 3π/2)

        if min_sinc < abs(a) || min_sinc < b # is a or b further than the first minimum (lowest point of the function)

            println("0 <= y <= b or 0 <= y <= abs(a) | sinc(y) -> lowes point of sinc")

            return sin(min_sinc)/min_sinc, res_max
        else
            fa = (a == 0) ? 1 : sin(a)/a
            fb = (b == 0) ? 1 : sin(b)/b

            res_min = min(fa, fb)
            
            return res_min, res_max
        end
    else
        # Cas 4: a ≤ 0 et b ≤ 0 ------------------------------------------------------------------------------------
        if ( a <= 0 && b <= 0)
            println("INIT NEG: a = ", a, " et b = ", b)

            a, b = abs(b), abs(a) # transformation en positif

            println("APRES NEG: a = ", a, " et b = ", b)
        end
        
        # Cas 5: a ≥ 0 et b ≥ 0 ------------------------------------------------------------------------------------
        
        println("POS: a = ", a, " et b = ", b)

        (b - a > 2π) && (b = a + 2π)

        ka::Float64 = ceil((a-π)/π)
        kb::Float64 = ceil((b-π)/π)
        
        

        if ka == kb # a and b are both in the range of 1 extremum
            extremum::Float64 = newton_method(((1+2ka)π)/2, threshold)

            fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
            fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

            if a <= extremum <= b # a and b are from both sides of the extremum 
                println("ka == kb; a <= extremum <= b")

                res_min = (ka%2 == 0) ? min(fa, fb) : sin(extremum)/extremum # if k even then extremum is a maximum, hence min = min(f(a), f(b)). Otherwise min = f(extremum).
                res_max = (ka%2 == 0) ? sin(extremum)/extremum : max(fa, fb) # if k even then extremum is a maximum, hence max = f(extremum). Otherwise max = max(f(a), f(b)). 

                return res_min, res_max
            else
                println("ka == kb; a, b <= extremum OR extremum <= a, b")

                res_min = min(fb, fa)
                res_max = max(fb, fa)

                return res_min, res_max
            end
        else
            extremum_ka::Float64 = newton_method(((1+2ka)π)/2, threshold)
            extremum_kb::Float64 = newton_method(((1+2kb)π)/2, threshold)

            if (a <= extremum_ka <= b <= extremum_kb) && (ka + 1 == kb)

                println("(a <= extremum_ka <= b <= extremum_kb) && (ka + 1 == kb)")

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
                fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

                res_min = (ka%2 == 0) ? fb : sin(extremum_ka)/extremum_ka
                res_max = (ka%2 == 0) ? ((ka == 0) ? fa : sin(extremum_ka)/extremum_ka) : fb

                return res_min, res_max
            elseif (a <= extremum_ka) && (((extremum_kb <= b) && (ka + 1 == kb)) || (ka + 2 == kb))

                println("(a <= extremum_ka) && (((extremum_kb <= b) && (ka + 1 == kb)) || (ka + 2 == kb))")

                (ka + 2 == kb) && (kb -= 1; extremum_kb = newton_method((((1+2kb)π)/2), threshold))

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)

                res_min = (ka%2 == 0) ? sin(extremum_kb)/extremum_kb : sin(extremum_ka)/extremum_ka
                res_max = (ka%2 == 0) ? ((ka == 0) ? fa : sin(extremum_ka)/extremum_ka) : sin(extremum_kb)/extremum_kb

                return res_min, res_max

            elseif (extremum_ka <= a <= b <= extremum_kb ) && (ka + 1 == kb) # a and b from different range but no extremum in between

                println("(extremum_ka <= a <= b <= extremum_kb ) && (ka + 1 == kb)")

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
                fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

                res_min = min(fb, fa)
                res_max = max(fb, fa)

                return res_min, res_max

            elseif (extremum_ka <= a <= extremum_kb <= b) && (ka + 1 == kb)

                println("(extremum_ka <= a <= extremum_kb <= b) && (ka + 1 == kb)")

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
                fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

                res_min = (kb%2 == 0) ? fa : sin(extremum_kb)/extremum_kb
                res_max = (kb%2 == 0) ? sin(extremum_kb)/extremum_kb : fa

                return res_min, res_max

            elseif (extremum_ka <= a <= b <= extremum_kb) && (ka + 2 == kb)

                println("(extremum_ka <= a <= b <= extremum_kb) && (ka + 2 == kb)")

                fa = (a == 0) ? 1 : sin(a)/a # compute f(a)
                fb = (b == 0) ? 1 : sin(b)/b # compute f(b)

                kc = ka + 1 # range between ka and kb
                extremum_kc = newton_method(((1+2kc)π)/2, threshold)

                res_min = (kc%2 == 0) ? min(fa, fb) : sin(extremum_kc)/extremum_kc
                res_max = (kc%2 == 0) ? sin(extremum_kc)/extremum_kc : max(fa, fb)

                return res_min, res_max

            elseif (extremum_ka <= a <= extremum_kb <= b) && (ka + 2 == kb)

                println("(extremum_ka <= a <= extremum_kb <= b) && (ka + 2 == kb)")

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




# ============================================== MAIN()==============================================================

function main()

    # Affectation des valeurs au bornes x de l'intevalle
    a::Float64 = 2.0
    b::Float64 = 17.0

    borneInf, borneSup = approximation_fct_sinus_cardinal(a,b)

    if (borneInf == -1) && (borneSup == -1)
        println("Impossible")
    else
        println("L'image des bornes [", a, ",", b, "] par la fonction sinus cardinale est l'intevalle [",borneInf, ",", borneSup,"]." )
    end

    return nothing
end
