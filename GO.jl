using JuMP, HiGHS, Printf


function approximation_fct_sinus_cardinal(a::Float64, b::Float64)

    min::Float64 = -1.0
    max::Float64 = -1.0
    
    # Cas 1: a > b ------------------------------------------------------------------------------------
    if a > b
        return min, max # -1.0 pour l'instant est une valeur sentinelle pour indiquer impossible 

    # Cas 2: a = b ------------------------------------------------------------------------------------
    elseif a==b
        if (a==0) && (b==0)
            return 1, 1
        end
    
        min = sin(a)/a
        max = min
        return min, max
        
    # Cas 3: 0 ∈ [a,b] ------------------------------------------------------------------------------------
    elseif (0 >= a && 0 <= b)
        max = 1

        #chercher le min
        # min = chercher_min(a,b)

        return -2, max # pour l'instant -2 remplace la routine pour chercher la valeur min
    else
        # Cas 4: a ≤ 0 et b ≤ 0 ------------------------------------------------------------------------------------
        if ( a <=0 && b <= 0)
            println("INIT NEG: a = ", a, " et b = ", b)

            a, b = abs(b), abs(a) # transformation en positif

            println("APRES NEG: a = ", a, " et b = ", b)
        end
        
        # Cas 5: a ≥ 0 et b ≥ 0 ------------------------------------------------------------------------------------
        
        println("POS: a = ", a, " et b = ", b)

        # min, max = approximation(a,b)
    

    end
    
    
    #@show min, max

    borneInf, borneSup = min, max #temporaire

    return borneInf, borneSup
end

function main()

    # Affectation des valeurs au bornes x de l'intevalle
    a::Float64 = -8.0
    b::Float64 = -6.0

    # 
    borneInf, borneSup = approximation_fct_sinus_cardinal(a,b)

    #borneInf, borneSup = a, b

    #println("hello world")

    #TODO si bornes = -1 alors impossible
    if (borneInf == -1) && (borneSup == -1)
        println("Impossible")
    else
        println("L'image des bornes [", a, ",", b, "] par la fonction sinus cardinale est l'intevalle [",borneInf, ",", borneSup,"]." )
    end

    return nothing
end