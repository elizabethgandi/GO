# ========================================== PACKAGES ===============================================================
using Printf

# ========================================== FONCTIONS ==============================================================
function approximation_fct_sinus_cardinal(a::Float64, b::Float64)

    min::Float64   = -1.0
    max::Float64   = -1.0
    borne::Float64 = 0.0
    
    # Cas 1: a > b --------------------------------------------------------------------------------------------------
    if a > b
        return min, max # -1.0 pour l'instant est une valeur sentinelle pour indiquer impossible 

    # Cas 2: a = b --------------------------------------------------------------------------------------------------
    elseif a==b
        if (a==0) && (b==0)
            return 1, 1
        end
    
        min = sin(a)/a
        max = min
        return min, max
        
    # Cas 3: 0 ∈ [a,b] ----------------------------------------------------------------------------------------------
    elseif (0 >= a && 0 <= b)
        max = 1

        if (b >= 3π/2) || (a <= -(3π/2))
            min = (sin((3π)/2))/((3π)/2)
        else
            println("En traitement ...")
            # Ici b < 3π/2 et a soit negatifs soit vaut 0

            min = -2
            #chercher le min
            # min = chercher_min(a,b)
        end

        return min, max # pour l'instant -2 remplace la routine pour chercher la valeur min
    else
        # Cas 4: a ≤ 0 et b ≤ 0 -------------------------------------------------------------------------------------
        if ( a <=0 && b <= 0)
            println("INIT NEG: a = ", a, " et b = ", b)

            a, b = abs(b), abs(a) # transformation en positif

            println("APRES NEG: a = ", a, " et b = ", b)
        end
        
        # Cas 5: a ≥ 0 et b ≥ 0 -------------------------------------------------------------------------------------

        # → dans ce cas 5 a, b ≠ 0 et strictement pos 
        
        println("POS: a = ", a, " et b = ", b)

        if (b-a <= 2π)
            nothing
        elseif (3π/2 >= a && 3π/2 <= b)
            min = (sin((3π)/2))/((3π)/2)
            max = 2 # valeur arbitraire
            # chercher le max
            # max = chercher_max(a,b)

        else
            b     = a+2π
            k     = ceil((a-(π/2))/(2π)) 
            borne = ((1+2k)π)/2

            (k%2 ==0) ? println("k paire!") : println("k impaire!")
            println("b = ", b, "k = ", k, " || borne = ", borne)

            if (a <= borne && b <= borne) || (a <= borne+(3π/2) && b <= borne+(3π/2) && a >= borne+π && b >= borne+π)
                max = sin(a)/a
                min = sin(b)/b
            
            elseif (a <= borne+π && b <= borne+π && a >= borne && b >= borne)
                min = sin(a)/a
                max = sin(b)/b

            else
                println("En traitement ...")
                # traitement de l'horreur 
                # min, max = approximation(a,b)
                nothing
            end

             println("min = ", min, " || max = ", max)
        end
    end

    borneInf, borneSup = min, max #temporaire

    return borneInf, borneSup
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