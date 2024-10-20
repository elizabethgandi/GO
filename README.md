# Optimisation globale

### Vue d'ensemble

Le fichier principal, **GO.jl**, est conçu pour développer et implémenter des algorithmes d'optimisation, à commencer par l'algorithme de Newton. L'objectif est de déterminer les valeurs minimales et maximales de la fonction sinc(x) (sinus cardinal) sur un intervalle donné [a, b]. L'algorithme final vise à fournir une estimation aussi précise que possible de l'image de cet intervalle par la fonction sinc(x), tout en optimisant le temps de calcul.

Dans **tInterval.jl**, nous implémentons la classe `tInterval` du language Julia pour effectuer des opérations fiables sur des intervalles et ainsi réduire les erreurs d'arrondi tout en contrôlant manuellement la précision des calculs.

### Execution

Depuis un REPL Julia, vous pouvez exécuter le programme en entrant la commande suivante : include("GO.jl") 