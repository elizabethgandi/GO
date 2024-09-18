
# function sinus cardinal
sinc(x::Float64)::Float64 = return (x == 0) ? (1) : (sin(x) / x)

# Derivee sinus cardinal
sinc_derive_premiere(x::Float64)::Float64 = return (x == 0) ? (0) : ((x * cos(x) - sin(x)) / (x^2))

# Derivee seconde sinus cardinal
sinc_derive_seconde(x::Float64)::Float64 = return (x == 0) ? (-0.33333333315479) : (((2 - x^2) * x * sin(x) - 2 * x^2 * cos(x)) / (x^4))

# Methode de newton
N(x::Float64)::Float64 = return (x - (sinc_derive_premiere(x) / sinc_derive_seconde(x)))

