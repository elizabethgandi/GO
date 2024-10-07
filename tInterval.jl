struct tInterval{T}
    l::Union{T, Nothing}
    u::Union{T, Nothing}
end

const verbose = false

# ====================< Base operation >====================

Base.:+(x::tInterval{T}, v::U) where {T, U}             = return tInterval{T}(x.l + v, x.u + v)         # [a, b] + c = [a + c, b + c]
Base.:+(v::U, x::tInterval{T}) where {T, U}             = return x + v                                  # c + [a, b] = [a + c, b + c]
Base.:+(x::tInterval{T}, y::tInterval{T}) where T       = return tInterval{T}(x.l + y.l, x.u + y.u)     # [a, b] + [c, d] = [a + c, b + d]

Base.:-(x::tInterval{T}, v::U) where {T, U}             = return tInterval{T}(x.l - v, x.u - v)         # [a, b] - c = [a - c, b - c]
Base.:-(v::U, x::tInterval{T}) where {T, U}             = return tInterval{T}(v - x.u, v - x.l)         # c - [a, b] = [c - b, c - a]
Base.:-(x::tInterval{T}, y::tInterval{T}) where T       = return tInterval{T}(x.l - y.u, x.u - y.l)     # [a, b] - [c, d] = [a - d, b - c]

function Base.:*(x::tInterval{T}, v::U) where {T, U}                                                    # [a, b] * c
    tmp::Tuple{Float64, Float64} = (x.l * v, x.u * v)                                                   # if c > 0 then [a * c, b * c]
    return tInterval(minimum(tmp), maximum(tmp))                                                        # else [b * c, a * c]
end
Base.:*(v::U, x::tInterval{T}) where {T, U} = return x * v                                              # c * [a, b] = [a, b] * c
function Base.:*(x::tInterval{T}, y::tInterval{T}) where T                                              # [a, b] * [b * c] = [min(ac, ad, bc, bd), max(ac, ad, bc, bd)]
    tmp::NTuple{4, T} = (x.l * y.l, x.l * y.u, x.u * y.l, x.u * y.u)
    return tInterval{T}(minimum(tmp), maximum(tmp))
end

function Base.:inv(x::tInterval{T}) where T                                                             # [a, b]^-1 = 
    if x.l == x.u == 0                                                                                  # if        a = b = 0       then ∅
        return tInterval(nothing, nothing)                                                              # elseif    b < 0           then [b^-1, a^-1]
    elseif x.u < 0                                                                                      # elseif    0 < a           then [b^-1, a^-1]
        return tInterval((x.l)^-1, (x.u)^-1)                                                            # elseif    a < 0 AND b = 0 then [-∞, a^-1]
    elseif x.l > 0                                                                                      # elseif    a = 0 AND b > 0 then [b^-1, +∞]
        return tInterval((x.u)^-1, (x.l)^-1)                                                            # else                           [-∞, +∞]
    elseif x.l < 0 && x.u == 0
        return tInterval(-Inf, (x.l)^-1)
    elseif x.l == 0 && x.u > 0
        return tInterval((x.u)^-1, Inf)
    else
        return tInterval(-Inf, Inf)
    end
end

function Base.:^(x::tInterval{T}, p::Int64) where T                                                     # [a, b]^c =  
    if p < 0                                                                                            # if        c < 0           then ([a, b]^-1)^c
        return (x^-1)^-p                                                                                # else
    else                                                                                                #  |  if   c ≡ 0 mod(1) AND c ≥ 3
        if p%2 == 1                                                                                     #  |  then [a^c, b^c]
            return tInterval((x.l)^p, (x.u)^p)                                                          #  |  else
        else                                                                                            #  |   |  if        b < 0 then [b^p, a^p]
            if x.u < 0                                                                                  #  |   |  else if   0 < a then [a^p, b^p]
                return tInterval((x.u)^p, (x.l)^p)                                                      #  |   |  else                 [0, max(a^p, b^p)]
            elseif x.l > 0                                                                              #  |  endif
                return tInterval((x.l)^p, (x.u)^p)                                                      # endif
            else
                return tInterval(0., max(x.l^p, x.u^p))
            end
        end
    end
end

Base.:/(x::tInterval{T}, y::tInterval{T}) where T = x * (y^(-1))                                        # [a, b] / [c, d] = [a, b] * ([c, d]^-1)
Base.:/(v::U, x::tInterval{T}) where {T, U} = v * x^-1                                                  # c / [a, b] =  c * ([a, b]^-1)
function Base.:/(x::tInterval{T}, v::U) where {T, U}                                                    # [a, b] / c = [min(a/c, b/c), max(a/c, b/c)]
    @assert v == 0 "Error in Base.:/(x::tInterval{T}, v::U) -> v == 0 undefined!"
    tmp::Tuple{Float64, Float64} = (x.l / v, x.u / v)
    return tInterval(minimum(tmp), maximum(tmp))
end

# ====================< Operation >====================

w(x::tInterval{T}) where T = x.u - x.l
m(x::tInterval{T}) where T = (x.u + x.l)/2
r(x::tInterval{T}) where T = (x.u - x.l)/2

# ====================< Set operation >===================='

∈(c::U, x::tInterval{T}) where {U, T} = x.l <= c <= x.u
∉(c::U, x::tInterval{T}) where {U, T} = !c∈x
Base.:==(x::tInterval{T}, y::tInterval{T}) where T = (x.l == y.l) && (x.u == y.u)
Base.:!=(x::tInterval{T}, y::tInterval{T}) where T = (x.l != y.l) || (x.u != y.u)
⊆(x::tInterval{T}, y::tInterval{T}) where T = y.l <= x.l <= x.u <= y.u
⊂(x::tInterval{T}, y::tInterval{T}) where T = y.l < x.l < x.u < y.u
∩(x::tInterval{T}, y::tInterval{T}) where T = tInterval{T}(max(x.l, y.l), min(x.u, y.u))

# ====================< Hausdorff distance >====================

d(x::tInterval{T}, y::tInterval{T}) where T = max(abs(x.l - y.l), abs(x.u - y.u))

# ====================< Elementary functions >====================

Base.:exp(x::tInterval{T}) where T = tInterval{T}(exp(x.l), exp(x.u))

function Base.:log(x::tInterval{T}) where T
    if x.u <= 0
        return tInterval(nothing, nothing) 
    elseif a > 0
        return tInterval(log(x.l), log(x.u))
    else
        return tInterval(-Inf, log(x.u))
    end
end

function Base.:cos(x::tInterval{T}) where T
    if w(x) ≥ 2π
        return tInterval{T}(-1, 1)
    else
        cos_l::Float64 = cos(x.l)
        cos_u::Float64 = cos(x.u)

        k_l::Float64 = ceil(x.l/π)
        k_u::Float64 = ceil(x.u/π)

        if k_l == k_u
            # ̲x and ̄x are on the same slope
            return tInterval{T}(min(cos_l, cos_u), max(cos_l, cos_u))
        elseif k_l + 1 == k_u
            # ̲x and ̄x have one extremum in between 
            extremum        ::Float64 = π * k_l
            cos_extremum    ::Float64 = cos(extremum)

            return tInterval{T}(min(cos_l, cos_u, cos_extremum), max(cos_l, cos_u, cos_extremum))
        else
            # ̲x and ̄x have two extremum in between (a max and a min)
            return tInterval{T}(-1, 1)
        end 
    end
end

function Base.:sin(x::tInterval{T}) where T
    if w(x) ≥ 2π
        return tInterval{T}(-1, 1)
    else
        sin_l::Float64 = sin(x.l)
        sin_u::Float64 = sin(x.u)

        k_l::Float64 = ceil((x.l - (π/2))/π)
        k_u::Float64 = ceil((x.u - (π/2))/π)

        if k_l == k_u
            # ̲x and ̄x are on the same slope
            return tInterval{T}(min(sin_l, sin_u), max(sin_l, sin_u))
        elseif k_l + 1 == k_u
            # ̲x and ̄x have one extremum in between 
            extremum        ::Float64 = π * k_l + π/2
            sin_extremum    ::Float64 = sin(extremum)

            return tInterval{T}(min(sin_l, sin_u, sin_extremum), max(sin_l, sin_u, sin_extremum))
        else
            # ̲x and ̄x have two extremum in between (a max and a min) 
            return tInterval{T}(-1, 1)
        end
    end
end

function sinc(x::tInterval{T}) where T
    verbose && println("Warning sinc(x::tInterval{T}): this method rely on basic interval operators wich tend to provide very pesimistics results")
    return sin(x)/x
end

function sinc_d_1st(x::tInterval{T}) where T
    verbose && println("Warning sinc_d_1st(x::tInterval{T}): this method rely on basic interval operators wich tend to provide very pesimistics results")
    return (x * cos(x) - sin(x)) / x^2
end

function sinc_d_2nd(x::tInterval{T}) where T
    verbose && println("Warning sinc_d_2nd(x::tInterval{T}): this method rely on basic interval operators wich tend to provide very pesimistics results")
    return ((2 - x^2) * x * sin(x) - 2 * x^2 * cos(x))/x^4
end
