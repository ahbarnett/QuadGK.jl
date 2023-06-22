# AHB added notes/docs 6/18/23

# integration segment (a,b), estimated integral I, and estimated error E.
# Segment is templated by 3 types, that of (a,b), that of I, and that of E:
struct Segment{TX,TI,TE}
    a::TX
    b::TX
    I::TI
    E::TE
    #fvals::Vector{TI}    # ahb  for speeding up later once works ?
end
# overload some comparison ops on two Segments:
Base.@pure Base.promote_rule(::Type{Segment{TX,TI,TE}}, ::Type{Segment{TX′,TI′,TE′}}) where {TX,TI,TE,TX′,TI′,TE′} =
    Segment{promote_type(TX,TX′), promote_type(TI,TI′), promote_type(TE,TE′)}
Base.convert(::Type{T}, s::Segment) where {T<:Segment} = T(s.a,s.b,s.I,s.E)
# make segments sort by error in a heap...
Base.isless(i::Segment, j::Segment) = isless(i.E, j.E)

# Internal routine: approximately integrate f(x) over the interval (a,b)
# by evaluating the single-segment integration rule (x,w,gw).
# nrm is a norm used to measure size of Ik-Ig (allows non-scalar).
# Return a Segment which contains (a,b,Ik,E), all of x's type.
# type F of f appears unused.
function evalrule(f::F, a,b, x,w,gw, nrm) where {F}
    # Ik and Ig are integrals via Kronrod and Gauss rules, respectively
    # Type of x controls all the types in this routine.
    # s is factor to scale [-1,1] to our interval [a,b]...
    s = convert(eltype(x), 0.5) * (b-a)
    n1 = 1 - (length(x) & 1) # 0 if even order, 1 if odd order
    # apply Kronrod (all nodes, w), and embedded Gauss rule (even nodes, gw)
    if n1 == 0 # even: Gauss rule does not include x == 0
        Ik = f(a + s) * w[end]
        Ig = zero(Ik)
    else # odd: don't count x==0 twice in Gauss rule
        f0 = f(a + s)
        Ig = f0 * gw[end]
        Ik = f0 * w[end] +
            (f(a + (1+x[end-1])*s) + f(a + (1-x[end-1])*s)) * w[end-1]
    end
    for i = 1:length(gw)-n1
        fg = f(a + (1+x[2i])*s) + f(a + (1-x[2i])*s)
        fk = f(a + (1+x[2i-1])*s) + f(a + (1-x[2i-1])*s)
        Ig += fg * gw[i]
        Ik += fg * w[2i] + fk * w[2i-1]
    end
    # Ik,Ig have type of f right now
    Ik_s, Ig_s = Ik * s, Ig * s # new variable since this may change the type
    E = nrm(Ik_s - Ig_s)
    if isnan(E) || isinf(E)
        throw(DomainError(a+s, "integrand produced $E in the interval ($a, $b)"))
    end
    return Segment(oftype(s, a), oftype(s, b), Ik_s, E)
end

# compute result = f(x1) + f(x2) in-place
function eval2x!(result, f::InplaceIntegrand{F}, x1, x2) where {F}
    f.f!(result, x1)
    f.f!(f.fx, x2)
    result .+= f.fx
end

# as above, but call assume a mutable result type (e.g. an array) and
# act in-place using `f!(result, x)`.
function evalrule(f::InplaceIntegrand{F}, a,b, x,w,gw, nrm) where {F}
    # Ik and Ig are integrals via Kronrod and Gauss rules, respectively
    s = convert(eltype(x), 0.5) * (b-a)
    n1 = 1 - (length(x) & 1) # 0 if even order, 1 if odd order
    fg, fk, Ig, Ik = f.fg, f.fk, f.Ig, f.Ik # pre-allocated temporary arrays
    f.f!(f.fx, a + s)
    if n1 == 0 # even: Gauss rule does not include x == 0
        Ik .= f.fx .* w[end]
        Ig .= zero.(Ik)
    else # odd: don't count x==0 twice in Gauss rule
        Ig .= f.fx .* gw[end]
        f.f!(f.fg, a + (1+x[end-1])*s)
        f.f!(f.fk, a + (1-x[end-1])*s)
        Ik .= f.fx .* w[end] .+ (f.fg + f.fk) .* w[end-1]
    end
    for i = 1:length(gw)-n1
        eval2x!(fg, f, a + (1+x[2i])*s, a + (1-x[2i])*s)
        eval2x!(fk, f, a + (1+x[2i-1])*s, a + (1-x[2i-1])*s)
        Ig .+= fg .* gw[i]
        Ik .+= fg .* w[2i] .+ fk .* w[2i-1]
    end
    Ik_s = Ik * s # new variable since this may change the type
    f.Idiff .= Ik_s .- Ig .* s
    E = nrm(f.Idiff)
    if isnan(E) || isinf(E)
        throw(DomainError(a+s, "integrand produced $E in the interval ($a, $b)"))
    end
    return Segment(oftype(s, a), oftype(s, b), Ik_s, E)
end
