# play around with quadgk, including plotting, variant for meromorphic.
# Barnett 6/18/23
using QuadGK
import QuadGK.Segment
using Printf
# linux-only, maybe?
using Gnuplot


a = 0.0
b = 2π
I,E, segs, nev = quadgk_count(x -> exp(x), a,b, rtol=1e-16);
Ie = exp(2π)-1.0
abs((I-Ie)/Ie)

plot(seg::Segment) = @gp real([seg.a seg.b]) imag([seg.a seg.b]) "w lp pt 1 lc rgb '#000000' notit"
plot!(seg::Segment) = @gp :- real([seg.a seg.b]) imag([seg.a seg.b]) "w lp pt 1 lc rgb '#000000' notit"
function plot(segs::Vector{Segment{TX,TI,TE}}) where {TX,TI,TE}
    for (j,s) in enumerate(segs)
        # starts new fig or adds to current fig
        j==1 ? plot(s) : plot!(s)
        mid = (s.a+s.b)/2
        x,y = real(mid),imag(mid)   # hack to use gnuplot set label (text x,y)
        #@gp :- "set label '$j' at $x,$y"
        #Imag = abs(s.I)
        #@gp :- "set label '$Imag' at $x,$y"
    end
    @gp :- "set size ratio -1"
end

function showprob(pts,segs)
    plot(segs)
    # add interval endpts
    @gp :- real(pts) imag(pts) "w p pt 7 notit"
end
showprob([a b],segs)

#Gnuplot.quit(:default)

# test find_near_poles
s = Segment(0.0,2π,0.0,0.0)    # [0,2π] interval
d=1e-1;      # pole dist
z0 = π+1im*d;    # pole at pi+i.d
f(x) = 1/(x-z0);
tr, zr = find_near_poles(f, s)
@printf "single pole err %.3g\n" abs(zr[1]-z0)
tr, zr = find_near_poles(f, s, 0.03)       # rho so small cuts out that pole
@assert isempty(tr)

f(x) = (x+3)/(x-z0);   # zero at -3, starts to make it suck
tr, zr = find_near_poles(f, s)
@printf "nearby zero case: single pole err %.3g\n" abs(zr[1]-z0)

tr, zr = find_near_poles( x->1/(x-π+1im*d), s)



# usual adative way...
I,E, segs, nev = quadgk_count(x -> 1/(x-z0), a,b, rtol=1e-15);
Ie = log(b-z0) - log(a-z0); abs((I-Ie)/Ie)
showprob([a b z0],segs)
# new way...
opts = QuadGK.Meroopts(0.5);
I,E, segs, nev = quadgk_count(x -> 1/(x-z0), a,b, rtol=1e-12, mero=opts);
Ie = log(b-z0) - log(a-z0); abs((I-Ie)/Ie)
showprob([a b z0],segs)

d = 1e-3
z0 = π-1im*d;    # pole at pi+i.d
I,E, segs, nev = quadgk_count(x -> 1/(x-z0), a,b, rtol=1e-12);
length(segs)
I,E, segs, nev = quadgk_count(x -> 1/(x-z0), a,b, rtol=1e-12, mero=opts);
length(segs)
Ie = log(b-z0) - log(a-z0); abs((I-Ie)/Ie)
showprob([a b z0],segs)

