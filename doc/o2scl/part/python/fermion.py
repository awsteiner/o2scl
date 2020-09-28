import sympy

t=sympy.Symbol('t')
k=sympy.Symbol('k')
x=sympy.Symbol('x')

print('------------------------------------------------------')
print('Non-degenerate expansion')
print('')

print('Pressure for particles:')
Pk=t**2*(-1)**(k+1)/k/k*sympy.exp(k*(x+1)/t)*sympy.besselk(2,k/t)
Pka=(2*t**2*(-1)**(k+1)/k/k*sympy.cosh(k*(x+1)/t)*
     sympy.besselk(2,k/t))
print('Pk = '+str(Pk))
print('')

Pka2=(Pk+t**2*(-1)**(k+1)/k/k*sympy.exp(k*(-x-1)/t)*
      sympy.besselk(2,k/t))
Pka3=Pka2.subs(sympy.exp(k*(x+1)/t),2*sympy.cosh(k*(x+1)/t)-
               sympy.exp(-k*(x+1)/t))
print(sympy.simplify(Pka3)-Pka)
print('')
print('Pressure for antiparticles:')
print('Pka = '+str(Pka))
print('')

print('Density for particles:')
nk=sympy.simplify(sympy.diff(Pk,x))
print('nk/Pk = '+str(sympy.simplify(nk/Pk)))
print('')

print('Density for antiparticles:')
nka=sympy.simplify(sympy.diff(Pka,x))
print('nka/Pka = '+str(sympy.simplify(nka/Pka)))
print('')

print('Entropy for particles:')
sk=sympy.simplify(sympy.diff(Pk,t))
print('sk = '+str(sk))
print('')

print('Alternate expression for entropy:')
# Temporarily undefine nk
nk_old=nk
nk=sympy.Symbol('nk')
skalt=(nk*(4*t-k*x-k)/k/t+
       (-1)**(k+1)/k*sympy.exp(k*(x+1)/t)*
       sympy.besselk(1,k/t))
print('skalt = '+str(skalt))
print('')

# Restore the definition of nk
nk=nk_old
skalt=skalt.subs('nk',nk)
print('skalt-sk = '+str(sympy.simplify(skalt-sk)))
print('')

print('Entropy for antiparticles:')
ska=sympy.simplify(sympy.diff(Pka,t))
print('ska = '+str(ska))
print('')

print('Alternate expression for entropy of antiparticles:')
nka_old=nka
nka=sympy.Symbol('nka')
skaalt=(-nka*(1+x)/t+2*(-1)**(k+1)/k*
        sympy.cosh(k*(x+1)/t)*sympy.besselk(3,k/t))
print('skaalt = '+str(skaalt))
print('')

nka=nka_old
skaalt=skaalt.subs('nka',nka)
print('skaalt-ska = '+str(sympy.simplify(skaalt-ska)))
print('')

print('Energy density for particles:')
edk=sympy.simplify(-Pk+t*sk+(x+1)*nk)
# Replace K_0 with a combination of K_1 and K_2
edk=edk.subs(sympy.besselk(0,k/t),sympy.besselk(2,k/t)-
             2*t/k*sympy.besselk(1,k/t))
print('edk = '+str(edk))
print('')

print('Energy density for antiparticles:')
edka=sympy.simplify(-Pka+t*ska-(x+1)*nka)
edka=edka.subs(sympy.besselk(0,k/t),sympy.besselk(2,k/t)-
               2*t/k*sympy.besselk(1,k/t))
edka=edka.subs(sympy.besselk(1,k/t),sympy.besselk(3,k/t)-
               4*t/k*sympy.besselk(2,k/t))
t1=edka.subs(sympy.besselk(2,k/t),sympy.Symbol('b2'))
t2=t1.subs(sympy.besselk(3,k/t),sympy.Symbol('b3'))
print('edka = '+str(sympy.simplify(t2).collect('b2')))
print('')

