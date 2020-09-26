import sympy

N1=50
N2=10

alpha=sympy.Symbol('alpha')
eta=sympy.Symbol('eta')
# Basic expression to approximate
cbt=(alpha**(sympy.Rational(-1,6))*(-1+sympy.sqrt(1+alpha))**
     (sympy.Rational(1,3)))
f=(1/cbt-cbt)
# Rewrite in terms of eta = alpha^(1/6)
f2=f.subs(alpha,eta**6)

if True:
    # Compute the series
    ser=sympy.series(f2,eta,0,N1)
    # Construct symbolic versions of the cbrt(2) and 2^(2/3)
    to3=sympy.Rational(2,1)**sympy.Rational(1,3)
    tt3=sympy.Rational(2,1)**sympy.Rational(2,3)

    # Just manually format the first two terms
    print('fp_t term=t03/eta;')
    print('fp_t term=-eta/to3;')

    # Print out the code for each term
    for i in range(5,N1):
        coeff=ser.coeff(eta,i)
        if i%6==1:
            coeff=sympy.simplify(coeff*to3)
            numer,denom=sympy.fraction(coeff)
            print('fp_t numer='+str(numer)+';')
            print('fp_t denom='+str(denom)+';')
            print('fp_t term=numer/denom*pow(eta,'+str(i)+')/to3;')
        else:
            coeff=sympy.simplify(coeff*tt3)
            numer,denom=sympy.fraction(coeff)
            if numer>0:
                print('fp_t numer='+str(numer)+';')
                print('fp_t denom='+str(denom)+';')
                print('fp_t term=numer/denom*pow(eta,'+str(i)+')/tt3;')
    print(' ')

if True:

    # Rewrite in terms of eta = alpha^(1/2)
    f2=f.subs(alpha,eta**2)
    
    ser=sympy.series(f2,eta,sympy.oo,N2)
    print(ser)

    # Print out the code for each term
    for i in range(-1,-N2,-1):
        coeff=ser.coeff(eta,i)
        numer,denom=sympy.fraction(coeff)
        if numer!=0:
            print('fp_t numer='+str(numer)+';')
            print('fp_t denom='+str(denom)+';')
            print('fp_t term=numer/denom*pow(eta,'+str(i)+');')
        
