#!/usr/bin/env python3

from math import atan2, acos, asin, factorial
from numpy.polynomial import Polynomial
from numpy.polynomial.polynomial import polydiv, polymul
from numpy import real, sqrt, sin, cos, pi, fabs, array
from scipy.integrate import quadrature, quad, fixed_quad
from scipy.special import ellipk, ellipkinc

rod = 180 / pi


def mprint(*arg):  # *arg - dowolna liczba argumentow
    # arg - string zawierajacy nazwe zmiennej
    for var in arg:
        if isinstance(var, str):
            if var in globals():
                print(var, '=', eval(var), end='; ')
            else:
                print(var)
        else:
            print(var)
    print()


def sgn(x):
    s = 0
    if x > 0:
        s = 1
    if x < 0:
        s = -1
    return s


# from numpy import array, append

if __name__ == '__main__':

    # fndamental
    r0 = 6378.137  # [km]
    J2 = 1.082626683553151e-3
    J3 = -2.5326564853322355e-6
    # GM = 398600.4418 #[km**3/s**2]
    GM = 398600.5  # [km**3/s**2]

    # J2,J3 coefficients
    C20 = 0.484165371736e-3
    C30 = -0.957254173792e-6
    m = 0
    n = 2
    k = 1
    J2 = C20 / sqrt(factorial(n + m) / (factorial(n - m) * (2 * n + 1) * k))
    mprint('J2')
    m = 0
    n = 3
    k = 1
    J3 = C30 / sqrt(factorial(n + m) / (factorial(n - m) * (2 * n + 1) * k))
    mprint('J3')

    # epoch
    # x0 = 10419.62404124620 #[km]
    # y0 = 2934.58461181878
    # z0 = 23096.26074220000
    # x0dot = -1.04945573496
    # y0dot = 3.81110422447
    # z0dot = -0.00542259216

    # epoch
    x0 = 18693.056970  # [km]
    y0 = -3373.018460
    z0 = 18420.184627
    x0dot = 2.053622091
    y0dot = 2.928143773
    z0dot = -1.526508464

    a0 = sqrt(x0 ** 2 + y0 ** 2 + z0 ** 2)
    mprint('a0')

    # aksenov 2007 str. -37-
    tmp = sqrt(J2 - (J3 / (2 * J2)) ** 2)
    c = r0 * tmp

    # c = 0
    # c = 209.729063

    mprint('c')
    c2 = c * c
    sigma = J3 / (2 * J2 * tmp)

    # sigma = 0
    # sigma = -0.03557155

    mprint('sigma')

    # aksenov 2007 str. -62-
    r0dash2 = x0 ** 2 + y0 ** 2 + (z0 - c * sigma) ** 2
    mprint('r0dash2')
    V02 = x0dot ** 2 + y0dot ** 2 + z0dot ** 2
    V0 = sqrt(V02)
    mprint('V02', 'V0')
    r01 = x0 * x0dot + y0 * y0dot + (z0 - c * sigma) * z0dot
    mprint('r01')
    xi02 = ((r0dash2 - c2) / 2) * (1 + sqrt(1 + (4 * c2 * (z0 - c * sigma) ** 2) / ((r0dash2 - c2) ** 2)))
    xi0 = sqrt(xi02)
    mprint('xi0')
    eta0 = (z0 - c * sigma) / xi0
    eta02 = eta0 * eta0
    mprint('eta0')
    cosw0 = x0 / sqrt(x0 ** 2 + y0 ** 2)
    sinw0 = y0 / sqrt(x0 ** 2 + y0 ** 2)
    mprint('cosw0')
    mprint('sinw0', sinw0)
    w0 = atan2(sinw0, cosw0)
    mprint('w0')

    # aks 2007 str. -55-
    J0 = xi02 + c2 * eta02
    mprint('J0')

    # aksenov 2007 str. -63-
    xi0dot = (xi0 * r01 - c2 * eta0 * z0dot) / (xi02 + c2 * eta02)
    mprint('xi0dot')
    eta0dot = (z0dot - xi0dot * eta0) / xi0
    mprint('eta0dot')
    w0dot = (x0 * y0dot - y0 * x0dot) / ((xi02 + c2) * (1 - eta02))
    mprint('w0dot')
    alpha1 = (V02 / 2.0) - (GM * (xi0 - c * sigma * eta0)) / (xi02 + c2 * eta02)
    alpha22 = r0dash2 * V02 - r01 ** 2 - c2 * z0dot ** 2 + (2 * GM * xi0 * eta0 * (c2 * eta0 + c * sigma * xi0)) / (
                xi02 + c2 * eta02)
    alpha2 = sqrt(alpha22)
    alpha3 = -x0dot * y0 + y0dot * x0
    mprint('alpha1')
    mprint('alpha2', 'alpha22')
    mprint('alpha3')

    # xi2, xi3, eta1, eta2
    # wielomian Phi(xi)
    cp4 = 2 * alpha1
    cp3 = 2 * GM
    cp2 = 2 * alpha1 * c2 - alpha22
    cp1 = 2 * GM * c2
    cp0 = c2 * (alpha3 ** 2 - alpha22)

    poly = Polynomial([cp0, cp1, cp2, cp3, cp4])
    proots = poly.roots()
    mprint('proots')
    xiroots = real(proots)
    mprint('xiroots')


    def Phi(xi):
        return cp4 * xi ** 4 + cp3 * xi ** 3 + cp2 * xi ** 2 + cp1 * xi + cp0


    def dPhi(xi):
        return 4 * cp4 * xi ** 3 + 3 * cp3 * xi ** 2 + 2 * cp2 * xi + cp1


    xiroot1 = float(xiroots[2])
    mprint('xiroot1')
    xiroot2 = float(xiroots[3])
    mprint('xiroot2')

    eps = 1e-8  # 0.01 mm
    # root polish newton

    i = 0
    m = 1

    while fabs(m) > eps:
        i = i + 1
        print('-- ', i, ' --')
        mprint('xiroot1', Phi(xiroot1), dPhi(xiroot1))

        m = Phi(xiroot1) / dPhi(xiroot1)
        mprint('m')

        xiroot1 = xiroot1 - m

    mprint('xiroot1', Phi(xiroot1))

    # root polish euler
    # i = 0
    # xiroot11 = float(xiroots[3])
    # xiroot12 = float(xiroots[3])+0.01
    # xiroot10 = xiroot11

    # while fabs(xiroot11-xiroot12) > eps:
    #    i = i+1
    #    mprint('-- ',i,' --')
    #    mprint(xiroot11, Phi(xiroot11))
    #    mprint(xiroot12, Phi(xiroot12))

    #    m = (Phi(xiroot12)-Phi(xiroot11))/(xiroot12-xiroot11)
    #    mprint(m)

    #    n = Phi(xiroot11)-m*xiroot11
    #    mprint(n)
    #    n = Phi(xiroot12)-m*xiroot12
    #    mprint(n)

    #    xiroot10 = -n/m
    #    if sgn(Phi(xiroot10)) == sgn(Phi(xiroot11)):
    #        xiroot11 = xiroot10
    #    else:
    #        xiroot12 = xiroot10

    # xiroot2 = xiroot12

    i = 0
    m = 1

    while fabs(m) > eps:
        i += 1
        print('-- ', i, ' --')
        mprint('xiroot2', Phi(xiroot2), dPhi(xiroot2))

        m = Phi(xiroot2) / dPhi(xiroot2)
        mprint('m')

        xiroot2 -= m

    mprint('xiroot2', Phi(xiroot2))

    mprint('xiroot1', 'xiroot2')

    a = (xiroot1 + xiroot2) / 2
    e = (xiroot2 - xiroot1) / (2 * a)
    e2 = e * e

    mprint('a')
    mprint('e')

    # wielomian F(eta)
    cf4 = -2 * alpha1 * c2
    cf3 = 2 * GM * c * sigma
    cf2 = 2 * alpha1 * c2 - alpha22
    cf1 = -2 * GM * c * sigma
    cf0 = -(alpha3 ** 2 - alpha22)

    poly = Polynomial([cf0, cf1, cf2, cf3, cf4])
    proots = poly.roots()
    mprint('proots')
    etaroots = real(proots)
    mprint('etaroots')


    def F(eta):
        return cf4 * eta ** 4 + cf3 * eta ** 3 + cf2 * eta ** 2 + cf1 * eta + cf0


    def dF(eta):
        return 4 * cf4 * eta ** 3 + 3 * cf3 * eta ** 2 + 2 * cf2 * eta + cf1


    if len(etaroots) == 4:
        etaroot1 = float(etaroots[1])
        mprint('etaroot1')
        etaroot2 = float(etaroots[2])
        mprint('etaroot2')

    if len(etaroots) == 2:
        etaroot1 = float(etaroots[0])
        mprint('etaroot1')
        etaroot2 = float(etaroots[1])
        mprint('etaroot2')

    eps = 1e-8  # 0.01 mm

    # root polish newton
    i = 0
    m = 1

    while fabs(m) > eps:
        i += 1
        print('-- ', i, ' --')
        mprint('etaroot1', F(etaroot1), dF(etaroot1))

        m = F(etaroot1) / dF(etaroot1)
        mprint('m')

        etaroot1 -= m

    mprint('etaroot1', F(etaroot1))

    i = 0
    m = 1

    while fabs(m) > eps:
        i += 1
        print('-- ', i, ' --')
        mprint('etaroot2', F(etaroot2), dF(etaroot2))

        m = F(etaroot2) / dF(etaroot2)
        mprint('m')

        etaroot2 -= m

    mprint('etaroot2', F(etaroot2))

    # etaroot2 = etaroot1

    ##root polish euler
    # i = 0
    # etaroot11 = float(etaroots[2])
    # etaroot12 = float(etaroots[2])+0.1
    # etaroot10 = etaroot11

    # while fabs(etaroot11-etaroot12) > eps:
    #    i = i+1
    #    mprint('-- ',i,' --')
    #    mprint(etaroot11, F(etaroot11))
    #    mprint(etaroot12, F(etaroot12))

    #    m = (F(etaroot12)-F(etaroot11))/(etaroot12-etaroot11)
    #    mprint(m)

    #    n = F(etaroot11)-m*etaroot11
    #    mprint(n)
    #    n = F(etaroot12)-m*etaroot12
    #    mprint(n)

    #    etaroot10 = -n/m
    #    if sgn(F(etaroot10)) == sgn(F(etaroot11)):
    #        etaroot11 = etaroot10
    #    else:
    #        etaroot12 = etaroot10

    # etaroot2 = etaroot12

    mprint(etaroot1, etaroot2)

    delta = etaroot2
    deltastar = etaroot1
    mprint('delta')
    mprint('deltastar')

    # sprawdzenie z wzorami aksenov 2007 str. -80-
    eps = c / (a * (1 - e2))
    Q = (1 - 2 * eps * sigma * delta - eps ** 2 * delta ** 2 * (1 - e2)) / (
                1 + 2 * eps ** 2 * delta ** 2 * (1 + e ** 2) + eps ** 4 * delta ** 4 * (1 - e2) ** 2)

    alp1 = -GM * (1 - eps ** 2 * (1 - e2) * (1 - delta ** 2) * Q) / (2 * a)
    mprint('alp1')
    alp22 = GM * a * (1 - e2) * (1 + 2 * eps ** 2 * (1 + e2) * (1 - delta ** 2) * Q + eps ** 4 * (1 - e ** 2) ** 2 * (
                1 - delta ** 2) * Q)
    alp2 = sqrt(alp22)
    mprint('alp2')
    alp32 = GM * a * (1 - e2) * (1 - delta ** 2) * (1 + 2 * eps ** 2 * (1 + e2) + eps ** 4 * (1 - e2) ** 2) * Q
    alp3 = sqrt(alp32)
    mprint('alp3')

    mprint('różnice alpha', alpha1 - alp1, alpha2 - alp2, alpha3 - alp3)

    # aksenov 2007 str. -99-
    p = -GM / (2 * alpha1) - a
    mprint('p')
    q2 = (c2 * (alpha3 ** 2 - alpha22)) / (2 * alpha1 * a * a * (1 - e2)) - p * p
    mprint('q2')
    if q2 < 1e-16:
        q2 = 0
    q = sqrt(q2)
    mprint('q')

    xi2 = xiroot2
    xi1 = xiroot1

    # aksenov 2007 str. -100-
    n1 = sqrt(p * p + q * q - 2 * p * xi2 + xi2 ** 2)
    n11 = sqrt(p * p + q * q - 2 * p * xi1 + xi1 ** 2)
    ebar = (n1 - n11) / (n1 + n11)
    k22 = ((xi2 - xi1) ** 2 - (n1 - n11) ** 2) / (4 * n1 * n11)
    kbar22 = -k22 / (1 - k22)
    # kbar2 = sqrt(kbar22)

    sigmabar2 = sqrt(-2 * alpha1 * n1 * n11 * (1 - k22))
    sigma2 = sqrt(-2 * alpha1 * n1 * n11)

    mprint('n1')
    mprint('n11')
    mprint('ebar')
    mprint('k22')
    mprint('kbar22')
    mprint('sigmabar2')
    mprint('sigma2')
    # aksenov 2007 str. -101-
    eta1 = deltastar
    eta2 = delta

    p1 = GM * c * sigma / (2 * alpha1) - (eta1 + eta2) * c2 / 2
    mprint('p1')
    q12 = -alpha22 / (2 * alpha1) + c2 * (1 - eta2 ** 2 - eta1 * eta2 - eta1 ** 2) + 2 * GM * c * sigma * (
                eta1 + eta2) / (2 * alpha1) + p1 ** 2
    q1 = sqrt(q12)
    mprint('q1')
    q12 = (alpha22 - alpha3 ** 2) / (2 * alpha1 * eta1 * eta2) + p1 ** 2
    q1 = sqrt(q12)
    mprint('q1')

    # aksenov 2007 str. -103-
    m1 = sqrt(q12 - p1 ** 2 + 2 * p1 * eta2 - c2 * eta2 ** 2)
    m11 = sqrt(q12 - p1 ** 2 + 2 * p1 * eta1 - c2 * eta1 ** 2)
    d = (m11 - m1) / (m1 + m11)
    gamma = (m11 * eta2 + m1 * eta1) / (m1 + m11)
    s = (m11 * eta2 - m1 * eta1) / (m1 + m11)
    khat12 = (c2 * (eta2 - eta1) ** 2 + (m1 - m11) ** 2) / (4 * m1 * m11)
    k12 = khat12 / (1 + khat12)
    # k1 = sqrt(k12)
    sigma1 = sqrt(-2 * alpha1 * m1 * m11 * (1 + khat12))

    # aksenov 1977 str. 72
    m1x = sqrt(q12 - (eta2 - p1) ** 2)
    m11x = sqrt(q12 - (eta1 - p1) ** 2)
    khat12x = ((eta2 - eta1) ** 2 + (m1x - m11x) ** 2) / (4 * m1x * m11x)
    k12x = khat12x / (1 + khat12x)
    sigma1x = sqrt(-2 * alpha1 * c2 * m1x * m11x * (1 + khat12x))

    mprint('m1')
    mprint('m11')
    mprint('d')
    mprint('gamma')
    mprint('s')
    mprint('khat12')
    mprint('k12')
    mprint('sigma1')


    # definicja n!!
    def dfact(n):
        df = 1.0
        if n > 0:
            for x in range(n, 0, -2):
                df = df * x
        return df


    # for x in range(0,20):
    #    mprint(x, dfact(x))

    # definicja K(kk)
    # def K(kk):
    #    K1 = 0
    #    for n in range(1,5):
    #        K1 = K1+(dfact(2*n-1)/dfact(2*n))**2*kk*(n)
    #    return pi*(1+K1)/2

    def K(m):
        K1 = 0
        for n in range(1, 12):
            K1 = K1 + (dfact(2 * n - 1) / dfact(2 * n)) ** 2 * m ** (n)
        return pi * (1 + K1) / 2


    nu = sigma1 * K(kbar22) / (sigmabar2 * K(k12)) - 1
    mprint('nu', repr(nu))

    nu = sigma1 * ellipk(kbar22) / (sigmabar2 * ellipk(k12)) - 1
    mprint('nu', repr(nu))

    I22 = ellipk(kbar22)
    # mprint(repr(I22))
    I12 = ellipk(k12)
    # mprint(repr(I12))
    aa = sigma1 * I22
    bb = sigmabar2 * I12
    # mprint(repr(aa),repr(bb))
    # nu = (aa - bb)/(bb)
    # mprint('nu',repr(nu))

    I22 = K(kbar22)
    # mprint(repr(I22))
    I12 = K(k12)
    # mprint(repr(I12))
    aa = sigma1 * I22
    bb = sigmabar2 * I12
    # mprint(repr(aa),repr(bb))
    # nu = (aa - bb)/(bb)
    # mprint('nu',repr(nu))

    #########
    print('#' * 20)
    mprint('K(k12)', repr(K(k12)))
    mprint('ellip(k12)', repr(ellipk(k12)))


    #########

    # definicja Ckn - kombinacje
    def ckn(k, n):
        return factorial(n) / (factorial(k) * factorial(n - k))


    # definicja kappa2(j)
    def kappa2(j):
        kp0 = pi * sigma1 / (2 * sigmabar2 * K(k12) * j)
        kp = 0
        for n in range(j, 12):
            kp = kp + (dfact(2 * n - 1) * ckn(n + j, 2 * n) * kbar22 ** n) / (2 ** (2 * n) * dfact(2 * n))
        return kp0 * kp


    def kappabar2(j):
        kp0 = -pi / (2 * K(k12) * j)
        kp = 0
        for n in range(j, 12):
            kp = kp + (dfact(2 * n - 1) * ckn(n + j, 2 * n) * k12 ** n) / (2 ** (2 * n) * dfact(2 * n))
        return kp0 * kp


    # aksenov 2007 str. -119-
    A = 1 - gamma ** 2
    B = 2 * (s * gamma - d)
    C = d * d - s * s

    # aksenov 2007 str. -121-
    gammabar = 4 / ((sqrt(A - B + C) + sqrt(A + B + C)) ** 2 - 4 * C)
    alpha = 0.5 * sqrt(gammabar) * sgn(alpha3) * (sqrt(A - B + C) + sqrt(A + B + C))
    beta = 0.5 * sqrt(gammabar) * sgn(alpha3) * (sqrt(A - B + C) - sqrt(A + B + C))

    mprint('gammabar')
    mprint('alpha')
    mprint('beta')

    # aksenov 2007 str. -126-
    n0 = sqrt((-2 * alpha1) ** 3) / GM
    estar = -2 * alpha1 * a * e / GM

    mprint('n0')
    mprint('estar')

    # aksenov 2007 str. -64-
    cosE = (a - xi0) / (a * e)
    sinE = sgn(xi0dot) * sqrt(1 - cosE ** 2)
    E = atan2(sinE, cosE)
    mprint('E', E * rod)

    cospsi = (cosE - ebar) / (1 - ebar * cosE)
    sinpsi = (sqrt(1 - ebar ** 2) * sinE) / (1 - ebar * cosE)
    psi = atan2(sinpsi, cospsi)
    mprint('psi', psi * rod)

    cosphitilde = (eta0 - gamma) / (-s + eta0 * d)
    sinphitilde = sgn(eta0dot) * sqrt(1 - cosphitilde ** 2)
    phitilde = atan2(sinphitilde, cosphitilde)
    mprint('phitilde', phitilde * rod + 360)

    xi = a * (1 - e * ebar + (ebar - e) * cospsi) / (1 + ebar * cospsi)
    mprint('xi', 'xi-xi0', xi - xi0)

    eta = (-s * cosphitilde + gamma) / (1 - d * cosphitilde)
    mprint('eta', 'eta-eta0', eta - eta0)

    # aksenov str -124-
    wtilde0 = atan2(sinphitilde, alpha * cosphitilde - beta)
    mprint('wtilde0', wtilde0 * rod + 360)


    # aksenov str - 104 -

    ################# test
    def iF(eta):
        return 1 / sqrt(F(eta))


    def iPhi(xi):
        return 1 / sqrt(Phi(xi))


    def f34(theta, m):
        f0 = 2 * ellipk(m) * theta / pi
        f1 = 0.0
        for n in range(1, 12):
            si = 0.0
            for j in range(1, n + 1):
                si = si + ckn(n + j, 2 * n) * sin(2 * j * theta) / j
            f1 = f1 + (dfact(2 * n - 1) * m ** n) * si / (2 ** (2 * n) * dfact(2 * n))
        return f0 + f1


    tau_c3 = -quad(iF, eta1, eta, epsabs=1e-20)[0]
    mprint('tau_c3 quad', 'tau_c3')

    tau_c3 = (ellipkinc(phitilde - pi / 2, k12) + K(k12)) / sigma1
    mprint('tau_c3 ellip', repr(tau_c3))

    tau_c3 = f34(phitilde, k12) / sigma1
    mprint('tau_c3 akse', repr(tau_c3))

    tau_c4 = quad(iPhi, xi1, xi, limit=200, epsabs=1e-20)[0]
    mprint('tau_c4 quad', 'tau_c4')

    tau_c4 = ellipkinc(psi, k22) / sigma2
    mprint('tau_c4 ellip', repr(tau_c4))

    tau_c4 = f34(psi, kbar22) / sigmabar2
    mprint('tau_c4 akse', repr(tau_c4))

    c3_c4 = tau_c3 - tau_c4
    mprint('c3_c4')

    omegatilde0 = c3_c4 * sigma1 * pi / (2 * K(k12))
    mprint('omegatilde0', repr(omegatilde0 * rod + 360))

    sk = 0
    skbar = 0
    for j in range(1, 12):
        sk = sk + kappa2(j) * sin(2 * j * psi)
        skbar = skbar + kappabar2(j) * sin(2 * j * phitilde)

    omegatilde0 = phitilde - (1 + nu) * psi - sk - skbar
    mprint('omegatilde0', repr(omegatilde0 * rod + 360))

    #### test ####

    # beta2 = alpha2*2*K(k12)*omegatilde0/(pi*sigma1)
    # mprint('beta2', beta2)

    # I1 = alpha2*quad(iPhi, xi1,  xi,  epsabs = 1e-20)[0]
    # I2 = alpha2*quad(  iF, eta1, eta, epsabs = 1e-20)[0]
    # mprint('I1', 'I2')

    # beta2 = -I1+I2
    # mprint('beta2')
    # omegatilde0 = beta2*sigma1*pi/(alpha2*2*K(k12))
    # mprint('omegatilde0', repr(omegatilde0*rod))

    # aksenov -119-
    # dzielenie wielomianow

    polyl = (1, -2 * d, d * d)
    polym = (A, B, C)

    p1 = polydiv(polyl, polym)
    mprint('p1')

    g0 = d * d / C
    d0 = 1 - d * d * A / C
    d1 = -2 * d - d * d * B / C

    mprint('g0', 'd0', 'd1')


    def mkdtau(par, k, j):
        dtau = [1, 0]
        for n in range(1, j):
            # print(n, dfact(2*n-1), dfact(2*n), k**(2*n))
            dtau.append(dfact(2 * n - 1) * k ** (2 * n) * par ** (2 * n) / dfact(2 * n))
            dtau.append(0)
        return dtau


    dtau = mkdtau(1, sqrt(k12), 10)
    mprint('dtau')

    poly2 = polymul(polyl, dtau)
    pol = polydiv(poly2, polym)

    mprint('poly2')
    mprint('pol')

    R = pol[1][1]
    S = pol[1][0]
    mprint('R', 'S')

    mprint('R/S', R / S)
    mprint('beta/alpha', beta / alpha)
    mprint(alpha3 * gammabar * S / (sigma1 * alpha))


    # al = 0.985
    # ka = 0.091

    # dtt1 = 1/sqrt(1-ka*ka*cos(al)**2)
    # dtt2 = mkdtau(cos(al), ka, 18)

    # print(repr(dtt1), repr(sum(dtt2)))
    # print(dtt1-sum(dtt2))

    def cos2n(phi, n):
        c0 = ckn(n, 2 * n) / (2 ** (2 * n))
        c1 = 0
        for j in range(1, n + 1):
            c1 = c1 + ckn(n + j, 2 * n) * cos(2 * j * phi)
        return c0 + c1 / (2 ** (2 * n - 1))


    def cos2n1(phi, n):
        c1 = 0
        for j in range(1, n + 1):
            c1 = c1 + ckn(n + j - 1, 2 * n - 1) * cos((2 * j - 1) * phi)
        return c1 / (2 ** (2 * n - 2))


    def cosn(phi, n):
        if n % 2:
            return cos2n1(phi, (n + 1) // 2)
        else:
            return cos2n(phi, n // 2)


    x = 0.332
    n = 6
    c1 = cos(x) ** (n)
    c2 = cosn(x, n)
    mprint('c1', 'c2')

    mprint('n')
    # mprint('ckn(n, 2*n)',ckn(n, 2*n))
    z = 2 ** (2 * n)


    # mprint('z')

    def polyc2n(n):
        polc = [ckn(n, 2 * n) / (2 ** (2 * n))]
        polc.append(0)
        for j in range(1, n + 1):
            polc.append(ckn(n + j, 2 * n) / (2 ** (2 * n - 1)))
            polc.append(0)
        return polc


    def polyc2n1(n):
        polc = [0]
        # print(polc)
        for j in range(1, n + 1):
            polc.append(ckn(n + j - 1, 2 * n - 1) / (2 ** (2 * n - 2)))
            polc.append(0)
        return polc


    def polyc2(n):
        if n % 2:
            return polyc2n1((n + 1) // 2)
        else:
            return polyc2n(n // 2)


    pc2ncoef = array(polyc2(n))
    mprint('pc2ncoef')
    # mprint(z*pc2ncoef)

    c = 0
    k = 0
    for i in pc2ncoef:
        c = c + i * cos(k * x)
        k = k + 1

    mprint('c')

    for i in range(1, 12):
        # print(k, i)
        pc2ncoef = array(polyc2(i))
        mprint('i', 'pc2ncoef')

    print(pc2ncoef[11])
