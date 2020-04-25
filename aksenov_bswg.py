#!/usr/bin/env python3

from numpy import sqrt, arctan2, arccos, sign, arcsin, sin, cos, arctan, pi, tan, arange


def m_print(*arg):  # *arg - dowolna liczba argumentow
    for var in arg:
        if isinstance(var, str):
            if var in globals():
                print(var, '=', eval(var), end='; ')
            else:
                print(var)
        else:
            print(var)
    print()


# fndamental const
ae = 6378.137  # [km]
J2 = 1.082626684e-3
J3 = -2.53265649e-6
GM = 398600.4418  # [km**3/s**2]

C20 = 0.484165371736e-3
C30 = -0.957254173792e-6

tmp = sqrt(J2 - (J3 / (2 * J2)) ** 2)
c = 209.729063022334
c2 = c * c

sigma = -0.0355715502674694


def spheroidal_coords(state_vec_0, c, sigma):
    x0 = state_vec_0[0]
    y0 = state_vec_0[1]
    z0 = state_vec_0[2]
    x0dot = state_vec_0[3]
    y0dot = state_vec_0[4]
    z0dot = state_vec_0[5]

    v02 = x0dot ** 2 + y0dot ** 2 + z0dot ** 2
    r02 = x0 ** 2 + y0 ** 2 + (z0 - c * sigma) ** 2
    r01 = x0 * x0dot + y0 * y0dot + (z0 - c * sigma) * z0dot

    xi02 = ((r02 - c2) / 2) * (1 + sqrt(1 + (4 * c2 * (z0 - c * sigma) ** 2) / (r02 - c2) ** 2))
    xi0 = sqrt(xi02)
    eta0 = (z0 - c * sigma) / xi0
    w0 = arctan2(y0, x0)

    return xi0, eta0, w0, v02, r02, r01


def first_integ(GM, v02, r02, r01, xi0, eta0):
    xi02 = xi0 * xi0
    eta02 = eta0 * eta0
    alpha1 = (v02 / 2.0) - (GM * (xi0 - c * sigma * eta0)) / (xi02 + c2 * eta02)
    alpha22 = r02 * v02 - r01 ** 2 - c2 * z0dot ** 2 + (2 * GM * xi0 * eta0 * (c2 * eta0 + c * sigma * xi0)) / (
            xi02 + c2 * eta02)
    alpha2 = sqrt(alpha22)
    alpha3 = -x0dot * y0 + y0dot * x0

    return alpha1, alpha2, alpha3


def orbit_axis(GM, c, sigma, alpha1, alpha2, alpha3):
    # series from BS + WG page 37
    eps_n = GM * c / alpha2 ** 2
    A = GM / sqrt(-2 * alpha1)

    a2_A = alpha2 / A
    a3_A = alpha3 / A
    a3_a2 = alpha3 / alpha2

    a = -GM / (2 * alpha1) * (1 - eps_n ** 2 * a3_A ** 2 -
                              eps_n ** 4 * (a3_A ** 2) * ((-4 + 8 * a3_a2 ** 2) +
                                                          a2_A ** 2 * (2 - 3 * a3_a2 ** 2)) -
                              eps_n ** 6 * (a3_A ** 2) * ((16 - 96 * a3_a2 ** 2 + 112 * a3_a2 ** 4) +
                                                          a2_A ** 2 * (-16 + 80 * a3_a2 ** 2 - 80 * a3_a2 ** 4) +
                                                          a2_A ** 4 * (3 - 12 * a3_a2 ** 2 + 10 * a3_a2 ** 4)))

    e2 = 1 - (a2_A ** 2) * (1 - eps_n ** 2 * a3_a2 ** 2 * (4 - 3 * a2_A ** 2) -
                            eps_n ** 4 * a3_a2 ** 2 * ((-16 + 32 * a3_a2 ** 2) +
                                                       a2_A ** 2 * (20 - 28 * a3_a2 ** 2) +
                                                       a2_A ** 4 * (-5 + 2 * a3_a2 ** 2)) -
                            eps_n ** 6 * a3_a2 ** 2 * ((64 - 384 * a3_a2 ** 2 + 448 * a3_a2 ** 4) +
                                                       a2_A ** 2 * (-112 + 544 * a3_a2 ** 2 - 528 * a3_a2 ** 4) +
                                                       a2_A ** 4 * (56 - 192 * a3_a2 ** 2 + 136 * a3_a2 ** 4) +
                                                       a2_A ** 6 * (-7 + 9 * a3_a2 ** 2 - 3 * a3_a2 ** 4)))

    s2 = 1 - a3_a2 ** 2 * (1 + eps_n ** 2 * a2_A ** 2 * (1 - a3_a2 ** 2) +
                           eps_n ** 2 * sigma ** 2 * (6 - 7 * a3_a2 ** 2) +
                           eps_n ** 4 * a2_A ** 4 * (1 - a3_a2 ** 2) * (1 - 2 * a3_a2 ** 2) +
                           2 * eps_n ** 4 * sigma ** 2 * a2_A ** 2 * (9 - 33 * a3_a2 ** 2 + 25 * a3_a2 ** 4) +
                           eps_n ** 6 * a2_A ** 6 * (1 - a3_a2 ** 2) * (1 - 5 * a3_a2 ** 2 + 5 * a3_a2 ** 4))

    return a, sqrt(e2), sqrt(s2)


def orbit_axis2(GM, c, sigma, alpha1, alpha2, alpha3):
    eps_n = GM * c / alpha2 ** 2
    eps_n2 = eps_n ** 2
    eps_n4 = eps_n ** 4
    eps_n6 = eps_n ** 6
    A = GM / sqrt(-2 * alpha1)

    a2_A = alpha2 / A
    a2_A2 = a2_A ** 2
    a2_A4 = a2_A ** 4
    a2_A6 = a2_A ** 6

    a3_A2 = (alpha3 / A) ** 2

    a3_a2 = alpha3 / alpha2
    a3_a22 = a3_a2 ** 2
    a3_a24 = a3_a2 ** 4

    a = -GM / (2 * alpha1) * (1 - eps_n2 * a3_A2 -
                              eps_n4 * a3_A2 * ((-4 + 8 * a3_a22) +
                                                a2_A2 * (2 - 3 * a3_a22)) -
                              eps_n6 * a3_A2 * ((16 - 96 * a3_a22 + 112 * a3_a24) +
                                                a2_A2 * (-16 + 80 * a3_a22 - 80 * a3_a24) +
                                                a2_A4 * (3 - 12 * a3_a22 + 10 * a3_a24)))

    e2 = 1 - a2_A2 * (1 - eps_n2 * a3_a22 * (4 - 3 * a2_A2) -
                      eps_n4 * a3_a22 * ((-16 + 32 * a3_a22) +
                                         a2_A2 * (20 - 28 * a3_a22) +
                                         a2_A4 * (-5 + 2 * a3_a22)) -
                      eps_n6 * a3_a22 * ((64 - 384 * a3_a22 + 448 * a3_a24) +
                                         a2_A2 * (-112 + 544 * a3_a22 - 528 * a3_a24) +
                                         a2_A4 * (56 - 192 * a3_a22 + 136 * a3_a24) +
                                         a2_A6 * (-7 + 9 * a3_a22 - 3 * a3_a24)))

    s2 = 1 - a3_a22 * (1 + eps_n2 * a2_A2 * (1 - a3_a22) +
                       eps_n2 * sigma ** 2 * (6 - 7 * a3_a22) +
                       eps_n4 * a2_A4 * (1 - a3_a22) * (1 - 2 * a3_a22) +
                       2 * eps_n4 * sigma ** 2 * a2_A2 * (9 - 33 * a3_a22 + 25 * a3_a24) +
                       eps_n6 * a2_A6 * (1 - a3_a22) * (1 - 5 * a3_a22 + 5 * a3_a24))

    return a, sqrt(e2), sqrt(s2)


def orbit_angle(GM, sigma, c, a, e, s, xi0, eta0, r01, z0dot, w0):
    eps = c / (a * (1 - e ** 2))
    ebar = e * (1 + eps ** 2 * (1 - e ** 2) * (1 - 2 * s ** 2)
                + eps ** 4 * (1 - e ** 2) * ((3 - 16 * s ** 2 + 14 * s ** 4)
                                             - 2 * e ** 2 * (1 - s ** 2) ** 2))

    sigma1 = sqrt(GM * a * (1 - e ** 2)) * (1 + (eps ** 2 / 2) * (1 - s ** 2) * (3 + e ** 2)
                                            + (eps ** 2 * sigma ** 2 / 2) * (6 - 7 * s ** 2)
                                            - (eps ** 4 / 8) * (1 - s ** 2)
                                            * ((9 + 11 * s ** 2) + e ** 2 * (6 + 34 * s ** 4) + e ** 4
                                               * (1 + 3 * s ** 2)))

    sigma2 = sqrt(GM * a * (1 - e ** 2)) * (1 - (eps ** 2 / 2) * (3 - 4 * s ** 2 - e ** 2)
                                            # + 4 * eps**3 * sigma * s * (1 - s**2)
                                            - (eps ** 4 / 8) * ((8 - 72 * s ** 2 + 64 * s ** 4)
                                                                + e ** 2 * (2 - 40 * s ** 2 + 48 * s ** 4) + e ** 4))

    k12 = eps ** 2 * s ** 2 * (1 + sigma ** 2 - e ** 2 - 4 * eps ** 2 * (1 - s ** 2) * (1 - e ** 2))
    k22 = eps ** 2 * e ** 2 * (s ** 2 - eps ** 2 * (1 - 10 * s ** 2 + 11 * s ** 4 + e ** 2 * s ** 4))

    cos_psi0 = (a * (1 - e * ebar) - xi0) / (xi0 * ebar - a * (ebar - e))

    xi0dot = (xi0 * r01 + c * c * eta0 * z0dot) / (xi0 * xi0 + c * c * eta0 * eta0)

    sin_psi0 = (a * e * (1 - ebar ** 2) * (xi0 ** 2 + c ** 2 * eta0 ** 2) * xi0dot) \
               / (sigma2 * (xi0 * ebar - a * (ebar - e)) ** 2 * sqrt(1 - k22 * (1 - cos_psi0 ** 2)))

    cos_E = (a - xi0) / (a * e)
    cos_psi1 = (cos_E - ebar) / (1 - ebar * cos_E)
    # print('cos_psi1', cos_psi1)
    psi0 = arctan2(sin_psi0, cos_psi0)

    d = eps * sigma * s * (1 - eps ** 2 * ((5 - 6 * s ** 2) - e ** 2 * (1 - 2 * s ** 2)))
    gamma = -eps * sigma * (1 - 2 * s ** 2 - eps ** 2 * ((3 - 12 * s ** 2 + 10 * s ** 4) + e ** 2 * (1 - 2 * s ** 4)))

    eta0dot = (z0dot - eta0 * xi0dot) / xi0
    sin_theta0 = (eta0 - gamma) / (s - eta0 * d)
    cos_theta0 = ((s - gamma * d) * (xi0 ** 2 + c ** 2 * eta0 ** 2) * eta0dot) \
                 / (sigma1 * (s - eta0 * d) ** 2 * sqrt(1 - k12 * sin_theta0 ** 2))

    theta0 = arctan2(sin_theta0, cos_theta0)

    # print(d, gamma, eta0dot, sin_theta0, cos_theta0, theta0)

    nu = (eps ** 2 / 4) * (1 + sigma ** 2) * (12 - 15 * s ** 2) + (eps ** 4 / 64) * (288 - 1296 * s ** 2 + 1035 * s ** 4
                                                                                     - e ** 2 * (
                                                                                             144 + 288 * s ** 2 - 510 * s ** 4))
    omega0 = theta0 - (1 + nu) * psi0 - (k12 / 8) * (1 + k12 / 2) * sin(2 * theta0) + (k22 / 8) * (
            1 + nu + k22 / 2) * sin(2 * psi0)
    omega1 = nu * psi0 + omega0

    cos_i = sqrt(1 - s ** 2)
    # print(nu, omega0, omega1, cos_i)

    mu = -(3 * cos_i / 2) * (eps ** 2 * (1 + sigma ** 2) + (eps ** 4 / 8) * (6 - 17 * s ** 2 - 24 * e ** 2 * s ** 2))
    mu1 = -2 * eps ** 2 * e * cos_i * (1 + (eps ** 2 / 8) * ((4 - 28 * s ** 2) - e ** 2 * (6 + 7 * s ** 2)))
    mu2 = (-eps ** 2 * e ** 2 * cos_i / 4) * (1 - (eps ** 2 / 4) * ((22 + s ** 2) + e ** 2 * (2 + s ** 2)))
    mu3 = eps ** 4 * e ** 3 * cos_i * (2 - s ** 2) / 4
    mu4 = eps ** 4 * e ** 4 * cos_i * (2 - s ** 2) / 64
    mu11 = eps ** 3 * sigma * cos_i * s * (1 - e ** 2)
    mu21 = eps ** 4 * s ** 2 * cos_i * (1 - e ** 2) ** 2 / 32
    # print('mu', mu, mu1, mu2, mu3, mu4, mu11, mu21)

    beta = 2 * eps * sigma * s * cos_i * (1 - eps ** 2 * (4 - 5 * s ** 2 + e ** 2 * s ** 2))
    # print('beta', beta)

    # print('>>>')
    # print(cos_i, theta0, sin(theta0), beta, cos(theta0))

    omegahat = arctan2((cos_i * sin(theta0) + beta), cos(theta0))
    if omegahat < 0:
        omegahat += 2 * pi
    # print(w0 - omegahat)

    OMEGA0 = w0 - omegahat - mu * psi0 - mu1 * sin(psi0) - mu2 * sin(2 * psi0) - mu3 * sin(3 * psi0) \
             - mu4 * sin(4 * psi0) - mu11 * cos(psi0 + omega1) - mu21 * sin(2 * (psi0 + omega1))
    # print('OMEGA0', OMEGA0)

    E0 = 2 * arctan(sqrt((1 - ebar) / (1 + ebar)) * tan(psi0 / 2))
    # print(E0)

    lm = -(3 / 16) * eps ** 4 * (1 - e ** 2) ** (3 / 2) * (8 - 32 * s ** 2 + 25 * s ** 4)
    lm1 = -(1 / 4) * eps ** 4 * s ** 2 * e * (4 - 5 * s ** 2) * (1 - e ** 2) ** (3 / 2)
    lm2 = (3 / 32) * eps ** 4 * s ** 4 * e ** 2 * (1 - e ** 2) ** (3 / 2)
    lm11 = 0.5 * eps ** 3 * sigma * s * (4 - 5 * s ** 2) * (1 - e ** 2) ** (3 / 2)
    lm21 = -(eps ** 2 / 4) * s ** 2 * (1 - e ** 2) ** (3 / 2) * (
            1 - (eps ** 2 / 4) * ((12 - 13 * s ** 2) - e ** 2 * (4 - 5 * s ** 2)))
    lm31 = -(eps ** 3 / 6) * sigma * s ** 3 * (1 - e ** 2) ** (3 / 2)
    lm41 = -(eps ** 4 / 64) * s ** 4 * (1 - e ** 2) ** (5 / 2)
    lm221 = (eps ** 4 / 16) * s ** 4 * e ** 2 * (1 - e ** 2) ** (3 / 2)

    # print('lm:', lm, lm1, lm2, lm11)
    # print('lm:', lm21, lm31, lm41, lm221)

    estar = e * (1 - eps ** 2 * (1 - e ** 2) * (1 - s ** 2) + eps ** 4 * s ** 2 * (1 - e ** 2) * (3 + e ** 2))
    # print('estar', estar)

    M0 = E0 - estar * sin(E0) - lm * psi0 + lm1 * sin(psi0) + lm2 * sin(2 * psi0) + lm11 * cos(psi0 + omega1) \
         + lm21 * sin(2 * (psi0 + omega1)) + lm31 * cos(3 * (psi0 + omega1)) + lm41 * sin(4 * (psi0 + omega1)) \
         + lm221 * sin(2 * psi0) * cos(2 * (psi0 + omega1))

    return eps, ebar, sigma1, sigma2, xi0dot, psi0, theta0, M0, OMEGA0, omega0


def s_vec(dt, alpha1, alpha3, a, e, s, OMEGA0, M0, omega0, c, sigma, GM):
    eps = c / (a * (1 - e ** 2))
    n0 = (-2 * alpha1) ** (3 / 2) / GM
    estar = e * (1 - eps ** 2 * (1 - e ** 2) * (1 - s ** 2) + eps ** 4 * s ** 2 * (1 - e ** 2) * (3 + e ** 2))
    ebar = e * (1 + eps ** 2 * (1 - e ** 2) * (1 - 2 * s ** 2)
                + eps ** 4 * (1 - e ** 2) * ((3 - 16 * s ** 2 + 14 * s ** 4)
                                             - 2 * e ** 2 * (1 - s ** 2) ** 2))
    nu = (eps ** 2 / 4) * (1 + sigma ** 2) * (12 - 15 * s ** 2) + (eps ** 4 / 64) * (288 - 1296 * s ** 2 + 1035 * s ** 4
                                                                                     - e ** 2 * (
                                                                                             144 + 288 * s ** 2 - 510 * s ** 4))

    # print('n0', n0)

    lm = -(3 / 16) * eps ** 4 * (1 - e ** 2) ** (3 / 2) * (8 - 32 * s ** 2 + 25 * s ** 4)
    lm1 = -(1 / 4) * eps ** 4 * s ** 2 * e * (4 - 5 * s ** 2) * (1 - e ** 2) ** (3 / 2)
    lm2 = (3 / 32) * eps ** 4 * s ** 4 * e ** 2 * (1 - e ** 2) ** (3 / 2)
    lm11 = 0.5 * eps ** 3 * sigma * s * (4 - 5 * s ** 2) * (1 - e ** 2) ** (3 / 2)
    lm21 = -(eps ** 2 / 4) * s ** 2 * (1 - e ** 2) ** (3 / 2) * (
            1 - (eps ** 2 / 4) * ((12 - 13 * s ** 2) - e ** 2 * (4 - 5 * s ** 2)))
    lm31 = -(eps ** 3 / 6) * sigma * s ** 3 * (1 - e ** 2) ** (3 / 2)
    lm41 = -(eps ** 4 / 64) * s ** 4 * (1 - e ** 2) ** (5 / 2)
    lm221 = (eps ** 4 / 16) * s ** 4 * e ** 2 * (1 - e ** 2) ** (3 / 2)

    M = n0 * (dt) + M0

    E = 0
    E0 = M
    psi_i = 0
    omega_i = omega0

    while abs(E - E0) > 1e-16:
        E = E0
        psi_i = 2 * arctan(sqrt((1 + ebar) / (1 - ebar)) * tan(E / 2))

        omega_i = nu * psi_i + omega0

        E0 = M + estar * sin(E) + lm * psi_i - lm1 * sin(psi_i) - lm2 * sin(2 * psi_i) - lm11 * cos(psi_i + omega_i) \
             - lm21 * sin(2 * (psi_i + omega_i)) - lm31 * cos(3 * (psi_i + omega_i)) - lm41 * sin(4 * (psi_i + omega_i)) \
             - lm221 * sin(2 * psi_i) * cos(2 * (psi_i + omega_i))
        # print('iter E', E, E0, E - E0, psi_i, omega_i)

    E = E0
    psi = psi_i
    omega = omega_i

    cos_i = sqrt(1 - s ** 2)
    mu = -(3 * cos_i / 2) * (eps ** 2 * (1 + sigma ** 2) + (eps ** 4 / 8) * (6 - 17 * s ** 2 - 24 * e ** 2 * s ** 2))
    mu1 = -2 * eps ** 2 * e * cos_i * (1 + (eps ** 2 / 8) * ((4 - 28 * s ** 2) - e ** 2 * (6 + 7 * s ** 2)))
    mu2 = (-eps ** 2 * e ** 2 * cos_i / 4) * (1 - (eps ** 2 / 4) * ((22 + s ** 2) + e ** 2 * (2 + s ** 2)))
    mu3 = eps ** 4 * e ** 3 * cos_i * (2 - s ** 2) / 4
    mu4 = eps ** 4 * e ** 4 * cos_i * (2 - s ** 2) / 64
    mu11 = eps ** 3 * sigma * cos_i * s * (1 - e ** 2)
    mu21 = eps ** 4 * s ** 2 * cos_i * (1 - e ** 2) ** 2 / 32

    OMEGA = mu * psi + OMEGA0 + mu1 * sin(psi) + mu2 * sin(2 * psi) + mu3 * sin(3 * psi) + mu4 * sin(4 * psi) \
            + mu11 * cos(psi + omega) + mu21 * sin(2 * (psi + omega))

    # print('E', E)
    # print('psi', psi)
    # print('omega', omega)
    # print('OMEGA', OMEGA)

    k12 = eps ** 2 * s ** 2 * (1 + sigma ** 2 - e ** 2 - 4 * eps ** 2 * (1 - s ** 2) * (1 - e ** 2))
    k22 = eps ** 2 * e ** 2 * (s ** 2 - eps ** 2 * (1 - 10 * s ** 2 + 11 * s ** 4 + e ** 2 * s ** 4))

    theta = psi + omega - (k22 / 8) * (1 + nu + k22 / 2) * sin(2 * psi) + (k12 / 8) * (1 + k12 / 2) * sin(
        2 * (psi + omega)) \
            + (3 / 256) * k22 ** 2 * sin(4 * psi) + (k12 ** 2 / 256) * sin(4 * (psi + omega)) - (k12 * k22 / 32) * sin(
        2 * psi) * cos(2 * (psi * omega))

    # print('theta', theta)
    ksi = a * (1 - e * cos(E))
    # print('ksi', ksi)

    d = eps * sigma * s * (1 - eps ** 2 * ((5 - 6 * s ** 2) - e ** 2 * (1 - 2 * s ** 2)))

    ro = sqrt((1 - eps ** 2 * sigma ** 2) * (ksi ** 2 + c ** 2)) / (1 + d * sin(theta))
    ro1 = ksi / (1 + d * sin(theta))
    beta = 2 * eps * sigma * s * cos_i * (1 - eps ** 2 * (4 - 5 * s ** 2 + e ** 2 * s ** 2))
    gamma = -eps * sigma * (1 - 2 * s ** 2 - eps ** 2 * ((3 - 12 * s ** 2 + 10 * s ** 4) + e ** 2 * (1 - 2 * s ** 4)))

    x = ro * (cos(theta) * cos(OMEGA) - cos_i * sin(theta) * sin(OMEGA) - beta * sin(OMEGA))
    y = ro * (cos(theta) * sin(OMEGA) + cos_i * sin(theta) * cos(OMEGA) + beta * cos(OMEGA))
    z = c * sigma + ro1 * (s * sin(theta) + gamma)

    w = OMEGA + arctan2(cos_i * sin(theta) + beta, cos(theta))
    eta = (s * sin(theta) + gamma) / (1 + d * sin(theta))

    sigma1 = sqrt(GM * a * (1 - e ** 2)) * (1 + (eps ** 2 / 2) * (1 - s ** 2) * (3 + e ** 2)
                                            + (eps ** 2 * sigma ** 2 / 2) * (6 - 7 * s ** 2)
                                            - (eps ** 4 / 8) * (1 - s ** 2)
                                            * ((9 + 11 * s ** 2) + e ** 2 * (6 + 34 * s ** 4) + e ** 4
                                               * (1 + 3 * s ** 2)))

    sigma2 = sqrt(GM * a * (1 - e ** 2)) * (1 - (eps ** 2 / 2) * (3 - 4 * s ** 2 - e ** 2)
                                            # + 4 * eps**3 * sigma * s * (1 - s**2)
                                            - (eps ** 4 / 8) * ((8 - 72 * s ** 2 + 64 * s ** 4)
                                                                + e ** 2 * (2 - 40 * s ** 2 + 48 * s ** 4) + e ** 4))

    ksidot = a * e * sigma2 * (1 - ebar ** 2) * sin(psi) * sqrt(1 - k22 * sin(psi) ** 2) \
             / ((ksi ** 2 + c ** 2 * eta ** 2) * (1 + ebar * cos(psi)) ** 2)

    etadot = (s - gamma * d) * sigma1 * cos(theta) * sqrt(1 - k12 * sin(theta) ** 2) \
             / ((ksi ** 2 + c ** 2 * eta ** 2) * (1 + d * sin(theta)) ** 2)

    wdot = alpha3 / ((ksi ** 2 + c ** 2) * (1 - eta ** 2))

    # print('dot:', ksidot, etadot, wdot)

    xdot = (x * ksi * ksidot) / (ksi ** 2 + c ** 2) - (x * eta * etadot) / (1 - eta ** 2) - y * wdot
    ydot = (y * ksi * ksidot) / (ksi ** 2 + c ** 2) - (y * eta * etadot) / (1 - eta ** 2) + x * wdot
    zdot = eta * ksidot + ksi * etadot

    return x, y, z, xdot, ydot, zdot


def orb_params(state_vector):
    x0, y0, z0, x0dot, y0dot, z0dot = state_vector
    xi0, eta0, w0, v02, r02, r01 = spheroidal_coords(state_vector, c, sigma)
    alpha1, alpha2, alpha3 = first_integ(GM, v02, r02, r01, xi0, eta0)
    a, e, s = orbit_axis2(GM, c, sigma, alpha1, alpha2, alpha3)
    eps, ebar, sigma1, sigma2, xi0dot, psi0, theta0, M0, OMEGA0, omega0 = \
        orbit_angle(GM, sigma, c, a, e, s, xi0, eta0, r01, z0dot, w0)

    return (alpha1, alpha3, a, e, s, OMEGA0, M0, omega0)


def propagate(dt, orb_params):
    alpha1, alpha3, a, e, s, OMEGA0, M0, omega0 = orb_params
    sv = s_vec(dt, alpha1, alpha3, a, e, s, OMEGA0, M0, omega0, c, sigma, GM)

    return sv


if __name__ == '__main__':

    x0 = 3.98131471724008E+0003
    y0 = -1.37962749950665E+0004
    z0 = 2.10809814453000E+0004
    x0dot = 3.51016302791983E+0000
    y0dot = -1.14609191177291E+0000
    z0dot = -1.41139984131000E+0000

    op = orb_params([x0, y0, z0, x0dot, y0dot, z0dot])

    for dt in arange(0, 1801, 60):
        [x, y, z, xdot, ydot, zdot] = propagate(dt, op)
        print(dt, [x, y, z, xdot, ydot, zdot])

    # print('diff [m]', 30 * '-')
    #
    # print('x-x0', x - x0)
    # print('y-y0', y - y0)
    # print('z-z0', z - z0)
    # print('xdot-x0dot', xdot - x0dot)
    # print('ydot-y0dot', ydot - y0dot)
    # print('zdot-z0dot', zdot - z0dot)
