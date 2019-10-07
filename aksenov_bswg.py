from numpy import sqrt, arctan2

# fndamental
ae = 6378.137  # [km]
# J2 = 1.082626683553151e-3
J2 = 1.082626684e-3
J3 = -2.53265649e-6
# GM = 398600.4418 #[km**3/s**2]
GM = 398600.5  # [km**3/s**2]

# J2,J3 coefficients
C20 = 0.484165371736e-3
C30 = -0.957254173792e-6

tmp = sqrt(J2 - (J3 / (2 * J2)) ** 2)
c = ae * tmp
c2 = c * c

# c = 0
# c = 209.729063

sigma = J3 / (2 * J2 * tmp)


# sigma = 0
# sigma = -0.03557155

def spheroidal_coords(state_vec_0, c, sigma):
    x0 = state_vec_0[0]
    y0 = state_vec_0[1]
    z0 = state_vec_0[2]
    x0dot = state_vec_0[3]
    y0dot = state_vec_0[4]
    z0dot = state_vec_0[5]

    # zc = z0 - c * sigma
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


def axis_a(GM, c, alpha1, alpha2, alpha3):
    # series from BS + WG page 37
    eps_n = GM * c / alpha2 ** 2
    A = GM / sqrt(-2 * alpha1)

    a = -GM / (2 * alpha1) * (1 - eps_n ** 2 * alpha3 ** 2 / A ** 2 -
                              eps_n ** 4 * (alpha3 ** 2 / A ** 2) * ((-4 + 8 * alpha3 ** 2 / alpha2 ** 2) +
                                                                     alpha2 ** 2 / A ** 2 * (
                                                                             2 - 3 * alpha3 ** 2 / alpha2 ** 2)) -
                              eps_n ** 6 * (alpha3 ** 2 / A ** 2) * (
                                      (16 - 96 * alpha3 ** 2 / alpha2 ** 2 + 112 * alpha3 ** 4 / alpha2 ** 4) +
                                      alpha2 ** 2 / A ** 2 * (
                                              -16 + 80 * alpha3 ** 2 / alpha2 ** 2 - 80 * alpha3 ** 4 / alpha2 ** 4) +
                                      alpha2 ** 4 / A ** 4 * (
                                              3 - 12 * alpha3 ** 2 / alpha2 ** 2 + 10 * alpha3 ** 4 / alpha2 ** 4)))

    e2 = 1 - (alpha2 ** 2 / A ** 2) * (1 - eps_n ** 2 * alpha3 ** 2 / alpha2 ** 2 * (4 - 3 * alpha2 ** 2 / A ** 2) -
                                       eps_n ** 4 * alpha3 ** 2 / alpha2 ** 2 * ((-16 + 32 * alpha3 ** 2 / A ** 2) +
                                                                                 alpha2 ** 2 / A ** 2 * (
                                                                                         20 - 28 * alpha3 ** 2 / A ** 2) +
                                                                                 alpha2 ** 4 / A ** 4 * (
                                                                                         -5 + 2 * alpha3 ** 2 / A ** 2)) -
                                       eps_n ** 6 * ((64 - 384 * alpha3 ** 2 / A ** 2 + 448 * alpha3 ** 4 / A ** 4) +
                                                     alpha2 ** 2 / A ** 2 * (
                                                             -112 + 544 * alpha3 ** 2 / A ** 2 - 528 * alpha3 ** 4 / A ** 4) +
                                                     alpha2 ** 4 / A ** 4 * (
                                                             56 - 192 * alpha3 ** 2 / A ** 2 + 136 * alpha3 ** 4 / A ** 4) +
                                                     alpha2 ** 6 / A ** 6 * (
                                                             -7 + 9 * alpha3 ** 2 / A ** 2 - 3 * alpha3 ** 4 / A ** 4)))

    return a, sqrt(e2)


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
                            eps_n ** 4 * a3_a2 ** 2 * ((-16 + 32 * a3_A ** 2) +
                                                       a2_A ** 2 * (20 - 28 * a3_A ** 2) +
                                                       a2_A ** 4 * (-5 + 2 * a3_A ** 2)) -
                            eps_n ** 6 * ((64 - 384 * a3_A ** 2 + 448 * a3_A ** 4) +
                                          a2_A ** 2 * (-112 + 544 * a3_A ** 2 - 528 * a3_A ** 4) +
                                          a2_A ** 4 * (56 - 192 * a3_A ** 2 + 136 * a3_A ** 4) +
                                          a2_A ** 6 * (-7 + 9 * a3_A ** 2 - 3 * a3_A ** 4)))

    s2 = 1 - a3_a2 ** 2 * (1 + eps_n ** 2 * a2_A ** 2 * (1 - a3_a2 ** 2) +
                           eps_n ** 2 * sigma ** 2 * (6 - 7 * a3_a2 ** 2) +
                           eps_n ** 4 * a2_A ** 4 * (1 - a3_a2 ** 2) * (1 - 2 * a3_a2 ** 2) +
                           2 * eps_n ** 4 * sigma ** 2 * a2_A ** 2 * (9 - 33 * a3_a2 ** 2 + 25 * a3_a2 ** 4) +
                           eps_n ** 6 * a2_A ** 6 * (1 - a3_a2 ** 2) * (1 - 5 * a3_a2 ** 2 + 5 * a3_a2 ** 4))

    return a, sqrt(e2), sqrt(s2)


if __name__ == '__main__':
    x0 = 18693.056970  # [km]
    y0 = -3373.018460
    z0 = 18420.184627
    x0dot = 2.053622091
    y0dot = 2.928143773
    z0dot = -1.526508464

    state_vec = (x0, y0, z0, x0dot, y0dot, z0dot)

    xi0, eta0, w0, v02, r02, r01 = spheroidal_coords(state_vec, c, sigma)
    alpha1, alpha2, alpha3 = first_integ(GM, v02, r02, r01, xi0, eta0)

    print(alpha1, alpha2, alpha3)

    a, e, s = orbit_axis(GM, c, sigma, alpha1, alpha2, alpha3)

    print(a, e, s)
