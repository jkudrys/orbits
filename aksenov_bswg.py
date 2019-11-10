#!/usr/bin/env python3


from numpy import sqrt, arctan2, arccos


def print_var(*vars):
    for var in vars:
        print(f'{var} = {eval(var)}')


# fndamental
ae = 6378.137  # [km]
# J2 = 1.082626683553151e-3
J2 = 1.082626684e-3
J3 = -2.53265649e-6
GM = 398600.4418  # [km**3/s**2]
# GM = 398600.5  # [km**3/s**2]

# J2,J3 coefficients
C20 = 0.484165371736e-3
C30 = -0.957254173792e-6

tmp = sqrt(J2 - (J3 / (2 * J2)) ** 2)
# c = ae * tmp
# c = 0
c = 209.7290630223342
c2 = c * c

# sigma = J3 / (2 * J2 * tmp)
sigma = -0.035571550267469


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

    # bs+wg mono
    # e2 = 1 - (a2_A ** 2) * (1 - eps_n ** 2 * a3_a2 ** 2 * (4 - 3 * a2_A ** 2) -
    #                         eps_n ** 4 * a3_a2 ** 2 * ((-16 + 32 * a3_A ** 2) +
    #                                                    a2_A ** 2 * (20 - 28 * a3_A ** 2) +
    #                                                    a2_A ** 4 * (-5 + 2 * a3_A ** 2)) -
    #                         eps_n ** 6 * a3_a2 ** 2 * ((64 - 384 * a3_A ** 2 + 448 * a3_A ** 4) +
    #                                                     a2_A ** 2 * (-112 + 544 * a3_A ** 2 - 528 * a3_A ** 4) +
    #                                                     a2_A ** 4 * (56 - 192 * a3_A ** 2 + 136 * a3_A ** 4) +
    #                                                     a2_A ** 6 * (-7 + 9 * a3_A ** 2 - 3 * a3_A ** 4)))

    # bs
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


def orbit_angle(GM, sigma, c, a, e, s, xi0, eta0, r01, z0dot):
    eps = c / (a * (1 - e ** 2))
    ebar = e * (1 + eps ** 2 * (1 - e ** 2) * (1 - 2 * s ** 2)
                # - 4 * eps**3 * sigma * s * (1 - e**2) * (1 - s**2)
                + eps ** 4 * (1 - e ** 2) * ((3 - 16 * s ** 2 + 14 * s ** 4)
                                             - 2 * e ** 2 * (1 - s ** 2) ** 2))
    sigma1 = sqrt(GM * a * (1 - e ** 2)) * (1 + eps ** 2 * (1 - s ** 2) * (3 + e ** 2) / 2
                                            + eps ** 2 * sigma ** 2 * (6 - 7 * s ** 2) / 2
                                            - eps ** 4 * (1 - s ** 2) * ((9 + 11 * s ** 2)
                                                                         + e ** 2 * (6 + 34 * s ** 4)
                                                                         + e ** 4 * (1 + 3 * s ** 2)) / 8)

    sigma2 = sqrt(GM * a * (1 - e ** 2)) * (1 - (eps ** 2 / 2) * (3 - 4 * s ** 2 - e ** 2)
                                            # + 4 * eps**3 * sigma * s * (1 - s**2)
                                            - (eps ** 4 / 8) * ((8 - 72 * s ** 2 + 64 * s ** 4)
                                                                + e ** 2 * (2 - 40 * s ** 2 + 48 * s ** 4) + e ** 4))

    k22 = eps ** 2 * e ** 2 * (s ** 2 - eps ** 2 * (1 - 10 * s ** 2 + 11 * s ** 4 + e ** 2 * s ** 4))
    k12 = eps ** 2 * s ** 2 * (1 + sigma ** 2 - e ** 2 - 4 * eps ** 2 * (1 - s ** 2) * (1 - e ** 2))
    # xi0dot = a*e*sigma2*(1-ebar**2)

    cos_psi0 = (a * (1 - e * ebar) - xi0) / (xi0 * ebar - a * (ebar - e))
    j0 = xi0 ** 2 + c ** 2 * eta0 ** 2
    # sin_psi0 = a*e*(1-ebar**2)*j0*xi0dot

    xi0dot = (xi0 * r01 + c * c * eta0 * z0dot) / (xi0 * xi0 + c * c * eta0 * eta0)

    cos_E = (a - xi0) / (a * e)
    # cos_psi =
    # xi0dot = (xi0*(x0*x0dot + y0*y0dot + (z0 -c*sigma)*z0dot) - c**2*eta0*z0dot)/(xi0**2 + c**2*eta0**2)

    return eps, ebar, sigma1, sigma2, xi0dot, k12, k22, cos_psi0


if __name__ == '__main__':
    x0 = 18693.056970  # [km]
    y0 = -3373.018460
    z0 = 18420.184627
    x0dot = 2.053622091
    y0dot = 2.928143773
    z0dot = -1.526508464

    state_vec = (x0, y0, z0, x0dot, y0dot, z0dot)

    print(30 * '*')
    xi0, eta0, w0, v02, r02, r01 = spheroidal_coords(state_vec, c, sigma)
    print_var('xi0', 'eta0', 'w0', 'v02', 'r02', 'r01')

    alpha1, alpha2, alpha3 = first_integ(GM, v02, r02, r01, xi0, eta0)
    print_var('alpha1', 'alpha2', 'alpha3')

    a, e, s = orbit_axis(GM, c, sigma, alpha1, alpha2, alpha3)
    print_var('a', 'e', 's')

    eps, ebar, sigma1, sigma2, xi0dot, k12, k22, cos_psi0 = orbit_angle(GM, sigma, c, a, e, s, xi0, eta0, r01, z0dot)
    print_var('eps', 'ebar', 'sigma1', 'sigma2', 'xi0dot', 'k12', 'k22', 'arccos(cos_psi0)')
