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
