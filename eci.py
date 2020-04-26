#!/usr/bin/env python3

from numpy import array, dot, pi, sin, cos, hstack

rod = 180 / pi


def jd(y, m, d, t, julian=0):
    if m < 3:
        y = y - 1
        m = m + 12
    a = int(y / 100)
    b = 2 - a + int(a / 4)
    if julian:
        b = 0
    j = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + b - 1524.5 + t / 24.0
    return j


def ERA(y, m, d, utc, ut1_utc=0):
    # utc in hours [0-24)
    # ut1_utc in seconds
    ut1 = utc + ut1_utc / 3600.0
    tu = jd(y, m, d, ut1) - 2451545
    era = 2 * pi * (0.7790572732640 + 1.00273781191135448 * tu)
    era = era % (2 * pi)
    return era


def GMST(y, m, d, utc, ut1_utc=0):
    # utc in hours [0-24)
    # ut1_utc in seconds
    ut1 = utc + ut1_utc / 3600.0
    t = (jd(y, m, d, ut1) - 2451545) / 36525.0
    sg = 24110.54841 + 8640184.812866 * t + 0.093104 * t * t - 0.0000062 * t * t * t
    sg = (sg % 86400.0) / 3600.0 + ut1
    if sg > 24:
        sg = sg - 24
    return (sg * 15.0) / rod


def GMSTE(y, m, d, utc, ut1_utc=0):
    # utc in hours [0-24)
    # ut1_utc in seconds
    ut1 = utc + ut1_utc / 3600.0
    th = ERA(y, m, d, ut1)
    tt = ut1 + (67.184 / 3600.0)  # leap seconds
    t = (jd(y, m, d, tt) - 2451545) / 36525.0
    sg = th + (
            0.014506 + 4612.156534 * t + 1.3915847 * t * t - 0.00000044 * t ** 3 - 0.000029956 * t ** 4 - 0.0000000368 * t ** 5) / (
                 rod * 3600)
    sg = sg % (2 * pi)
    return sg


def ecef2eci(vec, theta):
    # we = 7.2921151467e-5
    we = 7.292115e-5
    R = array([[cos(theta), -sin(theta), 0],
               [sin(theta), cos(theta), 0],
               [0, 0, 1]])

    Rdot = array([[-sin(theta), -cos(theta), 0],
                  [cos(theta), -sin(theta), 0],
                  [0, 0, 0]])

    x_ecef = array(vec[0:3])
    # print(x_ecef)
    x_eci = dot(R, x_ecef)
    # print(x_eci)

    xdot_ecef = array(vec[3:])
    # print(xdot_ecef)
    xdot_eci = dot(R, xdot_ecef) + dot(Rdot, x_ecef) * we
    # print(xdot_eci)

    return list(hstack((x_eci, xdot_eci)))


def eci2ecef(vec, theta):
    # we = 7.2921151467e-5
    we = 7.292115e-5
    R = array([[cos(theta), sin(theta), 0],
               [-sin(theta), cos(theta), 0],
               [0, 0, 1]])

    Rdot = array([[-sin(theta), cos(theta), 0],
                  [-cos(theta), -sin(theta), 0],
                  [0, 0, 0]])

    x_eci = array(vec[0:3])
    # print(x_ecef)
    x_ecef = dot(R, x_eci)
    # print(x_eci)

    xdot_eci = array(vec[3:])
    # print(xdot_ecef)
    xdot_ecef = dot(R, xdot_eci) + dot(Rdot, x_eci) * we
    # print(xdot_eci)

    return list(hstack((x_ecef, xdot_ecef)))
