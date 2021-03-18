#pragma once

#include <cstddef>

static size_t const negamma_knots = 24;
static double const egamma_knots[negamma_knots] = {
    0.015,  0.025,  0.0375, 0.0525, 0.07,  0.09,  0.1125, 0.1375,
    0.1625, 0.1875, 0.225,  0.275,  0.325, 0.375, 0.45,   0.55,
    0.7,    0.9,    1.15,   1.45,   1.8,   2.25,  2.75,   3.5};

static size_t const nenu_knots = 37;
static double const enu_knots[nenu_knots] = {
    0.15, 0.175, 0.2,  0.225, 0.25, 0.275, 0.3,  0.325, 0.35, 0.375,
    0.4,  0.425, 0.45, 0.55,  0.65, 0.75,  0.85, 0.95,  1.15, 1.35,
    1.55, 1.75,  1.95, 2.1,   2.5,  2.9,   3.3,  3.7,   4.1,  4.75,
    5.5,  6.25,  7.,   7.75,  8.5,  9.25,  10.};