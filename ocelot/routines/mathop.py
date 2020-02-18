import warnings

import numpy as np


def ra_to_rb(ra_list):
    """
    convert reciprocal \AA to reciprocal Born

    :param ra_list: a list of floats in 1/\AA
    :return:
    """
    rb_list = []
    for i in ra_list:
        rb = float(i) * np.pi * 2 * (1.0 / 1.8897259885789)  # so vasp outcar kpt cart is 2pi/A ???
        rb_list.append(rb)
    return rb_list


def ev_to_ha(eigen_list):
    """
    convert eV to Hartree

    :param eigen_list:
    :return:
    """
    ha_list = []
    for i in eigen_list:
        ha = float(i) * 0.0367493
        ha_list.append(ha)
    return ha_list


def fd_reci_2ndder(y, x, x0_index, step=1, accuracy='high'):
    """
    https://en.wikipedia.org/wiki/Finite_difference_coefficient
    x should be uniformly sampled over a seg

    :param y:
    :param x: length should be = (order+2)
    :param x0_index: index of x0 at which derivative is calculated
    :param step: step size based on index
    :param accuracy:
    :return: 1/y''
    """
    if accuracy == 'low':
        cfb = [1.0, -2.0, 1.0]
        cc = [1.0, -2.0, 1.0]
    else:
        cfb = [2.0, -5.0, 4.0, -1.0]  # forward and backward
        # cfb = [35.0/12.0, -26.0/3.0, 19.0/2.0, -14.0/3.0, 11.0/12.0]
        cc = [-1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0]
    if len(x) < 2 * len(cfb):
        warnings.warn('W: too few points for finite diff at accuracy = ' + accuracy)
        return np.inf
    summ = 0.0
    if x0_index == 0 or x0_index == 1 or x0_index == 2:  # forward fd
        for k in range(len(cfb)):
            summ += cfb[k] * y[k * step]
    elif x0_index == len(x) - 1 or x0_index == len(x) - 2 or x0_index == len(x) - 3:
        for k in range(len(cfb)):
            summ += cfb[k] * y[-1 + (-1 * k * step)]
    elif (x0_index + 2 * step) in range(len(x)) and (x0_index - 2 * step) in range(len(x)):
        for k in range(len(cc)):
            summ += cc[k] * y[int(round(x0_index + k * step - ((len(cc) - 1) / 2) * step))]
    else:
        warnings.warn('W: center fd failed as this kpt is too close to seg boundary')
        return np.inf
    der2 = summ / ((x[0] - x[step]) ** 2)
    return 1.0 / der2
