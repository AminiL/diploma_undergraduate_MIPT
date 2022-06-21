import itertools

from sympy import *
from typing import List
from fractions import Fraction


MIN_K_VAL = 23
K = Symbol('k', nonnegative=True, integer=True)


class IndSetCntExpr:
    def __init__(self, power_bases: List[int], factors: List[Fraction]):
        assert len(power_bases) == len(factors)
        self.power_bases = power_bases.copy()
        self.factors = factors.copy()
        self.remove_zero_factors()

    def sympified(self):
        s = sympify(0)
        for p, f in zip(self.power_bases, self.factors):
            s += sympify(f) * (p**K)
        return s

    def remove_zero_factors(self):
        self.power_bases = [self.power_bases[i] for i in range(len(self.factors)) if self.factors[i] != 0]
        self.factors = [self.factors[i] for i in range(len(self.factors)) if self.factors[i] != 0]

    def is_number(self):
        if len(self.power_bases) > 1:
            return False
        if len(self.power_bases) == 1:
            return self.power_bases[0] == 1
        return True

    def substitute(self, val):
        ret = Fraction(0)
        for p, f in zip(self.power_bases, self.factors):
            ret += f * (p**val)
        return ret

    def __iadd__(self, other):
        for op, of in zip(other.power_bases, other.factors):
            if op in self.power_bases:
                i = self.power_bases.index(op)
                self.factors[i] += of
            else:
                self.power_bases.append(op)
                self.factors.append(of)
        assert len(self.power_bases) == len(self.factors)
        self.remove_zero_factors()
        return self

    def __add__(self, other):
        ret = IndSetCntExpr(self.power_bases, self.factors)
        ret += other
        return ret

    def __mul__(self, other):
        if isinstance(other, self.__class__):
            ret = IndSetCntExpr([], [])
            for (sp, sf) in zip(self.power_bases, self.factors):
                for (op, of) in zip(other.power_bases, other.factors):
                    ret += IndSetCntExpr([op * sp], [of * sf])
            return ret
        else:
            ret = IndSetCntExpr(self.power_bases, self.factors)
            other_num = Fraction(other)
            ret.factors = [f * other_num for f in ret.factors]
            ret.remove_zero_factors()
            return ret

    def __imul__(self, other):
        new_val = self * other
        self.power_bases = new_val.power_bases.copy()
        self.factors = new_val.factors.copy()
        return self

    def __sub__(self, other):
        return self + other * (-1)

    def __isub__(self, other):
        new_val = self - other
        self.power_bases = new_val.power_bases.copy()
        self.factors = new_val.factors.copy()
        return self

    def find_index_of_max_abs_power(self):
        m = max(self.power_bases, key=abs)
        return self.power_bases.index(m)

    def sign_on_limit(self):
        if len(self.power_bases) == 0:
            return 0, 0
        tmp = IndSetCntExpr(self.power_bases, self.factors)
        id = tmp.find_index_of_max_abs_power()

        assert len(tmp.factors) > 0
        if tmp.is_number():
            return (1 if tmp.factors[0] > 0 else 0 if tmp.factors[0] == 0 else -1), 0
        val_0 = tmp.substitute(0)
        last_sign = 1 if val_0 > 0 else 0 if val_0 == 0 else -1
        last_k = 0

        for k in itertools.count(start=1):
            sum = Fraction(0, 1)
            for i in range(len(tmp.factors)):
                if i != id:
                    sum += abs(tmp.factors[i])
            if sum < abs(tmp.factors[id]):
                return last_sign, last_k

            val_k = tmp.substitute(0)

            new_sign = 1 if val_k > 0 else 0 if val_k == 0 else -1
            if new_sign != last_sign:
                last_sign = new_sign
                last_k = k

            for i in range(len(tmp.power_bases)):
                tmp.factors[i] *= tmp.power_bases[i]

        raise NotImplementedError


class PreRootVertex:
    k = Symbol('k', nonnegative=True, integer=True)

    def __init__(self, mod: int, delta: int = 0):
        assert 0 <= mod < 7

        kd = PreRootVertex.k + delta
        self.base_size = (PreRootVertex.k + delta) * 7
        self.mod = mod
        self.delta = delta
        self.sz = self.base_size + mod + 1
        if mod == 0:
            incl_power_base = 27
            incl_factor = Fraction(27**delta) if delta >= 0 else Fraction(1, 27**(-delta))
            self.incl_ind_cnt = IndSetCntExpr([incl_power_base], [incl_factor])
            # self.incl_ind_cnt = 27**kd
            excl_power_base = 35
            excl_factor = Fraction(35**delta) if delta > 0 else Fraction(1, 35**(-delta))
            self.excl_ind_cnt = IndSetCntExpr([excl_power_base], [excl_factor])
            # self.excl_ind_cnt = 35**kd
        elif mod == 1:
            incl_power_base = 27
            incl_factor = (Fraction(27 ** delta) if delta >= 0 else Fraction(1, 27 ** (-delta))) * 3
            self.incl_ind_cnt = IndSetCntExpr([incl_power_base], [incl_factor])
            # self.incl_ind_cnt = 27**kd * 3
            excl_power_base = 35
            excl_factor = (Fraction(35 ** (delta - 5)) if delta - 5 >= 0 else Fraction(1, 35 ** (-delta + 5))) * 97**4
            self.excl_ind_cnt = IndSetCntExpr([excl_power_base], [excl_factor])
            # self.excl_ind_cnt = 35**(kd - 5) * 97**4
        elif mod == 2:
            incl_power_base = 27
            incl_factor = (Fraction(27 ** delta) if delta >= 0 else Fraction(1, 27 ** (-delta))) * 3
            self.incl_ind_cnt = IndSetCntExpr([incl_power_base], [incl_factor])
            # self.incl_ind_cnt = 27**kd * 3
            excl_power_base = 35
            excl_factor = (Fraction(35 ** (delta - 1)) if delta - 1 >= 0 else Fraction(1, 35 ** (-delta + 1))) * 97
            self.excl_ind_cnt = IndSetCntExpr([excl_power_base], [excl_factor])
            # self.excl_ind_cnt = 35**(kd - 1) * 97
        elif mod == 3:
            incl_power_base = 27
            incl_factor = (Fraction(27 ** delta) if delta >= 0 else Fraction(1, 27 ** (-delta))) * 9
            self.incl_ind_cnt = IndSetCntExpr([incl_power_base], [incl_factor])
            # self.incl_ind_cnt = 27**kd * 9
            excl_power_base = 35
            excl_factor = (Fraction(35 ** (delta - 6)) if delta - 6 >= 0 else Fraction(1, 35 ** (-delta + 6))) * 97**5
            self.excl_ind_cnt = IndSetCntExpr([excl_power_base], [excl_factor])
            # self.excl_ind_cnt = 35**(kd - 6) * 97**5
        elif mod == 4:
            incl_power_base = 27
            incl_factor = (Fraction(27 ** delta) if delta >= 0 else Fraction(1, 27 ** (-delta))) * 9
            self.incl_ind_cnt = IndSetCntExpr([incl_power_base], [incl_factor])
            # self.incl_ind_cnt = 27**kd * 9
            excl_power_base = 35
            excl_factor = (Fraction(35 ** (delta - 2)) if delta - 2 >= 0 else Fraction(1, 35 ** (-delta + 2))) * 97**2
            self.excl_ind_cnt = IndSetCntExpr([excl_power_base], [excl_factor])
            # self.excl_ind_cnt = 35**(kd - 2) * 97**2
        elif mod == 5:
            incl_power_base = 27
            incl_factor = (Fraction(27 ** delta) if delta >= 0 else Fraction(1, 27 ** (-delta))) * 27
            self.incl_ind_cnt = IndSetCntExpr([incl_power_base], [incl_factor])
            # self.incl_ind_cnt = 27**(kd + 1)
            excl_power_base = 35
            excl_factor = (Fraction(35 ** (delta - 7)) if delta - 7 >= 0 else Fraction(1, 35 ** (-delta + 7))) * 97**6
            self.excl_ind_cnt = IndSetCntExpr([excl_power_base], [excl_factor])
            # self.excl_ind_cnt = 35**(kd - 7) * 97**6
        elif mod == 6:
            incl_power_base = 27
            incl_factor = (Fraction(27 ** delta) if delta >= 0 else Fraction(1, 27 ** (-delta))) * 27
            self.incl_ind_cnt = IndSetCntExpr([incl_power_base], [incl_factor])
            # self.incl_ind_cnt = 27**(kd + 1)
            excl_power_base = 35
            excl_factor = (Fraction(35 ** (delta - 3)) if delta - 3 >= 0 else Fraction(1, 35 ** (-delta + 3))) * 97**3
            self.excl_ind_cnt = IndSetCntExpr([excl_power_base], [excl_factor])
            # self.excl_ind_cnt = 35**(kd - 3) * 97**3
        self.kd = kd


def get_forest_params(f: List[PreRootVertex]):
    szf = 0
    f_ind_cnt = IndSetCntExpr([1], [Fraction(1)])
    f_excl_ind_cnt = IndSetCntExpr([1], [Fraction(1)])
    for v in f:
        szf += v.sz
        f_excl_ind_cnt *= v.excl_ind_cnt
        f_ind_cnt *= (v.excl_ind_cnt + v.incl_ind_cnt)
    return szf, f_excl_ind_cnt, f_ind_cnt


def swap_forest(f: List[PreRootVertex], t: List[PreRootVertex]):
    szf, f_excl_ind_cnt, f_ind_cnt = get_forest_params(f)
    szt, t_excl_ind_cnt, t_ind_cnt = get_forest_params(t)
    assert simplify(szt - szf) == 0
    return f_ind_cnt - t_ind_cnt, t_excl_ind_cnt - f_excl_ind_cnt


def is_swap_edge(from_forest, to_forest):
    F12, F21_excl = swap_forest(from_forest, to_forest)
    if F12.is_number() and F12.substitute(0) == 0 and F21_excl.is_number() and F21_excl.substitute(0) == 0:  # do not add edge which equalize independent sets count
        return False, -1, F12, F21_excl
    max_since = MIN_K_VAL
    sign, since = F12.sign_on_limit()
    max_since = max(max_since, since)
    if sign > 0:
        sign, since = F21_excl.sign_on_limit()
        max_since = max(max_since, since)
        if sign > 0:
            sign, since = (F12 - F21_excl).sign_on_limit()
            max_since = max(max_since, since)
            if sign > 0: # because i_(F0) / i(F0) could equals ONE (if F0 is empty)
                return True, max_since, F12, F21_excl
            else:
                return False, -1, F12, F21_excl
        else:
            return True, max_since, F12, F21_excl
    return False, -1, F12, F21_excl
