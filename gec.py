from random import random

def gec(lgth, pa, pa0, pb0):
    seq = []
    r = random()
    state = 0 if r < pa else 1
    for i in range(lgth):
        p = pa0 if state == 0 else pb0
        r0 = random()
        s = 0 if r0 < p else 1
        seq.append(s)
        r = random()
        state = 0 if r < pa else 1
    return seq
