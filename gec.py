from random import random

def gec(lgth, q, Q, p0, p1):
    seq = []
    r = random()
    state = 0 if r < Q else 1
    for i in range(lgth):
        p = p0 if state == 0 else p1
        pq = Q if state == 0 else q
        r0 = random()
        s = 0 if r0 < p else 1
        seq.append(s)
        r = random()
        state = 0 | state if r < pq else 1 ^ state
    return seq
