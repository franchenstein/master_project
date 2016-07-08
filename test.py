#!/usr/bin
import graphgenerator
import probabilisticgraph
g = probabilisticgraph.ProbabilisticGraph(path='graphs/ternary_even_shift/rtp_L8_alpha0.95_omega.yaml')
h = g.expand_last_level(l=2, method='omega_inverted', synch_words=['0'])
print h
