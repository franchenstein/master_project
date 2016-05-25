#!/usr/bin
import graphgenerator
gg = graphgenerator.GraphGenerator('graphs/even_shift/rtp_L4_dmark.json', ['0'], 'graphs/even_shift/mk1_dmark.json')
gg.mk1('chi-squared', 0.95)
