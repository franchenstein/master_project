import probabilisticGraph as pg

'''
This class creates a probabilistic graph from a rooted tree with probabilities
g. It takes its states with length D, extracts them and connects them in a 
D-Markov fashion, creating a new graph.
'''
class DMarkov(pg.ProbabilisticGraph):
    def __init__(self, g, D):
        h1 = g.expand_last_level(D, 'dmark')
        h2 = h1.remove_unreachable_states()
        states = [x for x in h2.states if x.name_length() == D]
        pg.ProbabilisticGraph.__init__(self, states, g.alphabet)
