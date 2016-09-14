import graph
import probabilisticstate as pst
from scipy import stats
from random import random as rd
'''
Probabilistic version of a graph. It uses probabilistic states instead of 
common states. This allows for a method that computes statistical tests of two
states' morphs, create randomly generated sequences and use statistical criteria
to complete the last level of the graph.
'''


class ProbabilisticGraph(graph.Graph):
    def __init__(self, states=[], alphabet=[], path=''):
        p_states = []
        for s in states:
            if s:
                s = pst.ProbabilisticState(s.name, s.outedges)
                p_states.append(s)
        graph.Graph.__init__(self, p_states, alphabet, path)
    
    '''
    Name: compare_morphs
    Input:
        *morphs 1 & 2: morphs of two differents states contained in the current
         graph;
        *alpha: the precision parameter for the statistical test;
        *test: the type of statistical test that will be used. Right now, chi-
         squared and kolmogorv-smirnov are implemented.
    Output: 
        *[p >= alpha, p]: the first element returns True if the test passes for
         the specified alpha and False otherwise. The second one gives the exact
         p-value that returned from the test.
    Description:
        Performs a statistical test between two state's morphs, which are prob
        distributions of their output edges. First it extracts just the probs
        from one of the morphs and then it extracts it for the second morph, 
        making sure they are in the same order (i.e. the probs for the same 
        symbol are in the same position). It performs the test and returns its
        results.
    '''
    @staticmethod
    def compare_morphs(morph1, morph2, alpha, test):
        probs1 = [float(x[1]) for x in morph1]
        probs2 = []
        symbols_2 = [x[0] for x in morph2]
        #Loop to guarantee the probability distributions are in the same order:
        to_del = []
        i = 0
        for a in morph1:
            if a[0] in symbols_2:
                for b in morph2:
                    if b[0] == a[0]:
                        if float(b[1]) == float(a[1]) == 0.0:
                            to_del.append(i)
                        else:
                            probs2.append(float(b[1]))
            else:
                probs2.append(0.0)
            i += 1
        for d in to_del[::-1]:
            del probs1[d]
        if probs1 == probs2:
            return [True, 1.0]
        else: 
            if test == "chi-squared":
                [X, p] = stats.chisquare(probs1, probs2)         
            elif test == "ks":
                [KS, p] = stats.ks_2samp(probs1, probs2)
            return [p >= alpha, p] 
    
    '''
    Name: generate_sequence
    Input:
        *length: The desired length for the generated sequence;
        *ini_state: initial state from which the sequence-generation walk will
         begin.
    Output: 
        *data: Sequence of user-specified length.
    Description:
        Starts from an initial state and uses the randomstep function from the
        probabilistic state to retrieve a randomly generated symbol and a
        destination state for the next iteration. Stores the names of visited
        states, which might be useful in the future, to check if there are 
        states that are not reached during a certain 
    '''
    @staticmethod
    def generate_sequence(length, ini_state):
        data = ''
        s = ini_state
        visited_states = [s.name]
        for i in range(0, length):
            d, s = s.random_step()
            if s:
                if s.name not in visited_states:
                    visited_states.append(s.name)
            data += d
        return data, visited_states
    
    '''
    Name: expand_last_level
    Input:
        *L: The length of the labels in the level to be considered the last;
        *method: which algorithm will be used to expand the last level;
        *test: statistical test that is used in some of the methods (either chi-
         squared or KS);
        *alpha: statistical test precision parameter.
    Output: 
        *A new probabilistic graph with the last level connected by the chosen
         method.
    Description:
        Chooses from one of three methods of last level expansion to reconnect
        the states at the level L. The methods are:
            *D-Markov: ends the last level with D-Markov connections;
            *Old: Tests the destination state with all its possible suffixes and
             substitutes it by the  one that returns the highest p-value;
            *New: Similar to the old method, but checks only if the suffixes 
             pass the statistical test and return the one with the longest label 
    '''            
    def expand_last_level(self, l, method, alpha=0.95, test='chi-squared', synch_words=[]):
        if method == "dmark":
            h = self.dmark_expansion(l)
        else:
            h = self.expansion(l, alpha, test, method, synch_words)
        return h         
    
    '''
    Name: dmark_expansion
    Input:
        *L: The length of the labels in the level to be considered the last.
    Output: 
        *A new probabilistic graph with the last level connected by the chosen
         method.
    Description:
        Implements the D-Markov expansion. An example of a D-Markov expansion
        will suffice to understandhow it works. 
        Suppose a prob graph over binary alphabet. The state 0101 would go to 
        state 01010 when the 0-labeled edge is taken. If the graph is to be 
        ended by the length 4, the state labeled with the last 4 symbols of 
        01010 is the new destination, that is, 1010. 
        If such a state does not exist, the algorithm will drop the first digit
        of the label until it finds a state that exists. If none does, it will 
        connect to the root.
    '''       
    def dmark_expansion(self, l):
        last_level = [x for x in self.states if x.name_length() == l]
        new_last_level = []
        for x in last_level:
            new_outedges = []
            for e in x.outedges:
                for i in range(1, l+1):
                    if i < l:
                        dest = x.name[i:] + e[0]
                        next_state = self.state_named(dest)
                        if next_state or (e[2] == 0.0):
                            new_outedges.append((e[0], next_state, e[2]))
                            break
                    else:
                        new_outedges.append((e[0], self.root(), e[2]))
            x.outedges = new_outedges
            new_last_level.append(x)
        new_states = [x for x in self.states if x.name_length() < l]
        new_states.extend(new_last_level)
        return ProbabilisticGraph(new_states, self.alphabet)
    
    '''
    Name: expansion
    Input:
        *L: The length of the labels in the level to be considered the last;
        *method: which algorithm will be used to expand the last level;
        *test: statistical test that is used in some of the methods (either chi-
         squared or KS);
        *alpha: statistical test precision parameter.
    Output: 
        *A new probabilistic graph with the last level connected by the chosen
         method.
    Description:
        The new and old methods are very similar. They will make the same 
        statistical tests, but use different criteria to choose from the results
        This function applies the general method for both of them and calls
        an appropriate function to choose between the old and new methods' 
        criteria. 
    '''       
    def expansion(self, l, alpha, test, method, s_words=[]):
        last_level = [x for x in self.states if x.name_length() == l]
        new_last_level = []
        new_states = [x for x in self.states if x.name_length() < l]
        for s in last_level:
            new_outedges = []
            for edge in s.outedges:
                a = edge[0]
                next_name = s.name + a
                true_next = self.state_named(next_name)
                if true_next:
                    lgth = len(next_name)
                    results = []
                    for i in range(1, lgth+1):
                        if i < lgth:
                            candidate = self.state_named(next_name[i:])
                        else:
                            candidate = self.root()
                        if candidate:
                            r = self.compare_morphs(true_next.morph(), candidate.morph(),
                                                    alpha, test)
                        else:
                            r = [False, 0.0]
                        results.append([r, candidate])
                    if method == 'old':
                        new_next = self.old_method(results)
                    elif method == 'new':
                        new_next = self.new_method(results, next_name[1:])
                    elif method == 'omega':
                        new_next = self.omega_method(true_next, results, new_states + new_last_level, next_name[1:],
                                                     alpha, test)
                    elif method == 'omega_inverted':
                        synch_n_ext = []
                        if s_words:
                            for word in new_states + last_level:
                                for syn in s_words:
                                    s_len = len(syn)
                                    prefix = word.name[:s_len]
                                    if syn == prefix:
                                        synch_n_ext.append(word)
                                        break
                        new_next = self.omega_inverted_method(true_next, results, synch_n_ext, next_name[1:],
                                                              alpha, test)
                    new_outedge = (a, new_next, edge[2])
                    new_outedges.append(new_outedge)                     
                else:
                    new_outedges.append((a, None, '0.0'))
            s.outedges = new_outedges
            new_last_level.append(s)
        new_states.extend(new_last_level)
        h = ProbabilisticGraph(new_states, self.alphabet)
        return h
    
    '''
    Name: old_method
    Input:
        *results: A list of 2-tuples. The first element of the 2-tuple is the
         result of a statistical test and the second is a state (which is a 
         suffix of the real expanded state) for which the test was taken.
    Output: 
        *The state for which the test result was the highest.
    Description:
        Implements the old method criterion. 
    '''
    @staticmethod
    def old_method(results):
        w = [r[0][1] for r in results]
        arg = w.index(max(w))
        return results[arg][1]
    
    '''
    Name: new_method
    Input:
        *results: A list of 2-tuples. The first element of the 2-tuple is the
         result of a statistical test and the second is a state (which is a 
         suffix of the real expanded state) for which the test was taken.
    Output: 
        *The state with longest label that passed the test.
    Description:
        Implements the new method criterion. If no state passes the test, it
        simply applies the d-markov criterion for this specific state.
    '''
    @staticmethod
    def new_method(results, default_name):
        w = [c[1] for c in results if c[0][0]]
        if w:
            lens = [len(y.name) for y in w]
            arg = lens.index(max(lens))
            return w[arg]
        else:
            return [x[1] for x in results if x[1].name == default_name][0]

    def omega_method(self, exp, results, new_last_level, default_name, alpha, test):
        w = [c[1] for c in results if c[0][0]]
        if w:
            lens = [len(y.name) for y in w]
            arg = lens.index(max(lens))
            return w[arg]
        else:
            for s in new_last_level:
                r = self.compare_morphs(s.morph(), exp.morph(), alpha, test)
                if r[0]:
                    return s
            return [x[1] for x in results if x[1].name == default_name][0]

    def omega_inverted_method(self, exp, results, synch_words, default_name, alpha, test):
        for s in synch_words:
            r = self.compare_morphs(s.morph(), exp.morph(), alpha, test)
            if r[0]:
                return s
        w = [c[1] for c in results if c[0][0]]
        if w:
            lens = [len(y.name) for y in w]
            arg = lens.index(max(lens))
            return w[arg]
        else:
            return [x[1] for x in results if x[1].name == default_name][0]
    
    '''
    Name: open_graph_file
    Input:
        *path: file path where the graph file is saved.
    Description:
        Adapts super's open_graph_file in order to convert all states to prob
        states.
    '''
    def open_graph_file(self, path):
        aux = graph.Graph(path=path)
        states = []
        for s in aux.states:
            states.append(pst.ProbabilisticState(s.name, s.outedges))
        for ns in states:
            newedges = []
            for edge in ns.outedges:
                if edge[1]:
                    destname = edge[1].name
                    newdest = [x for x in states if x.name == destname][0]
                else:
                    newdest = None
                newedge = (edge[0], newdest, edge[2])
                newedges.append(newedge)
            ns.outedges = newedges
        self.states = states
        self.alphabet = aux.alphabet

    def irreducible(self, seq_len):
        d, v = self.generate_sequence(seq_len, self.states[0])
        irreducible_states = [s for s in self.states if s.name in v]
        self.states = irreducible_states

    def remove_unreachable_states(self):
        g = graph.Graph(self.states, self.alphabet)
        h = g.remove_unreachable_states()
        return ProbabilisticGraph(h.states, h.alphabet)
