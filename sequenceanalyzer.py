import yaml
import numpy as num

'''
This class opens a sequence saved in a yaml file and has methods to compute
various parameters based on said sequence.
'''


class SequenceAnalyzer():
    
    def __init__(self, path, probabilities=[], alphabet=[],
                 conditional_probabilities=[], conditional_entropy=[],
                 kldivergence=0, l1metric=0, autocorrelation=[]):
        self.probabilities = probabilities
        self.alphabet = alphabet
        self.conditional_probabilities = conditional_probabilities
        self.conditional_entropy = conditional_entropy
        self.kldivergence = kldivergence
        self.l1metric = l1metric
        self.autocorrelation = autocorrelation
        self.sequence_path = path
        with open(path, 'r') as file_:
            print "Sequence Analyzer opening sequence at: " + path
            self.seq = yaml.load(file_)
            print "*****************"
            print "Sequence loaded!"
            print "*****************"
    
    '''
    Name: calc_probs
    Input:
        *L: maximum sub-sequence length to be analyzed.
    Output: 
        *probabilities: a list of dictionaries. Each dictionary contains keys
         that are sequences of the same length. The value associated to a key
         is a probability of that subsequence appearing in the original sequence
        *alphebt: the unique symbols that appear in the sequence.
    Description:
        Checks the number of occurances of subsequences of lengths from 1 to L.
        Divides the number of occurances by the sequence's length in order to
        obtain relative frequencies. Creates a dictionary for subsequences of 
        each length. When checking for subsequences of length 1, the method 
        records each individual symbol that appears and stores it as the
        sequence's alphabet.
    '''             
    def calc_probs(self, L):
        self.probabilities = []
        self.alphabet = []
        print "Calculating subsequence probabilities for sequence at: " + self.sequence_path
        print "L = " + str(L)
        for l in range(1, L + 1):
            print "Sequence: " + self.sequence_path
            print "Calculating probabilities of subsequences of length: " + str(l)
            current_probs = {}
            for i in range(0, len(self.seq) - (l - 1)):
                current_value = ''.join(str(e) for e in self.seq[i:i+l])
                #Deduce alphabet:
                if l == 1:
                    if not (current_value in self.alphabet):
                        self.alphabet.append(current_value)
                if not current_value in current_probs.keys():
                    current_probs[current_value] = 1
                else:
                    current_probs[current_value] += 1
            for key in current_probs.keys():
                current_probs[key] /= float(len(self.seq))
            if not sum(current_probs.values()) == 1.0:
                dif = (1.0 - sum(current_probs.values()))/len(current_probs.values())
                for k in current_probs.keys():
                    current_probs[k] += dif
            self.probabilities.append(current_probs)
        print "*****************"
        print "Sequence: " + self.sequence_path
        print "Probabilities calculated!"
        print "*****************"
        return [self.probabilities, self.alphabet] 
    
    '''
    Name: calc_cond_probs
    Input:
        *L: maximum sub-sequence length to be analyzed.
    Output: 
        *conditional_probabilities: a list of dictionaries. Each dictionary 
         contains keys that are of the form:
         symbol|subsequence
         meaning the probability of "symbol" occuring after that subsequence.
         There is one dictionary for each length of subsequence.
    Description:
        Calculates the probability of each symbol in alphabet occuring each
        subsequence in probabilities and create a similiar dictionary for those
        conditional probabilities.
    '''     
    def calc_cond_probs(self, L):
        self.conditional_probabilities = []
        print "Calculating subsequence conditional probabilities for sequence at: " + self.sequence_path
        print "L = " + str(L)
        if self.probabilities:
            self.conditional_probabilities = [self.probabilities[0]]
            for l in range(0, L):
                print "Sequence: " + self.sequence_path
                print "Calculating conditional probabilities of subsequences of length: " + str(l)
                d = {}
                l1 = self.probabilities[l]
                l2 = self.probabilities[l+1]
                for s in l1:
                    acc = 0
                    conds = []
                    for a in self.alphabet:
                        cond = a + "|" + s
                        t = s + a
                        if t in l2.keys():
                            d[cond] = l2[t]/l1[s]
                            acc += d[cond]
                        else:
                            d[cond] = 0.0
                    if not acc == 1.0:
                        dif = (1.0 - acc)/len(self.alphabet)
                        for c in conds:
                            d[c] += dif
                self.conditional_probabilities.append(d)
        else:
            print "Probabilities not computed."
            print "Run calc_probs function before this one."
        print "*****************"
        print "Sequence: " + self.sequence_path
        print "Conditional probabilities calculated!"
        print "*****************"
        return self.conditional_probabilities
        
    '''
    Name: calc_cond_entropy
    Input:
        *L: maximum sub-sequence length to be analyzed.
    Output: 
        *conditional_entropy: An array of conditional entropies for subsequences
         of length up to L.
    Description:
        Calculates the conditional entropy of each symbol occuring given sub-
        sequences of lengths up to L.
    '''    
    def calc_cond_entropy(self, L):
        self.cond_entropy = []
        print "Calculating conditional entropy for sequence at: " + self.sequence_path
        print "L = " + str(L)
        if self.probabilities: 
            if self.conditional_probabilities:
                for l in range(0, L +1):
                    print "Sequence: " + self.sequence_path
                    print "Calculating conditional entropy of length: " + str(l)
                    acc = 0
                    p = self.probabilities[l]
                    pcond = self.conditional_probabilities[l]
                    for x in p.keys():
                        if l == 0:
                            acc -= p[x]*num.log2(p[x])
                        else:
                            y_given_x = x[-1] + '|' + x[0:-1]
                            if not pcond[y_given_x] == 0:
                                acc -= p[x]*num.log2(pcond[y_given_x])
                    self.cond_entropy.append(acc)
            else:
                print "Conditional probabilities not computed."
                print "Run calc_cond_probs function before this one." 
        else:
            print "Probabilities not computed."
            print "Run calc_probs function before this one."
        print "*****************"
        print "Sequence: " + self.sequence_path
        print "Conditional entropy calculated!"
        print "*****************"
        return self.cond_entropy        
    
    '''
    Name: calc_kldivergence
    Input:
        *base_probs: A list of probability dictionaries to which the 
         probabilities contained in this class will be compared.
        *K: The length/level of probabilities from each that will be compared.
    Output: 
        *kldivergence: The Kullback-Leibler Divergence between the probability
         distributions of sequences of length K from base_probs and
         self.probabilities.
    Description:
        Calculates the KL Divergence of prob distributions of K-length seqs.
    '''       
    def calc_kldivergence(self, base_probs, K):
        self.kldivergence = 0
        print "Calculating Kullback-Leibler divergence for sequence at: " + self.sequence_path
        print "K = " + str(K)
        if self.probabilities:
            for key in base_probs[K].keys():
                p = base_probs[K][key]
                if key in self.probabilities[K].keys():
                    q = self.probabilities[K][key]
                else:
                    q = 1e-15 #Default non-zero really small value
                self.kldivergence += p*num.log2(p/q)
        else:
            print "[error] Probabilities not computed."
            print "Run calc_probs function before this one."
        print "*****************"
        print "Sequence: " + self.sequence_path
        print "Kullback-Leibler divergence calculated!"
        print "*****************"
        return self.kldivergence              
     
    '''
    Name: calc_autocorrelation
    Input:
        *up_to: The maximum window to be used in the autocorrelation.
    Output: 
        *autocorrelation: the sequence's autocorrelation.
    Description:
        Calculates the correlation of the sequence with itself with a window 
        size up to the input variable.
    '''   
    def calc_autocorrelation(self, up_to):
        self.autocorrelation = []
        temp = []
        print "Calculating autocorrelation probabilities for sequence at: " + self.sequence_path
        print "Up to " + str(up_to)
        for i in range(0, up_to):
            acc = 0
            for j in range(0, len(self.seq) - i):
                acc += float(self.seq[i + j])*float(self.seq[j])
            temp.append(acc)
        m = max(temp)
        self.autocorrelation = [float(x)/(2*float(m)) for x in temp]
        print "*****************"
        print "Sequence: " + self.sequence_path
        print "Autocorrelation calculated!"
        print "*****************"
        return self.autocorrelation
    
    '''
    Name: calc_l1metric
    Input:
    	*base_probs: A list of probability dictionaries to which the 
         probabilities contained in this class will be compared.
        *up_to: The maximum window to be used in the metric calculation.
    Output: 
        *l1metric: the L1 metric to estimate divergence between sequences.
    Description:
        Calculates the L1 metric of the probability distributions of 
        subsequences of lengths up to the input parameter between the base
        sequence and the object's sequence.
    '''      
    def calc_l1metric(self, base_probs, up_to):
        self.l1metric = 0
        acc = 0
        print "Calculating l1-metric for sequence at: " + self.sequence_path
        print "Up to " + str(up_to)
        if self.probabilities:
            for i in range(0, up_to):
                for key in base_probs[i].keys():
                    if key in self.probabilities[i].keys():
                        v2 = self.probabilities[i][key]
                    else:
                        v2 = 0
                    acc += abs(base_probs[i][key] - v2)/(2**i)
            self.l1metric = acc
        else:
            print "Probabilities not computed."
            print "Run calc_probs function before this one."
        print "*****************"
        print "Sequence: " + self.sequence_path
        print "L1-Metric calculated!"
        print "*****************"
        return self.l1metric
    
    '''
    Name: save_cond_probs_as_graph
    Input:
    	*path: path where the output graph will be saved.
    Output:
    Description:
        Saves the conditional probabilities in the format described in graph.py
        The conditional probabilities are used to form a Rooted Tree With Proba-
        bilities where a given subsequence will be a state and its transitions
        will be the conditional probabilities of each symbol of the alphabet
        given the current state.
    '''    
    def save_cond_probs_as_graph(self, path):
        print "Generating rooted tree with probabilities for sequence at: " + self.sequence_path
        if self.conditional_probabilities:
            i = 0
            states = []
            for level in self.conditional_probabilities:
                if i > 0:
                    names = [x[(x.index('|') + 1):] for x in level.keys()]
                    state_names = list(set(names))
                    for name in state_names:
                        outedges = []
                        for a in self.alphabet:
                            prob = str(level[a + '|' + name])
                            dest = name + a
                            oedge = (a, dest, prob)
                            outedges.append(oedge)
                        state = (name, outedges)
                        states.append(state)
                else:
                    name = 'e'
                    outedges = []
                    for k in level.keys():
                        prob = str(level[k])
                        oedge = (k, k, prob)
                        outedges.append(oedge)
                    state = (name, outedges)
                    states.append(state)
                i += 1
            with open(path, 'w') as file_:
                yaml.dump([states, self.alphabet], file_)
            print "*****************"
            print "Sequence: " + self.sequence_path
            print "Rooted Tree With Probabilities Generated!"
            print "*****************"
        else:
            print "Conditional probabilities not computed."
            print "Run calc_cond_probs function before this one."            
