import json
import numpy as num

'''
This class opens a sequence saved in a json file and has methods to compute 
various parameters based on said sequence.
'''
class SequenceAnalyzer():
    
    def __init__(self, path, probabilities = [], alphabet = [], 
                 conditional_probabilities = [], conditional_entropy = [],
                 kldivergence = 0, l1metric = 0, autocorrelation = []):
        self.probabilities = probabilities
        self.alphabet = alphabet
        self.conditional_probabilities = conditional_probabilities
        self.conditional_entropy = conditional_entropy
        self.kldivergence = kldivergence
        self.l1metric = l1metric
        self.autocorrelation = autocorrelation
        with open(path, 'r') as file_:
            self.seq = json.load(file_)
    
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
        l = 1
        self.probabilities = []
        self.alphabet = []
        while l <= L:
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
            self.probabilities.append(current_probs)
            l += 1
        return [self.probabilities, self.alphabet] 
    
    '''
    Name: calc_cond_probs
    Input:
        *L: maximum sub-sequence length to be analyzed.
    Output: 
        *conditional)probabilities: a list of dictionaries. Each dictionary 
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
        if self.probabilities:
            self.conditional_probabilities = [self.probabilities[0]]
            l = 0
            while l < L:
                d = {}
                l1 = self.probabilities[l]
                l2 = self.probabilities[l+1]
                for s in self.alphabet:
                    for a in l1:
                        cond = s + "|" + a
                        t = a + s
                        if t in l2.keys():
                            d[cond] = l2[t]/l1[a]
                        else:
                            d[cond] = 0.0
                self.conditional_probabilities.append(d)
                l += 1            
        else:
            print "Probabilities not computed."
            print "Run calc_probs function before this one." 
        return self.conditional_probabilities
        
    def calc_cond_entropy(self, L):
        self.cond_entropy = [] 
        try: 
            try:   
                l = 0
                while l <= L:
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
                    l += 1
            except:
                print "Conditional probabilities not computed."
                print "Run calc_cond_probs function before this one." 
        except:
            print "Probabilities not computed."
            print "Run calc_probs function before this one."                 
        return self.cond_entropy        
        
    def calc_kldivergence(self, base_probs, K):
        self.kldivergence = 0
        try:
            for key in base_probs[K].keys():
                p = base_probs[K][key]
                if key in self.probabilities[K].keys():
                    q = self.probabilities[K][key]
                else:
                    q = 1e-15
                self.kldivergence += p*num.log2(p/q)
        except:
            print "[error] Probabilities not computed."
            print "Run calc_probs function before this one." 
        return self.kldivergence              
        
    def calc_autocorrelation(self, up_to):
        self.autocorrelation = []
        temp = []
        for i in range(0, up_to):
            acc = 0
            for j in range(0, len(self.seq) - i):
                acc += float(self.seq[i + j])*float(self.seq[j])
            temp.append(acc)
        m = max(temp)
        self.autocorrelation = [float(x)/(2*float(m)) for x in temp]
        return self.autocorrelation
        
    def calc_l1metric(self, base_probs, up_to):
        self.l1metric = 0
        acc = 0
        try:
            for i in range(0, up_to):
                for key in base_probs[i].keys():
                    if key in self.probabilities[i].keys():
                        v2 = self.probabilities[i][key]
                    else:
                        v2 = 0
                    acc += abs(base_probs[i][key] - v2)/(2**i)
            self.l1metric = acc
        except:
            print "Probabilities not computed."
            print "Run calc_probs function before this one." 
        return self.l1metric
        
    def save_cond_probs_as_graph(self, path):
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
                json.dump([states, self.alphabet], file_)
        else:
            print "Conditional probabilities not computed."
            print "Run calc_cond_probs function before this one."            
