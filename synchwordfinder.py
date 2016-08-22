import probabilisticgraph as pg
import copy


class SynchWordFinder:
    def __init__(self, graph_path, w, l, alpha, test):
        self.w = w
        self.alpha = alpha
        self.test = test
        self.path = graph_path
        self.s = pg.ProbabilisticGraph(path='graphs/' + graph_path + '/rtp_L' + str(l) + '.yaml')
        self.t = copy.deepcopy(self.s)
        for state in self.t.states:
            state.name = state.name[::-1]
        self.gamma = [[self.s.root(), True, False]]
        self.delta = [self.t.root()] + self.t.root().obtain_children()
        self.omega_syn = []
        self.theta = []
        self.psi = []

    def next_valid_state(self):
        candidates = [x for x in self.gamma if x[1] is True and x[2] is False]
        if candidates:
            return candidates[0]
        else:
            return []

    def find_synch_words(self):
        while True:
            c = self.next_valid_state()
            if c:
                l = c[0].name_length()
                if l < self.w:
                    lamda = [s for s in self.delta if s.name_length() > l and
                                                      s.name not in self.psi and
                                                      (c[0].name, s.name) not in self.theta]
                    count = 0
                    for el in lamda:
                        suf = self.is_suffix(c[0].name, el.name)
                        if suf:
                            self.psi.append(el.name)
                            p = self.s.compare_morphs(c[0].morph(), el.morph(), self.alpha, self.test)
                            self.theta.append((c[0].name, el.name))
                            if p[0]:
                                if count >= len(lamda):
                                    c[2] = True
                            else:
                                c[1] = False
                                self.expand_trees(c[0])
                                self.gamma.remove(c)
                                for g in self.gamma:
                                    g[2] = False
                                self.psi = []
                                break
                        count += 1
                        if count == len(lamda):
                            c[2] = True
                else:
                    self.gamma.remove(c)
            else:
                self.omega_syn = [x[0].name for x in self.gamma if x[1] is True]
                return self.omega_syn

    def is_suffix(self, candidate, full):
        if candidate == 'e':
            return True
        else:
            partial = full
            aux = self.s.root()
            for i in range(len(candidate)):
                aux = aux.next_state_from_edge(partial[i])
                if not aux:
                    return False
            n = [x for x in self.gamma if x[0].name == aux.name]
            if n:
                return (aux.name is candidate) and (n[0][1])
            else:
                return False

    def shortest_valid_suffix(self, name):
        n = name
        aux = self.gamma[0]
        i = 0
        while (aux[1] == False) and (i < len(n)):
            s = aux[0].next_state_from_edge(n[i])
            aux = [x for x in self.gamma if x[0].name == s.name][0]
            i += 1
            if not aux:
                return aux
        return aux[0].name

    def expand_trees(self, c):
        children = c.obtain_children()
        self.gamma.extend([[x, True, False] for x in children])
        rev = [x.name[::-1] for x in children]
        for r in rev:
            n = self.shortest_valid_suffix(r)
            if n == r:
                d = self.t.state_named(n)
                d_children = d.obtain_children()
                for dc in d_children:
                    if not dc in self.delta:
                        self.delta.append(dc)

