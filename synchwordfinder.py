import probabilisticgraph as pg
import copy


class SynchWordFinder:
    def __init__(self, graph_path, w, l, alpha, test):
        self.w = w
        self.alpha = alpha
        self.test = test
        self.path = graph_path
        self.s = pg.ProbabilisticGraph(path='graphs/' + graph_path + '/rtp_L' + str(l) + '.yaml')
        self.candidacy_flags = {}
        self.tested_flags = {}
        self.t = copy.deepcopy(self.s)
        self.psi = {}
        for state in self.t.states:
            state.name = state.name[::-1]
            self.psi[state.name] = []
            self.candidacy_flags[state.name] = True
            self.tested_flags[state.name] = False
        e = self.s.root()
        self.gamma = [[e, self.candidacy_flags[e.name], self.tested_flags[e.name]]]
        self.delta = [self.t.root()] + self.t.root().obtain_children()
        self.omega_syn = []
        self.theta = []

    def next_valid_state(self):
        candidates = [x for x in self.gamma
                      if self.candidacy_flags[x[0].name] is True and self.tested_flags[x[0].name] is False]
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
                                                      c[0].name not in self.psi[s.name] and
                                                      (c[0].name, s.name) not in self.theta]
                    count = 0
                    for el in lamda:
                        candidate = c[0].name
                        suf = self.is_suffix(candidate, el.name[0:len(candidate)], self.t.root())
                        if suf:
                            self.psi[el.name].append(c[0].name)
                            p = self.s.compare_morphs(c[0].morph(), el.morph(), self.alpha, self.test)
                            self.theta.append((c[0].name, el.name))
                            if p[0]:
                                if count >= len(lamda):
                                    self.tested_flags[c[0].name] = True
                            else:
                                self.candidacy_flags[c[0].name] = False
                                self.expand_trees(c[0])
                                self.gamma.remove(c)
                                for g in self.gamma:
                                    self.tested_flags[g[0].name] = False
                                break
                        count += 1
                        if count == len(lamda):
                            self.tested_flags[c[0].name] = True
                else:
                    self.gamma.remove(c)
            else:
                self.omega_syn = [x[0].name for x in self.gamma if self.candidacy_flags[x[0].name] is True and
                                  self.tested_flags[x[0].name] is True]
                return self.omega_syn

    def is_suffix(self, candidate, string, state):
        if candidate == 'e':
            return True
        elif not state:
            return False
        elif string == '':
            return state.name == candidate
        else:
            nxt = state.next_state_from_edge(string[0])
            rest = string[1:]
            return self.is_suffix(candidate, rest, nxt)

    def shortest_valid_suffix(self, name):
        n = name[::-1]
        aux = self.t.root()
        i = 0
        while (self.candidacy_flags[aux.name] == False) and (i < len(n)):
            aux = aux.next_state_from_edge(n[i])
            i += 1
            if not aux:
                return aux
        return aux.name

    def expand_trees(self, c):
        children = c.obtain_children()
        new_nodes = [[x, self.candidacy_flags[x.name], self.tested_flags[x.name]] for x in children
                     if x.name == self.shortest_valid_suffix(x.name)]
        self.gamma.extend(new_nodes)
        rev = [x[0].name[::-1] for x in new_nodes]
        for r in rev:
            d = self.t.state_named(r)
            self.delta.extend(d.obtain_children())

