import probabilisticgraph as pg
import copy


class SynchWordFinder:
    def __init__(self, graph_path, w, l, alpha, test, l2range=[1]):
        self.w = w
        self.alpha = alpha
        self.test = test
        self.l2range = l2range
        self.path = graph_path
        self.s = pg.ProbabilisticGraph(path='graphs/' + graph_path + '/rtp_L' + str(l) + '.yaml')
        self.candidacy_flags = {}
        self.tested_flags = {}
        self.t = copy.deepcopy(self.s)
        for state in self.t.states:
            state.name = state.name[::-1]
            self.psi[state.name] = []
            self.candidacy_flags[state.name] = True
            self.tested_flags[state.name] = False
        e = self.s.root()
        self.gamma = [[e, self.candidacy_flags[e.name], self.tested_flags[e.name]]]
        self.delta = [self.t.root()] + self.t.root().obtain_children()
        self.omega_syn = []
        self.valid_suffixes = {}
        self.suffixes = {}
        self.suffixes['e'] = self.t.root().obtain_children()

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
                    candidate = c[0].name
                    lamda = self.suffixes[candidate] if candidate in self.suffixes.keys() else []
                    if not lamda:
                            self.tested_flags[candidate] = True
                    else:
                        count = 0
                        final = len(lamda)
                        for el in lamda:
                            for l2 in self.l2range:
                                p = self.s.compare_morphs(c[0].extended_morph(l2), el.extended_morph(l2),
                                                          self.alpha, self.test)
                                if p[0] == False:
                                    break
                            self.suffixes[candidate].remove(el)
                            if p[0]:
                                if count >= final:
                                    self.tested_flags[c[0].name] = True
                            else:
                                self.candidacy_flags[c[0].name] = False
                                self.expand_trees(c[0])
                                self.gamma.remove(c)
                                for g in self.gamma:
                                    self.tested_flags[g[0].name] = False
                                break
                            count += 1
                            if count == final:
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

    def shortest_valid_suffix(self, state, name):
        if not state:
            return None
        elif self.candidacy_flags[state.name]:
            return state
        else:
            nxt = state.next_state_from_edge(name[0])
            rst = name[1:]
            return self.shortest_valid_suffix(nxt, rst)

    def expand_trees(self, c):
        rev = []
        gamma_children_states = c.obtain_children()
        if c.name in self.valid_suffixes.keys():
            gamma_children_states.extend(self.valid_suffixes[c.name])
            del self.valid_suffixes[c.name]
        for gcs in gamma_children_states:
            short = self.shortest_valid_suffix(self.t.root(), gcs.name[::-1])
            if short:
                if short.name == gcs.name:
                    self.gamma.append([gcs, self.candidacy_flags[gcs.name], self.tested_flags[gcs.name]])
                    rev.append(gcs.name[::-1])
                else:
                    if short.name in self.valid_suffixes.keys():
                        self.valid_suffixes[short.name].append(gcs)
                    else:
                        self.valid_suffixes[short.name] = [gcs]
        for r in rev:
            d = self.t.state_named(r)
            d_children = d.obtain_children()
            self.delta.extend(d_children)
            for d_c in d_children:
                suf = self.shortest_valid_suffix(self.t.root(), d_c.name)
                if suf:
                    if suf.name in self.suffixes.keys():
                        self.suffixes[suf.name].append(d_c)
                    else:
                        self.suffixes[suf.name] = [d_c]