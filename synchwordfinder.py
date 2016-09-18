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
        self.tested = []
        self.t = copy.deepcopy(self.s)
        for state in self.t.states:
            state.name = state.name[::-1]
            self.candidacy_flags[state.name] = True
        self.gamma = [self.s.root()]
        self.delta = [self.t.root()] + self.t.root().obtain_children()
        self.valid_suffixes = {}
        self.suffixes = {}
        self.suffixes['e'] = self.t.root().obtain_children()

    def find_synch_words(self):
        while self.gamma:
            c = self.gamma.pop(0)
            l = c.name_length()
            if l < self.w:
                candidate = c.name
                lamda = self.suffixes[candidate] if candidate in self.suffixes.keys() else []
                if not lamda:
                        self.tested.append(c)
                else:
                    count = 0
                    for el in lamda:
                        for l2 in self.l2range:
                            p = self.s.compare_morphs(c.extended_morph(l2), el.extended_morph(l2),
                                                      self.alpha, self.test)
                            if p[0] == False:
                                break
                        if not p[0]:
                            self.candidacy_flags[c.name] = False
                            self.expand_trees(c)
                            self.gamma += self.tested
                            self.tested = []
                            break
                        count += 1
                    if count == len(lamda):
                        self.suffixes[candidate] = []
                        self.tested.append(c)
        return [x.name for x in self.tested]

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
        rev = self.expand_gamma(c)
        self.expand_delta(rev)

    def expand_gamma(self, c):
        rev = []
        gamma_children_states = c.obtain_children()
        if c.name in self.valid_suffixes.keys():
            gamma_children_states.extend(self.valid_suffixes[c.name])
            del self.valid_suffixes[c.name]
        for el in self.suffixes[c.name]:
            new_suf = self.shortest_valid_suffix(self.t.root(), el.name)
            if new_suf and new_suf.name != el.name:
                if new_suf.name in self.suffixes.keys():
                    self.suffixes[new_suf.name].append(el)
                else:
                    self.suffixes[new_suf.name] = [el]
        for gcs in gamma_children_states:
            short = self.shortest_valid_suffix(self.t.root(), gcs.name[::-1])
            if short:
                if short.name == gcs.name:
                    self.gamma.append(gcs)
                    rev.append(gcs.name[::-1])
                else:
                    if short.name in self.valid_suffixes.keys():
                        self.valid_suffixes[short.name].append(gcs)
                    else:
                        self.valid_suffixes[short.name] = [gcs]
        return rev

    def expand_delta(self, rev):
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