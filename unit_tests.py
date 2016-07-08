import probabilisticgraph as pg
from os import listdir
import sys
import getopt

class TestPFSACompleteness():

    def __init__(self, graphs_path):
        self.graph_list = listdir(graphs_path)
        self.graph_list = [graphs_path + '/' + p for p in self.graph_list if 'rtp' not in p]
        self.graphs = []
        for p in self.graph_list:
            g = pg.ProbabilisticGraph(path=p)
            self.graphs.append(g)

    def test_morph_completeness(self):
        succ = True
        for g in self.graphs:
            fail = False
            fail_states = []
            for s in g.states:
                dist = sum([float(e[2]) for e in s.outedges])
                dif = abs(1.0 - dist)
                if dif > 1e-9:
                    fail = True
                    succ = False
                    fail_states.append(s.name)
            if fail:
                print "**********"
                print "test_morph_completeness fail: Morphs do not add up to 1"
                print "At PFSA: " + self.graph_list[self.graphs.index(g)]
                print "At states: "
                print fail_states
                print "**********"
        if succ:
            print "**********"
            print "test_morph_completeness success!"
            print "**********"
        else:
            print "**********"
            print "test_morph_completeness failure!"
            print "**********"

    def test_reachable_states(self):
        succ = True
        for g in self.graphs:
            names = [s.name for s in g.states]
            fail = False
            fail_states = []
            for s in g.states:
                dests = [e[1] for e in s.outedges if e[1]]
                dests_names = [d.name for d in dests]
                if not set(names).issuperset(set(dests_names)):
                    fail = True
                    succ = False
                    fail_states.append(s.name)
            if fail:
                print "**********"
                print "test_reachable_states fail: States have edges to inexistent states"
                print "At PFSA: " + self.graph_list[self.graphs.index(g)]
                print "At states: "
                print fail_states
                print "**********"
        if succ:
            print "**********"
            print "test_reachable_states success!"
            print "**********"
        else:
            print "**********"
            print "test_reachable_states failure!"
            print "**********"

    def test_valid_states(self):
        succ = True
        for g in self.graphs:
            fail = False
            fail_states = []
            for s in g.states:
                dests = [e[1] for e in s.outedges if e[1]]
                if not dests:
                    fail = True
                    succ = False
                    fail_states.append(s.name)
            if fail:
                print "**********"
                print "test_valid_states fail: States with no outgoing edges."
                print "At PFSA: " + self.graph_list[self.graphs.index(g)]
                print "At states: "
                print fail_states
                print "**********"
        if succ:
            print "**********"
            print "test_valid_states success!"
            print "**********"
        else:
            print "**********"
            print "test_valid_states failure!"
            print "**********"


    def test_valid_morphs(self):
        succ = True
        for g in self.graphs:
            fail = False
            fail_states = []
            for s in g.states:
                dests = [e[1] for e in s.outedges if float(e[2]) > 0.0]
                for d in dests:
                    if not d:
                        succ = False
                        fail = True
                        fail_states.append(s)
            if fail:
                print "**********"
                print "test_valid_morphs fail: States with non zero probability to Nowhere."
                print "At PFSA: " + self.graph_list[self.graphs.index(g)]
                print "At states: "
                print fail_states
                print "**********"
        if succ:
            print "**********"
            print "test_valid_morphs success!"
            print "**********"
        else:
            print "**********"
            print "test_valid_morphs failure!"
            print "**********"

    def test_suite(self):
        self.test_morph_completeness()
        self.test_reachable_states()
        self.test_valid_states()
        self.test_valid_morphs()


def read_input(argv):
    path = ''
    try:
        opts, args = getopt.getopt(argv, "hp:", ["path="])
    except getopt.GetoptError:
        print "unit_tests.py -p <path>"
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print "unit_tests.py -p <path>"
            sys.exit()
        if opt in ("-p", "--path"):
            path = arg
    return path


def main(path):
    tester = TestPFSACompleteness(path)
    tester.test_suite()


if __name__ == "__main__":
    path = read_input(sys.argv[1:])
    main(path)
