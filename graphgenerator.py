#!/usr/bin/env
import probabilisticgraph as pg
import partition as pt
import partitionset as ps
import moore as mr
import probabilisticstate as pst

'''
This class contains the methods for creating the final model for the original
sequence from the rooted tree with probabilities. It will take one such tree
and the desired parameters and apply one of the graph generation methods to
reduce the tree to a minimal graph. Any graph that can be reduced might have its
end level terminated with three types of termination methods:

I - The state in the last level goes through statistical tests with the states
of each of its suffixes and is substituted by the one with the highest p result.

II - The state in the last level goes through statistical tests with the states
of each of its suffixes and is substituted by the longest state which passes the
test.

III - A D-Markov termination.
'''
class GraphGenerator():
    #Initialization:
    def __init__(self, original_path, synch_words, save_path):
        self.original_graph = pg.ProbabilisticGraph(path=original_path)
        self.synch_words = self.set_synch_words(synch_words)
        self.save_path = save_path

    def set_synch_words(self, s_words):
        synch_words = []
        for w in s_words:
            synch_words.append(self.original_graph.state_named(w))
        return synch_words

    '''
    Name: mk1
    Inputs:
        *test: Statistical test type (chi-squared or Kolmogorov-Smirnov);
        *alpha: Statistical test quality parameter.
    Outputs:
        *reduced_graph = Minimal graph for the original sequence.
    Description: It starts with the synchronization word state as root and
    creates a class for it. In each iteration, the states expanded in the last
    iteration are then expanded. Each of the new states goes through a
    statistical test with the present classes. If the test passes, the new state
    is added to this class. If it fails, a new class is created with the new
    state. After this is done for all states, a graph reduction technique is
    applied to the partition set and a minimal graph is returned.
    '''
    def mk1(self, test, alpha):
        reduced_graph = self.apply_moore(test, alpha)
        reduced_graph.save_graph_file(self.save_path + '_mk1.yaml')
        return reduced_graph

    '''
    Name: mk2
    Inputs:
    Outputs:
    Description: This algorithm starts with all synchronization word states. It
    then expands each state. If a new state ends with a synchronization word, it
    is eliminated and edges connecting to it are then connected to the
    respective synchronization word.
    '''
    def mk2(self):
        synchlist = [x.name for x in self.synch_words]
        s = list(self.synch_words)
        newstates = []
        while True:
            if not s:
                break
            else:
                aux = s.pop(0)
                nexts = [x[1] for x in aux.outedges]
                i = 0
                new_outedges = []
                for n in nexts:
                    newdest = n
                    if n:
                        for w in synchlist:
                            if self.is_suffix(w, n.name):
                                newdest = [x for x in self.synch_words if x.name == w][0]
                                break
                    if newdest:
                        oedge = (aux.outedges[i][0], newdest, aux.outedges[i][2])
                    else:
                        oedge = (aux.outedges[i][0], None, aux.outedges[i][2])
                    new_outedges.append(oedge)
                    i += 1
                new_state = pst.ProbabilisticState(aux.name, new_outedges)
                newstates.append(new_state)
                states_names = [x.name for x in newstates]
                new_children = new_state.obtain_children()
                new_children = [x for x in new_children if x.name not in states_names]
                s.extend(new_children)
        new_graph = pg.ProbabilisticGraph(newstates, self.original_graph.alphabet)
        new_graph.save_graph_file(self.save_path + '_mk2.yaml')
        return new_graph

    def mk2_moore(self, test, alpha):
        self.original_graph = self.mk2()
        reduced_graph = self.apply_moore(test, alpha)
        reduced_graph.save_graph_file(self.save_path + '_mk2_moore.yaml')
        return reduced_graph

    def equivalence_classes(self, test, alpha):
        states = [s for s in self.original_graph.states if s not in self.synch_words]
        new_states = []
        while True:
            if states:
                s = states.pop(0)
            else:
                break
            found = False
            for x in self.synch_words:
                if self.is_prefix(x.name, s.name):
                    r = self.original_graph.compare_morphs(x.morph(), s.morph(), alpha, test)
                    if r:
                        dests_s = [(e[0], e[1]) for e in s.outedges if e[1]]
                        dests_x = [(e[0], e[1]) for e in x.outedges if e[1]]
                        if not set(dests_x).isdisjoint(dests_s):
                            for y in self.original_graph.states:
                                i = 0
                                for edge in y.outedges:
                                    if edge[1] == s:
                                        new_edge = (edge[0], x, edge[2])
                                        y.outedges[i] = new_edge
                                    i += 1
                            found = True
                            break
            if found:
                if s in states:
                    states.remove(s)
                if s in self.original_graph.states:
                    self.original_graph.states.remove(s)
            new_states.append(s)
        new_graph = self.apply_moore(test, alpha)
        return new_graph

    def mk3(self, test, alpha):
        self.original_graph = self.mk2()
        synch_names = [w.name for w in self.synch_words]
        self.synch_words = self.set_synch_words(synch_names)
        new_graph = self.equivalence_classes(test, alpha)
        reduced_graph = self.apply_moore(test, alpha)
        reduced_graph.save_graph_file(self.save_path + '_mk3.yaml')
        return new_graph

    '''
    Name: create_initial_partition
    Inputs:
        *init_state: Initial state where to start expansion;
        *test: Statistical test type (chi-squared or Solmogorov-Smirnov);
        *alpha: Statistical test quality parameter.
    Outputs:
        *P: Partition set where all states in a give partition have statistically
         similar morphs
    Description: Creates the initial partition to be used in the mk1 algorithm.
    It will start at a given state, expand it and compare its children to the
    available partitions and add the child to it if the test is succesful, but
    create a new partition for it if is not succesful for any partition.
    '''
    def create_initial_partition(self, init_state, alpha, test):
        partition_0 = pt.Partition(init_state)
        children = init_state.obtain_children()
        partitions = [partition_0]
        while True:
            if children:
                c = children.pop(0)
                if c:
                    fail_count = 0
                    for p in partitions:
                        fail_count += 1
                        pmorph = self.partition_morph(p.outedges[0])
                        result = self.original_graph.compare_morphs(pmorph, c.morph(), alpha, test)
                        if result[0]:
                            p.add_to_partition(c)
                            break
                        else:
                            if fail_count == len(partitions):
                                new_partition = pt.Partition(c)
                                partitions.append(new_partition)
                    new_states = c.obtain_children()
                    for n in new_states:
                        not_in_partitions = True
                        for p in partitions:
                            if n.name in p.name:
                                not_in_partitions = False
                        if not_in_partitions:
                            children.append(n)
            else:
                return partitions

    def apply_moore(self, test, alpha):
        p = self.create_initial_partition(self.synch_words[0], alpha, test)
        partition_set = ps.PartitionSet(p)
        reduced_classes = mr.moore(partition_set, self.original_graph)
        reduced_graph = reduced_classes.recover_graph(self.original_graph)
        reduced_graph = pg.ProbabilisticGraph(reduced_graph.states, reduced_graph.alphabet)
        reduced_graph.reassign_dest_edges(reduced_graph.states)
        reduced_graph = reduced_graph.remove_unreachable_states()
        return reduced_graph

    '''
    Name: is_suffix
    Inputs:
        *w: string, possible suffix of n;
        *n: string for which w might be a suffix.
    Outputs:
        *P: true if w is suffix of n, false otherwise.
    '''
    @staticmethod
    def is_suffix(w,n):
        if len(w) > len(n):
            return False
        else:
            nSuffix = n[-len(w):]
            return w == nSuffix

    @staticmethod
    def is_prefix(w, n):
        if len(w) > len(n):
            return False
        else:
            n_prefix = n[:len(w)]
            return w == n_prefix

    '''
    Name: partition_morph
    Inputs:
        *partition_oedges: the first element of the
         partition outedges.
    Outputs:
        *morph: the probability distribution associated
        to the partition outedges.
    Description:
        I don't think this method should be in this class.
        I need the morph from a partition, but the partition
        class is not specific to probabilistic states.
        I will leave this here until a better solution shows up.
    '''
    @staticmethod
    def partition_morph(partition_oedges):
        morph = []
        for edge in partition_oedges:
            morph.append((edge[0], edge[2]))
        return morph
