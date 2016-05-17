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
    #Class attributes:
    self.original_graph = pg.ProbabilisticGraph([],[])
    self.synch_words = []
    self.save_path = ''
    
    #Class methods:
    
    #Initialization:
    def __init__(self, original_path, synch_words, save_path):
        self.original_graph.open_graph_file(original_path)
        self.synch_words = synch_words
        self.save_path = save_path
        
    
    '''    
    Name: mk1
    Inputs:
        *test: Statistical test type (chi-squared or Solmogorov-Smirnov);
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
        p = self.create_initial_partition(self.synch_words[0], alpha, test)
        partition_set = ps.PartitionSet(p)
        reduced_classes = mr.moore(partition_set, self.original_graph)
        reduced_graph = reduced_classes.recover_graph(self.original_graph)
        reduced_graph.save_graph_file(self.save_path)
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
        s = self.synch_words   
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
                    for w in synchlist:
                        newdest = n
                        if self.is_suffix(w, n.name):
                            newdest = [x for x in synch_words 
                                       if x.name == w][0]                                       
                            break
                    oedge = (aux.outedges[i][0], newdest, aux.outedges[i][2])
                    new_outedges.append(oedge)
                    i += 1    
                new_state = pst.ProbabilisticState(aux.name, new_outedges)
                newstates.append(new_state)
                states_names = [x.name for x in newstates]
                new_children = new_state.obtain_children()
                s.extend(new_children)
        new_graph = pg.ProbabilisticGraph(newstates, 
                                          self.original_graph.alphabet)
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
         children = self.original_graph.extend(init_state)
         partitions = [Q]
         while True:
            if children:
                c = children.pop(0)
                fail_count = 0
                for p in partitions:
                    result = self.original_graph.compare_morphs(p.outedges[0], 
                                                                c.outedges,
                                                                alpha, test)
                    if result:
                        p.add_to_partition(c)
                        break
                    else:
                        if fail_count < len(partitions):
                            fail_count += 1
                        else:
                            new_partition = pt.Partition(c)
                            partitions.append(new_partition)
                new_states = self.original_graph_extend(c)
                for n in new_states:
                    not_in_partitions = True
                    for p in partitions:
                        if n.name in p.name:
                            not_in_partitions = False
                    if not_in_partitions:
                        children.append(n)                   
            else:
                return P
    
    '''
    Name: is_suffix
    Inputs: 
        *w: string, possible suffix of n;
        *n: string for which w might be a suffix.
    Outputs: 
        *P: true if w is suffix of n, false otherwise.
    '''             
    def __is_suffix(self,w,n):
        if len(w) > len(n):
            return False
        else:
            nSuffix = n[-len(w):]
            return w == nSuffix
