#!/usr/bin/env
import graph as gr
import partition as pt
import partitionset as ps
import state as st

'''
Name: default_initial_partition
Inputs: 
    *g: graph for which the initial partition has to be found.
Outputs: 
    *initial partition following the default rule. 
Description: The default rule for initial partitioning divides the original
graph in partitions containing states whose outgoing edges have the same
labels.
'''  
def default_initial_partition(g):
    partitions = []
    for s in g.states:
        s_edge_labels = set([oedge[0] for oedge in s.outedges])
        if not partitions:
            prt = pt.Partition(s)
            partitions.append(prt)
        else:
            fail_count = 0
            for p in partitions:
                p_edge_labels = set([oedge[0] for oedge in p.outedges])
                if p_edge_labels == s_edge_labels:
                    p.add_to_partition(s)
                else:
                    if fail_count < len(partitions):
                        fail_count += 1
                    else:
                        prt = pt.Partition(s)
                        partitions.append(prt)                
    return ps.PartitionSet(partitions)

'''
Name: moore
Inputs: 
    *graph: graph that has to be reduced;
    *inital_partition: the initial rough partition that moore will split.
Outputs: 
    *reduced_partitions: the minimal partitions that make a deterministic
     graph. 
Description: This function will call the moore_iterations while the current
partition is different from the partition obtained after a iteration. Once
they are equal, the reduced partitions are returned.
'''           
def moore(initial_partition, graph, simple=False):
    if simple:
        p = simplemoore(initial_partition, graph)
        return ps.PartitionSet(p)
    else:
        current_partition = initial_partition
        while True:
            new_partition = moore_iteration(graph, current_partition.partitions)
            if len(current_partition.partitions) == len(new_partition):
                return ps.PartitionSet(new_partition)
            else:
                current_partition = ps.PartitionSet(new_partition)

'''
Name: moore_by_parts
Inputs: 
    *graph: graph that has to be reduced;
    *inital_partition: the initial rough partition that moore will split;
    *n_iter: number of times the algorithm will iterate.
Outputs: 
    *reduced_partitions: the minimal partitions that make a deterministic
     graph. 
Description: Similar to the original moore, but it does not wait for the
new partition to be equal to the one before the iteration, but has a pre-
determined amount of iterations that will be run.
'''     
def moore_by_parts(graph, initial_partition, n_iter = -1):
    #By default, it applies the regular moore algorithm
    if niter == -1:
        current_partition = moore(graph, initial_partition)
    else:
        current_partition = initial_partition
        for i in range(0, n_iter):
            current_partition = moore_iteration(graph, current_partition)
    return current_partition                 

'''
Name: moore_iteration
Inputs: 
    *graph: graph that has to be reduced;
    *current_partition: the partition for the current iteration
Outputs: 
    *new_partition: the new partition produced after a iteration. 
Description: This method will apply only one iteration of the Moore
algorithm. An explanation for the algorithm can be found at: 
http://arxiv.org/abs/1010.5318
'''                 
def moore_iteration(graph, current_partition):
    partition_for_alphabet = []
    for a in graph.alphabet:
        splits = []
        for p in current_partition:
            new_splits = splitting(p, a, graph.states)
            valid_splits = [sp for sp in new_splits if sp.name]
            splits.append(valid_splits)      
        partition_for_letter = splits[0]
        for split in splits[1:]:
            partition_for_letter = coarsest_partition(partition_for_letter, 
                                                      split)
        partition_for_alphabet.append(partition_for_letter)        
    final_partition = partition_for_alphabet[0]
    for pb in partition_for_alphabet[1:]:
        final_partition = coarsest_partition(final_partition, pb)
    new_partition = coarsest_partition(current_partition, final_partition)
    return new_partition 
        
'''
Name: splitting
Inputs: 
    *partition: the partition that has to be split;
    *letter: letter that will be taken into consideration to split the
     partion;
    *states: the list of states from the original graph.
Outputs: 
    *[p1, p2]: the two partitions obtained. 
Description: The method iterated for each state in the states. First it 
checks if the current state has an outgoing edge labeled by the letter
parameter. If it does not, the current state is added to both splits. If
it does, it checks if that outgoing edge reaches a state in the current 
partition. If it does, the current state is added to the first split. If it
does not, it is added to the second split.
'''     
def splitting(partition, letter, states):
    p1 = pt.Partition(st.State("", []))
    p2 = pt.Partition(st.State("", []))
    for s in states:
        edge_labels = [edge[0] for edge in s.outedges]
        if letter in edge_labels:
            next_state = s.next_state_from_edge(letter)
            if next_state:
                if next_state.name in partition.name:
                    p1.add_to_partition(s)
                else:
                    p2.add_to_partition(s)
            else:
                p2.add_to_partition(s)
        else:
            p1.add_to_partition(s)
            p2.add_to_partition(s)
    return [p1, p2]

'''
Name: coarsest_partition
Inputs: 
    *partitions 1 & 2: partitions for which we wish to find the coarsest
     partition.
Outputs: 
    *coarse: the coarsest partition between partitions 1 & 2. 
Description: Applies the intersection function between each pair of sets in
partition 1 and partition 2. Checks for redundancy. The set of all the 
unique intersections is the coarsest partition.
'''
def coarsest_partition(partition1, partition2):
    coarse = []
    for p1 in partition1:
        for p2 in partition2:
            inter = intersection(p1, p2)
            if inter.name: #disregards empty intersections
                if coarse:
                    #to avoid redundancy:
                    names = [el.name for el in coarse]
                    if inter.name not in names:
                        coarse.append(inter)
                else:
                    coarse.append(inter)
    return coarse
    
'''
Name: intersection
Inputs: 
    *p1 and p2: sets of states for which we wish to find the intersection.
Outputs: 
    *p: intersection between p1 and p2
Description: Finds and returns the states in common between partitions p1
and p2.
'''      
def intersection(p1, p2):
    inter = []
    for a in p1.name:
        for b in p2.name:
            if a == b:
                inter.append(b)
                break
    if inter:
        p = pt.Partition(st.State(inter.pop(0), []))
        for i in inter:
            p.add_to_partition(st.State(i, []))
    else:
        p = pt.Partition(st.State('', []))
    return p

def simplemoore(pset, graph):
    partition = pset.partitions
    while True:
        old_partition = partition
        s0 = graph.state_named(partition[0].name[0])
        ed = []
        for e in s0.outedges:
            for p in partition:
                if e[1]:
                    if e[1].name in p.name:
                        ed.append([e[0], p])
                        break
        p0 = [[[s0], ed]]
        for p in partition[1:]:
            for nm in p.name:
                s = graph.state_named(nm)
                fail = True
                for q in p0:
                    for e in s.outedges:
                        for d in q[1]:
                            if e[0] == d[0]:
                                if e[1]:
                                    if e[1].name in d[1].name:
                                        q[0].append(s)
                                        fail = False
                                        break

                if fail:
                    ed = []
                    for e in s.outedges:
                        for o in partition:
                            if e[1]:
                                if e[1].name in o.name:
                                    ed.append([e[0], o])
                    p0.append([[s], ed])
        partition = []
        for p in p0:
            for i in range(len(p[0])):
                if i == 0:
                    eq_class = pt.Partition(p[0][i])
                else:
                    eq_class.add_to_partition(p[0][i])
            partition.append(eq_class)
        if len(partition) == len(old_partition):
            return partition
