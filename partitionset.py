import partition
import graph as gr
import state as st

'''
A Partition Set is a collection of partition, which, in turn are collections
of states which share some common feature. From a deterministic partition set
it is possible to recover a graph.
'''


class PartitionSet:
    def __init__(self, partitions):
        self.partitions = partitions

    # Not in current use:
    # def update_partitions_edges(self):
    #    for part in self.partitions:
    #        for other in self.partitions:
    #            part.update_edges(other)

    '''
    Name: recover_graph
    Input:
        *g: original graph that generated this partition set
    Output:
        *h: graph created by making the partition set into a graph
    Description:
    '''

    def recover_graph(self, g):
        states = [g.state_named(p.name[0]) for p in self.partitions]
        states = [x for x in states if x != None]  ##Just making sure no invalid states
        new_states = []
        for s in states:
            oedge = []
            for a in g.alphabet:
                t = s.next_state_from_edge(a)
                if t != None:
                    for p in self.partitions:
                        if t.name in p.name:
                            for e in s.outedges:
                                if e[0] == a:
                                    newedge = []
                                    i = 0
                                    for element in e:
                                        if i == 1:
                                            # The outedges are created pointing just
                                            # to a name, not to a state:
                                            newedge.append(p.name[0])
                                        else:
                                            newedge.append(element)
                                        i += 1
                                    newedge = tuple(newedge)
                                    oedge.append(newedge)
                                    break
                            break
            u = st.State(s.name, oedge)
            new_states.append(u)
        h = gr.Graph(new_states, g.alphabet)
        # As commented before, the edges were created pointing to names. They are
        # corrected here:
        corrected_states = h.reassign_dest_edges(h.states)
        h.states = corrected_states
        return h
