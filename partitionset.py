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

    def recover_graph(self, g, base_probs=None):
        states = []
        for p in self.partitions:
            if base_probs:
                probs = {k:v for k,v in base_probs.iteritems() if k in p.names}
                total_probs = sum(probs.values())
            nm = p.names[0]
            edges = []
            for a in g.alphabet:
                if base_probs:
                    transition = 0
                    for i in range(len(p.names)):
                        for e in p.outedges[i]:
                            if e[0] == a:
                                tr = e[2]
                                if i == 0:
                                    dest = e[1].name
                                break
                        transition += probs[p.names[i]]*tr
                    transition /= total_probs
                    edges.append([a, dest, transition])
                else:
                    for e in p.outedges[0]:
                        if e[0] == a:
                            dest = e[1].name
                    edges.append([a, dest])
            states.append(st.State(nm, edges))
        h = gr.Graph(states, g.alphabet)
        for s in h.states:
            for e in s.outedges:
                e[1] = h.state_named(e[1])
        return h
        # states = [g.state_named(p.name[0]) for p in self.partitions if g.state_named(p.name[0])]
        # new_states = []
        # for s in states:
        #     oedge = []
        #     for a in g.alphabet:
        #         t = s.next_state_from_edge(a)
        #         if t:
        #             for p in self.partitions:
        #                 if t.name in p.name:
        #                     dest = p.name[0]
        #                     break
        #                 else:
        #                     dest = ''
        #         else:
        #             dest = ''
        #         for e in s.outedges:
        #             if e[0] == a:
        #                 newedge = []
        #                 i = 0
        #                 for element in e:
        #                     if i == 1:
        #                         # The outedges are created pointing just
        #                         # to a name, not to a state:
        #                         newedge.append(dest)
        #                     else:
        #                         newedge.append(element)
        #                     i += 1
        #                 newedge = tuple(newedge)
        #                 oedge.append(newedge)
        #     u = st.State(s.name, oedge)
        #     new_states.append(u)
        # h = gr.Graph(new_states, g.alphabet)
        # return h
