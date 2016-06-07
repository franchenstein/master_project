import state as st
import json

class Graph:
    '''
    This class represents a graph by a list of states (its nodes) and a list of 
    the letters (the alphabet) from which it is constructed. As the states
    contain the information about the edges, a list of them fully describe a
    graph.
    There are methods to read a list of states from a file and to save it in the
    same format.
    '''
    
    def __init__(self, states, alphabet):
        self.states = states        #Graph's list of states
        self.alphabet = alphabet    #List of letters representing the alphabet
    
    '''
    Name: save_graph_file
    Input:
        *path: file path where the graph file will be saved.
    Description:
        Saves the graph as a json file. The file will contain a list of saved
        states and the graph alphabet. The saved state list will contain tuples
        of (state name, state outedges). The state outedges have to be modified
        as json cannot save the whole state. The destination state in the 
        outedges tuple is replaced just by the name of the state. The dest state
        will be recovered in the open graph file function.
    '''    
    def save_graph_file(self, path):
        savestates = []
        for s in self.states:
            savestates.append((s.name, s.outedges))
        with open(path, 'w') as file_:
            json.dump([savestates, self.alphabet], file_)
        return
    
    '''
    Name: open_graph_file
    Input:
        *path: file path where the graph file is saved.
    Description:
        Opens a graph saved in json format described above. To recover the dest
        state for each outedge, the state whose name is in the saved outedge
        will be searched for in the states list and it will substitute its name
        in the outedge.
    '''            
    def open_graph_file(self, path):
        with open(path, 'r') as file_:
            savedstates, alph = json.load(file_)
        self.alphabet = alph
        states = []
        for x in savedstates:
            s = st.State(x[0], tuple(x[1]))
            states.append(s)
        self.states = states
        return
        
    '''
    Name: reassign_dest_edges
    Input:
        *states: list of states with their outedges pointing only to a 
        destination state name, but not the whole state
    Output:
        *states: the corrected list, with the outedges correctly pointing to 
        the destination state.
    @staticmethod
    def reassign_dest_edges(states):
        for s in states:
            new_outedges = []
            for e in s.outedges:
                e_dest = [x for x in states if x.name == e[1]]
                if e_dest:
                    e_dest = e_dest[0]
                else:
                    e_dest = None
                new_e = []
                i = 0
                for element in e:
                    if i == 1:
                        new_e.append(e_dest)
                    else:
                        new_e.append(element)
                    i += 1
                new_e = tuple(new_e)
                new_outedges.append(new_e)
            s.outedges = new_outedges
        return states
    '''

    '''
    Name: root
    output:
        *The root of a rooted tree with probabilitites.
    Description:
        This function considers the element with label of length 0 as the root
        of a rooted tree with probabilities. The function will search for and
        return it.
    '''        
    def root(self):
        for s in self.states:
            if s.name == 'e':
                return s
            else:
                return st.State("", [])
                
    '''
    Name: state_named
    Input:
        *state_name: the name of the state which should be retrieved
    Output:
        *state: the state with the name given as input
    Description:
        Receives a state name as input, searches for a state with this name if
        it is found, returns None if it is not.
    '''    
    def state_named(self, state_name):
        state_names = [x.name for x in self.states]
        if state_name in state_names:
            i = state_names.index(state_name)
            return self.states[i]
        else:
            return None
            
    '''
    Name: remove_unreachable_states
    Input: A graph described by the class Graph
    Output: A graph where all unreachable states of the input graph are removed.
    Description: The algorithm goes through all the outedges from each state of 
    the graph. Based on the outedges, it creates a list with all destination 
    states. It then proceeds to run through the graph's states and create a list
    of them which only includes states whose names are on the list of reachable 
    states. The alphabet is updated accordingly. With those two elements, a new 
    reduced graph is created and returned. 
    '''        
    def remove_unreachable_states(self):
        
        old_size = len(self.states)
        reachable_states = [] #This will receive the reachable states' names
        #Creates a list of all states' outedges:
        aux = [x.outedges for x in self.states] 
        
        for outedges in aux: #Goes through each state's outedge list
            for outedge in outedges: #Goes through each outedge in the outedge list
                #Checks if the destination state of the current outedge is already
                #in the list:
                if outedge[1]:
                    s = self.state_named(outedge[1])
                    if s.name not in reachable_states:
                        #If it is not, it is considered as a new reachable state.
                        reachable_states.append(s.name)
                    
        #A new list of states is created only with states whose names are in 
        #reachableStates            
        new_states = [x for x in self.states if x.name in reachable_states]
        
        #List of outedges lists of the new states:
        aux = [x.outedges for x in new_states]
        new_alphabet = []  #Receives the new alphabet
        for outedges in aux:  #Goes through each state's outedge list
            for outedge in outedges:  #Goes through each outedge in the outedge list
                #Checks if the outedge label is already in the alphabet:
                if outedge[0] not in new_alphabet:
                    #If it's not, it is included to the new alphabet.
                    new_alphabet.append(outedge[0])
        
        #Creates a new graph, without previous unreachable states:
        reduced_graph = Graph(new_states, new_alphabet)
        newSize = len(reduced_graph.states)
        if (old_size != newSize):
            reduced_graph = reduced_graph.remove_unreachable_states()
        
        return reduced_graph

    def __str__(self):
        for s in self.states:
            print s
        r = '****************************************\n'
        r += 'Number of states: ' + str(len(self.states)) + '\n'
        return r

    def print_state_named(self, n):
        s = self.state_named(n)
        print s
