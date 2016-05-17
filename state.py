class State:
    '''
    This class represents a state or node in a graph. This representation is
    done using a name (or a label) to refer to the state and a list of its 
    outgoing edges. The outgoing edges list is a list of tuples where the first
    element is the letter from the automata/shift space alphabet labeling the
    edge and the second element is the name of the state to where that edge
    leads.
    There are methods to retrieve the destination of an outgoing edge based on
    the letter of the alphabet labeling that edge and, respectively, to retrieve
    a state name from the letter labeling an edge.
    '''
    
    def __init__(self, name, outedges):
        self.name = name            #The state's name/label
        self.outedges = outedges    #List of outgoing edges
        #An outgoing edge is a 2-tuple composed of:
        #(label, destination state)    
    
    '''
    Input: letter from the graph's alphabet
    Output: Destination from the edge containing the input letter as label.
    '''        
    def next_state_from_edge(self, label):    
        #Finds and returns a state name from outedges based on a letter.
        match = [x[1] for x in self.outedges if x[0] == label][0]
        return match
    
    '''
    Input: Destination state's name.
    Output: Label from the edge that goes to the desired state.
    '''    
    def edge_leads_to_state(self, state_name):    
        #Finds and returns an edge label from outedges based on a state name.
        match = [x[0] for x in self.outedges if x[1].name == state_name][0]
        return match
    '''
    Input: 
    Output: Returns all the destinations from all outgoing edges.
    '''      
    def obtain_children(self):
        children = [x[1] for x in self.outedges]
        return children
   
    '''
    Input: 
    Output: Returns the length of the state's name. It is redefined in order to
    return length 0 for the empty string.
    '''     
    def name_length(self):
        if self.name == 'e':
            return 0
        else:
            return len(self.name)
