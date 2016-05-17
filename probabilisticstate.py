import state
import random as rnd
'''
Probabilistic version of state class. The outedges include a third parameter,
which reflects the probability of that outgoing edge being taken. Methods 
regarding this new parameter are added.
'''
class ProbabilisticState(state.State):
    
    def __init__(self, name, outedges):
        state.State.__init__(self, name, outedges)
    
    '''
    Name: prob_to_next_state
    Input: 
        *state_name: The destination state's label.
    Output: 
        *match: probability of reaching this state from the current state.
    '''     
    def prob_to_next_state(self, state_name):
        match = [x[2] for x in self.outedges if x[1] == state_name][0]
        return float(match)
    
    '''
    Name: prob_to_next_letter
    Input:
        *letter: letter from the graph's alphabet.
    Output:
        *match: probability of outputing the input letter.
    '''     
    def prob_to_next_letter(self, letter):
        match = [x[2] for x in self.outedges if x[0] == letter][0]
        return float(match)
    
    '''
    Name: morph
    Input: 
    Output: 
        *m: The state's morph, i.e. the probability distribution of its outputs.
    Description: Statistically comparing state's morphs will be used many times.
    The morph is a 2-tuple with only the output letter and the probability that
    it occurs, disregarding the destination state.
    '''     
    def morph(self):
        m = [(x[0], x[2]) for x in self.outedges]
        return m
   
    '''
    Name: randomstep
    Input:
    Output: 
        *A 2-tuple of the randomly chosen output symbol and the destination 
         state.
    Description:
        Takes on step in a walk through the graph. It randomly generates a
        real number in [0,1] and compares it to the state's morph and chooses
        an output symbol and destination accordingly. If, by some error, the
        randomly generated number does not fall into the distribution, an error
        2-tuple is returned.
    '''     
    def randomstep(self):
        r = rnd.random()
        acc = 0
        for e in outedges:
            if acc <= r < (acc + float(e[2])):
                return e[:2]
            else:
                acc += float(e[2])
        return ('', None) #Error 2-tuple.
