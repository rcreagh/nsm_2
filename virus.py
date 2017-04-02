import random

import networkx as nx

def advance():
    """Advance the system by one time period.
    """
    
    #first loop captures changes made in a time period but does not apply them
    for V in G.nodes():
        # only living, infected nodes can affect population
        if infected[V] and not dead[V]:
            for N in G.neighbors(V):
                if not infected[N] and not immune[N]:
                    spam = random.random()
                    if spam < transmit * G.get_edge_data(V, N)['weight']:
                        newly_infected[N] = True
            eggs = random.random()
            # infected nodes can recover or die or stay infected
            if eggs < recover:
                newly_recovered[V] = True
            elif eggs < recover + die:
                newly_dead[V] = True
                
    # the second loop applies changes at the end of the time period            
    for V in G.nodes():
        if newly_infected[V]:
            infected[V] = True
            newly_infected[V] = False
        if newly_recovered[V]:
            infected[V] = False
            newly_recovered[V] = False
        if newly_dead[V]:
            dead[V] = True
            newly_dead[V] = False
            for E in G.edges(V):
                G.remove_edge(*E)

# define virus
transmit = 0.8
recover = 0.2
die = 0.1

#define population
G = nx.Graph()
E = (('A', 'B', 0.7),
     ('A', 'C', 0.8),
     ('A', 'D', 0.9),
     ('B', 'G', 0.7),
     ('B', 'C', 0.3),
     ('C', 'D', 0.6),
     ('C', 'E', 0.8),
     ('C', 'G', 0.4),
     ('C', 'F', 0.5),
     ('D', 'F', 0.9),
     ('G', 'E', 0.5),
     ('E', 'F', 0.3))
G.add_weighted_edges_from(E)
dead = {}
immune = {}
infected = {}
newly_infected = {}
newly_recovered = {}
newly_dead = {}
for V in G.nodes():
    dead[V] = False
    immune[V] = False
    infected[V] = False
    newly_infected[V] = False
    newly_recovered[V] = False
    newly_dead[V] = False    
infected['A'] = True
immune['C'] = True
immune['F'] = True

# run simulation:
for i in range(100):
    advance()
    
# get results
print("Person\tInfected  Immune  Dead")
for V in G.nodes():
    print("{}\t{}\t  {}\t  {}".format(V, infected[V], immune[V], dead[V]))

