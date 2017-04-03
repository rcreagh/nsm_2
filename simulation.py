#! /usr/bin/python3
"""This script simulates the spread of a computer virus across the machines of
employees in a company.
If a node is in state 0, it is not infected, but vulnerable.
If a node is in state 1, it is infected.
If a node is in state -1, it is immune to infection.
"""

import networkx as nx
import random

N_EMPLOYEES = 5000
EDGE_EXISTENCE_PROBABILITY = 0.1
PROBABILITY_OF_INFECTION = 0.5
N_TIMESTEPS = 10
WINDOWS_OS_USER_PROPORTION = 0.95
N_INIT_INFECTED_USERS = 1
PREVIOUSLY_INFECTED_USERS = set()

def infection_init(G):
    """Make a graph with one initial infected node.

    Args:
      G: (networkx.Graph) object
    Returns:
      None, modifies G in place.
    """
    for u in G.nodes():
        G.node[u]["state"] = 0
        rand = random.random()
        if rand <= WINDOWS_OS_USER_PROPORTION:
            G.node[u]["OS"] = "Windows"
        else:
            G.node[u]["OS"] = "Linux"

    init = random.sample(G.nodes(), N_INIT_INFECTED_USERS)
    for u in init:
        G.node[u]["state"] = 1
        G.node[u]["OS"] = "Windows"


def update_infection_status(G):
    """Update infection status of each node.
    Args:
      G: (networkx.Graph) object
    Returns:
      None, modifies G in place.
    """
    newly_infected_nodes = {}
    for u in G.nodes():
        if G.node[u]["state"] == 1:
            if u in PREVIOUSLY_INFECTED_USERS:
                continue
            PREVIOUSLY_INFECTED_USERS.add(u)
            for n in G.neighbors(u):
                if G.node[n]["state"] == 0:
                    rand = random.random()
                    if (G.node[u]["OS"] == "Windows" and
                        rand <= PROBABILITY_OF_INFECTION):
                        newly_infected_nodes[n] = 1
                    else:
                        G.node[n]["state"] = -1
    for key, value in newly_infected_nodes.items():
        user = key
        G.node[user]["state"] = 1


def main():
    G = nx.erdos_renyi_graph(N_EMPLOYEES, EDGE_EXISTENCE_PROBABILITY)
    infection_init(G)
    for i in range(N_TIMESTEPS):
        update_infection_status(G)
        number_of_infected_users = sum(
                G.node[i]["state"] > 0 for i in G.nodes())
        print("Number of infected users at step %d: %d" %
              (i, number_of_infected_users))


if __name__ == "__main__":
    main()
