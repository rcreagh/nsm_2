#!/usr/bin/python3
"""Priority Queue based simulation. The priority queue stores the events that
   are due to occur. Creates a virus simulation modeling the spread through the
   population."""
# Keep sorted alphabetically.
import collections
import enum
import networkx
import numpy
import pandas
import queue
import random


# For the random graph
N = 1000
M = 30
SEED = 13378


SimulationResult = collections.namedtuple("SimulationResult", (
    "MortalityRate", "TimeToDeath", "TimeToRecovery", "CommunicationTime",
    "ResistantProportion", "Time", "Infected", "Resistant", "Recovered",
    "Dead", "Communication"))


class Status(enum.Enum):
  INFECTED = 1
  RESISTANT = 2
  RECOVERED = 3
  DEAD = 4


class NoNodesLeftAliveException(Exception):
  """Exception raised when there are no nodes left in the passed in graph."""
  pass


# Only need to inherit from object in Python2.
class VirusSimulation:
  def __init__(self, graph, mortality_rate=.2, time_to_death=25,
               time_to_recovery=25, communication_time=1,
               resistant_proportion=.1, max_iterations=100000,
               modify_graph=False):
    self.graph = graph
    # Stores infection status of the nodes. Stored here to prevent modification
    # of passed in graph. Stops us from having to reinitalize or copy graphs.
    self.node_status = {}
    self.pq = queue.PriorityQueue()
    self.mortality_rate = mortality_rate
    self.time_to_death = time_to_death
    self.time_to_recovery = time_to_recovery
    self.communication_time = communication_time
    self.resistant_proportion = resistant_proportion
    self.max_iterations = max_iterations
    self.modify_graph = modify_graph
    self.time = 0
    # Stats recording
    self.resistant = 0
    self.deaths = 0
    self.recoveries = 0
    self.infections = 0
    self.communications = 0
    if len(self.graph.nodes()) == 0:
      raise NoNodesLeftAliveException()
    # We infect the first node.
    self.infect(random.choice(self.graph.nodes())
, True)

  def next_event(self, exponential_mean):
    return self.time + numpy.random.exponential(exponential_mean)

  def infect(self, node, patient_zero=False):
    if not patient_zero and random.random() < self.resistant_proportion:
      self.node_status[node] = Status.RESISTANT
      self.resistant += 1
    else:
      self.node_status[node] = Status.INFECTED
      # Enqueue communication
      if random.random() < self.mortality_rate:
        # Node will die
        death_time = self.next_event(self.time_to_death)
        self.pq.put((death_time, self.death, node))
      else:
        # Node will recover
        recovery_time = self.next_event(self.time_to_recovery)
        self.pq.put((recovery_time, self.recover, node))
      # Regardless we need to enque communications
      communication_time = self.next_event(self.communication_time)
      self.pq.put((communication_time, self.communicate, node))
      self.infections += 1

  def communicate(self, node):
    if self.node_status[node] == Status.INFECTED:
      # Select a random neighbor if you have any and infect them if possible.
      if len(self.graph.neighbors(node)) > 0:
        target_node = random.choice(self.graph.neighbors(node))
        if self.node_status.get(target_node) is None:
          self.infect(target_node)
      # Enqueue another communication from this infected node.
      communication_time = self.next_event(self.communication_time)
      self.pq.put((communication_time, self.communicate, node))
      self.communications += 1

  def recover(self, node):
    if self.node_status[node] == Status.INFECTED:
      self.node_status[node] = Status.RECOVERED
      self.recoveries += 1

  def death(self, node):
    if self.node_status[node] == Status.INFECTED:
      self.node_status[node] = Status.DEAD
      self.deaths += 1
      if self.modify_graph:
        self.graph.remove_node(node)

  def input_parameters(self):
    return [self.mortality_rate,
            self.time_to_death,
            self.time_to_recovery,
            self.communication_time,
            self.resistant_proportion]

  def run_virus(self):
    while self.time < self.max_iterations and not self.pq.empty():
      timestep, function, node = self.pq.get()
      self.time = timestep
      function(node)
    print("Final time = %s" % (self.time))
    print(("%d infections, %d resistant, %d recovered, %d dead "
           "%d communications") % (
        self.infections, self.resistant, self.recoveries, self.deaths,
        self.communications))
    # This returns all the data regarding this simulation.
    return (*self.input_parameters(),
            self.time,
            self.infections,
            self.resistant,
            self.recoveries,
            self.deaths,
            self.communications)


if __name__ == "__main__":
  graph = networkx.barabasi_albert_graph(N, M, seed=SEED)
  results = []
  for i in range(30):
    try:
      virus_simulation = VirusSimulation(graph, modify_graph=True)
      results.append(SimulationResult(*virus_simulation.run_virus()))
    except NoNodesLeftAliveException:
      # If the graph nodes are empty from delete in modify_graph mode.
      continue
  dataframe = pandas.DataFrame(results)
  print(dataframe)
