#!/usr/bin/python3
"""Priority Queue based simulation. The priority queue stores the events that
   are due to occur."""
import queue
import networkx
import random
import numpy
import enum
import collections
import pandas


# For the random graph
N = 1000
M = 30
SEED = 13378


Result = collections.namedtuple("SimulationResults", (
    "MortalityRate", "TimeToDeath", "TimeToRecovery", "EmailTime",
    "ResistantProportion", "Time", "Infected", "Resistant", "Recovered",
    "Dead", "Emails"))


class Status(enum.Enum):
  INFECTED = 1
  RESISTANT = 2
  RECOVERED = 3
  DEAD = 4


# Only need to inherit from object in Python2.
class VirusSimulation:
  def __init__(self, n_nodes, m_edges, seed, mortality_rate=.2,
               time_to_death=25, time_to_recovery=25, email_time=1,
               resistant_proportion=.1):
    self.graph = networkx.barabasi_albert_graph(n_nodes, m_edges, seed=seed)
    self.pq = queue.PriorityQueue()
    self.mortality_rate = mortality_rate
    self.time_to_death = time_to_death
    self.time_to_recovery = time_to_recovery
    self.email_time = email_time
    self.resistant_proportion = resistant_proportion
    self.time = 0
    # Stats recording
    self.resistant = 0
    self.deaths = 0
    self.recoveries = 0
    self.infections = 0
    self.emails = 0
    # We infect the first node, node 0.
    self.infect(0, True)

  def input_parameters(self):
    return [self.mortality_rate,
            self.time_to_death,
            self.time_to_recovery,
            self.email_time,
            self.resistant_proportion]

  def next_event(self, exponential_mean):
    return self.time + numpy.random.exponential(exponential_mean)

  def infect(self, node, patient_zero=False):
    if not patient_zero and random.random() < self.resistant_proportion:
      self.graph.node[node]["status"] = Status.RESISTANT
      self.resistant += 1
    else:
      self.graph.node[node]["status"] = Status.INFECTED
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
      email_time = self.next_event(self.email_time)
      self.pq.put((email_time, self.mail, node))
      self.infections += 1

  def mail(self, node):
    if self.graph.node[node]["status"] == Status.INFECTED:
      # Select a random neighbor and infect them if possible.
      target_node = random.choice(self.graph.neighbors(node))
      if self.graph.node[target_node].get("status") is None:
        self.infect(target_node)
      # Enqueue another communication from this infected node.
      email_time = self.next_event(self.email_time)
      self.pq.put((email_time, self.mail, node))
      self.emails += 1

  def recover(self, node):
    if self.graph.node[node]["status"] == Status.INFECTED:
      self.graph.node[node]["status"] = Status.RECOVERED
      self.recoveries += 1

  def death(self, node):
    if self.graph.node[node]["status"] == Status.INFECTED:
      self.graph.node[node]["status"] = Status.DEAD
      self.deaths += 1

  def run_virus(self, max_iterations=100000):
    while self.time < max_iterations and not self.pq.empty():
      timestep, function, node = self.pq.get()
      self.time = timestep
      function(node)
    print("Final time = %s" % (self.time))
    print("%d infections, %d resistant, %d recovered, %d dead %d emails" % (
        self.infections, self.resistant, self.recoveries, self.deaths,
        self.emails))
    # This returns all the data regarding this simulation.
    return (*self.input_parameters(),
            self.time,
            self.infections,
            self.resistant,
            self.recoveries,
            self.deaths,
            self.emails)


if __name__ == "__main__":
  results = []
  for i in range(30):
    virus_simulation = VirusSimulation(N, M, SEED)
    results.append(Result(*virus_simulation.run_virus()))
  dataframe = pandas.DataFrame(results)
  print(dataframe)
