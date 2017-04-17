#!/usr/bin/python3
"""Priority Queue based simulation. The priority queue stores the events that
   are due to occur. Creates a virus simulation modeling the spread through the
   population."""
# Keep sorted alphabetically.
import copy
import collections
import enum
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx
import numpy
import os
import pandas
import queue
import random
import urllib.request

SEED = 13378

# Generate plot on local machine rather than writing to file. For debugging.
SHOW_PLOT = False

SimulationResult = collections.namedtuple("SimulationResult", (
    "NodeCount", "ConnectedComponentCount", "MortalityRate", "TimeToDeath",
    "TimeToRecovery", "CommunicationTime", "ResistantProportion", "Time",
    "Infected", "Resistant", "Recovered", "Dead", "Communication", "Weighted"))


class Status(enum.Enum):
  """Enum for infection statuses."""
  UNEXPOSED = 1
  INFECTED = 2
  RESISTANT = 3
  RECOVERED = 4
  DEAD = 5


class NoNodesLeftAliveException(Exception):
  """Exception raised when there are no nodes left in the passed in graph."""
  pass


# Only need to inherit from object in Python2.
class VirusSimulation:
  def __init__(self, graph, mortality_rate=.2, time_to_death=25,
               time_to_recovery=25, communication_time=1,
               resistant_proportion=.1, max_iterations=100000,
               modify_graph=False, relationship_scores=None, weighted=False,
               plot_shape=False, plot_interval=25, plot_name=None):
    if weighted and relationship_scores is None:
      raise ValueError(
          "Cannot use weighted version without relationship scores.")
    self.graph = graph
    self.n_nodes = len(self.graph.nodes())
    self.relationship_scores = relationship_scores
    self.connected_components = networkx.number_connected_components(graph)
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
    self.weighted = weighted
    self.plot_shape = plot_shape(graph) if plot_shape else None
    self.plot_name = plot_name
    # At what interval to regenerate plots.
    self.plot_interval = plot_interval
    # Track what interval we are on.
    self.plot_counter = 0
    self.time = 0
    # Stats recording
    self.resistant = 0
    self.deaths = 0
    self.recoveries = 0
    self.infections = 0
    self.communications = 0
    if self.n_nodes == 0:
      raise NoNodesLeftAliveException()
    for node in self.graph.nodes():
      self.node_status[node] = Status.UNEXPOSED
    # We infect the first node.
    self.patient_zero = random.choice(self.graph.nodes())
    self.infect(self.patient_zero, True)
    self.plot_graph(force=self.plot_shape is not None)

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
      target_node = self.select_random_neighbor(node)
      if target_node is not None:
        if self.node_status.get(target_node) == Status.UNEXPOSED:
          self.infect(target_node)
        # Enqueue another communication from this infected node.
        communication_time = self.next_event(self.communication_time)
        self.pq.put((communication_time, self.communicate, node))
        self.communications += 1

  def select_random_neighbor(self, node):
    """Chooses which neighbor to communicate to next.

      Args:
        node: A node we need to chose a neighbor of.
      Returns:
        A neighboring node to target next. Returns None if there are no valid
        targets remaining.
      """
    if len(self.graph.neighbors(node)) == 0:
      return None
    if self.weighted:
      return self.weighted_random_neighbor(node)
    candidates = [
        neighbor for neighbor in self.graph.neighbors(node) if
        self.node_status.get(neighbor) != Status.DEAD]
    if candidates:
      return random.choice(candidates)
    return None

  def weighted_random_neighbor(self, node):
    """Choose who to meet based on our weighted criteria.

      Args:
        node: A node we need to chose a neighbor of.
      Returns:
        A neighboring node to target next. Returns None if there are no valid
        targets remaining.
      """
    candidates, weights = [], []
    neighbors = self.graph.neighbors(node)
    for neighbor in neighbors:
      if self.node_status.get(neighbor) != Status.DEAD:
        candidates.append(neighbor)
        weights.append(self.relationship_scores[node][neighbor])
    np_weights = numpy.array(weights)
    if candidates:
      return numpy.random.choice(
        candidates, 1, p=np_weights/np_weights.sum())[0]
    return None

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
    """Return selected input parameters."""
    return [self.n_nodes,
            self.connected_components,
            self.mortality_rate,
            self.time_to_death,
            self.time_to_recovery,
            self.communication_time,
            self.resistant_proportion]

  def run_virus(self):
    """Simulate virus until termination conditions reached."""
    while self.time < self.max_iterations and not self.pq.empty():
      timestep, function, node = self.pq.get()
      self.time = timestep
      function(node)
      self.plot_graph()
    print("Iteration final time = {}".format(self.time))
    print(("{} infections, {} resistant, {} recovered, {} dead "
           "{} communications").format(
        self.infections, self.resistant, self.recoveries, self.deaths,
        self.communications))
    self.plot_graph(force=self.plot_shape is not None)
    # This returns all the data regarding this simulation.
    return (*self.input_parameters(),
            self.time,
            self.infections,
            self.resistant,
            self.recoveries,
            self.deaths,
            self.communications,
            self.weighted)

  def plot_graph(self, force=False):
    """Decides wheter we need to generate plots.

    If necessary we generate them and increment our counter.
    """
    if force or (self.plot_shape and self.time >= self.plot_counter):
      node_status_color_map = {
        Status.INFECTED: "orange",
        Status.RESISTANT: "green",
        Status.RECOVERED: "cyan",
        Status.DEAD: "red",
        Status.UNEXPOSED: "blue"}
      node_colors = [node_status_color_map.get(node) for node in
                     self.node_status.values()]
      nodes = networkx.draw_networkx_nodes(self.graph, self.plot_shape,
                                           node_size=10, node_color=node_colors)
      nodes.set_edgecolor("k")
      # Plot patient zero node on top of all other nodes
      p0_color = node_status_color_map[self.node_status[self.patient_zero]]
      nodes = networkx.draw_networkx_nodes(self.graph, self.plot_shape,
                                           nodelist=[self.patient_zero],
                                           node_size=100, node_shape="*",
                                           node_color=p0_color)
      nodes.set_edgecolor("k")
      networkx.draw_networkx_edges(self.graph, self.plot_shape)
      legend_handles = []
      for key, color in node_status_color_map.items():
        legend_handles.append(mpatches.Patch(color=color, label=key))
      plt.legend(fontsize=8, handles=legend_handles)
      plt.tick_params(
          top='off', bottom='off', left='off', right='off', labelleft='off',
          labelbottom='off')
      plt.title("%s Graph @timestep %d" % (self.plot_name, self.time))
      plt.axis=("off")
      if SHOW_PLOT:
        plt.show()
      else:
        file_name = "%s%d.pdf" % (self.plot_name, self.time)
        print("Generating %s" % file_name)
        plt.savefig(file_name, format="pdf")
        plt.clf()
      # Increment counter.
      self.plot_counter += self.plot_interval


def relationship_score(graph, node, neighbor):
  """Compute the relationship score for two nodes in a graph.

  Args:
    node: A node from a networkx graph.
    neighbor: Another node from a networkx graph.

  Returns: (float) between 0 and 1."""
  a = set(graph.neighbors(node))
  b = set(graph.neighbors(neighbor))
  intersection = a.intersection(b)
  union = a.union(b)
  return (1 + len(intersection))/(1 + len(union))


def generate_relationship_scores(graph):
  """Compute the relaionship score for all pairs of linked nodes in a graph.

  Args:
    graph: (networkx.Graph) object.
  Returns:
    (dict) of dicts. Mapping each node to a dict of its neigbors mapped to
    float relationship scores.
  """
  relationship_score_map = {}
  for node in graph.nodes():
    relationship_score_map[node] = {}
    neighbors = graph.neighbors(node)
    for neighbor in neighbors:
      relationship_score_map[node].update({
          neighbor: relationship_score(graph, node, neighbor)})
  return relationship_score_map


def repeated_virus_simulation(graph, repetitions=30):
  """Repeatedly attack the same graph with multiple viruses.

  Each virus removes the dead nodes from the graph so the graph becomes sparser
  on subsequent iterations.

  Args:
    graph: (networkx.Graph) object.
    repetitions: (int) The number of generations to simulate for.
  Returns:
    (pandas.Dataframe) of the results.
  """
  # Prevent modification of external, passed in graph object.
  graph_copy = copy.deepcopy(graph)
  print("Simulating multiple sequential virus attacks on population")
  relationship_scores = generate_relationship_scores(graph_copy)
  results = []
  for i in range(repetitions):
    try:
      virus_simulation = VirusSimulation(
        graph_copy, relationship_scores=relationship_scores, modify_graph=True,
        weighted=True)
      results.append(SimulationResult(*virus_simulation.run_virus()))
    except NoNodesLeftAliveException:
      # If the graph nodes are empty from delete in modify_graph mode.
      continue
  dataframe = pandas.DataFrame(results)
  return dataframe


def weighting_system_analysis(graph, file_name, repetitions=30):
  relationship_scores = generate_relationship_scores(graph)
  results = []
  for i in range(repetitions):
    print("%s Complete" % (i/repetitions))
    weighted_virus_simulation = VirusSimulation(
        graph, relationship_scores=relationship_scores,
        weighted=True)
    unweighted_virus_simulation = VirusSimulation(
        graph, relationship_scores=relationship_scores,
        weighted=False)
    print("True"),
    results.append(SimulationResult(
        *weighted_virus_simulation.run_virus()))
    print("False"),
    results.append(SimulationResult(
        *unweighted_virus_simulation.run_virus()))
  dataframe = pandas.DataFrame(results)
  dataframe.to_csv(file_name)
  return dataframe


def get_facebook_graph():
  """Get real-world facebook data as a networkx graph.

  If the file is available locally we use that. Otherwise we fetch it online.

  Returns:
    A networkx Graph object.
  """
  file_name = "facebook_combined.txt.gz"
  if not os.path.exists(file_name):
    url = "https://snap.stanford.edu/data/%s" % file_name
    urllib.request.urlretrieve(url, file_name)
  return networkx.read_edgelist(file_name, nodetype=int)


# Program modes to run.
WEIGHTING_ANALYSIS = True
REPEATED_VIRUSES = False
PLOT_VIRUS = False

# Graph to use.
RANDOM_GRAPH = True
FACEBOOK_GRAPH = False

if __name__ == "__main__":

  # Generate our candidate graphs.
  facebook_graph = get_facebook_graph()
  nodes = len(facebook_graph.nodes())
  edges = len(facebook_graph.edges())
  edges_per_node = edges//nodes # Integer division.

  random_graph = networkx.barabasi_albert_graph(nodes, edges_per_node,
                                                seed=SEED)

  if WEIGHTING_ANALYSIS:
    if FACEBOOK_GRAPH:
      fb_dataframe = weighting_system_analysis(facebook_graph,
                                              "fb_weighting_analysis.csv", 10)
    if RANDOM_GRAPH:
      random_dataframe = weighting_system_analysis(
          random_graph, "random_weighting_analysis.csv", 10)

  if REPEATED_VIRUSES:
    if RANDOM_GRAPH:
      print("Random Data - Barabási Albert model:")
      result = repeated_virus_simulation(random_graph)
      print(result)

    if FACEBOOK_GRAPH:
      # Get real-world data
      print("Real-world Data - Facebook graph:")
      result = repeated_virus_simulation(facebook_graph)
      print(result)

  if PLOT_VIRUS:
    print("Generating plots - this may take some time")
    if RANDOM_GRAPH:
      print("Random Data - Barabási Albert model:")
      relationship_scores = generate_relationship_scores(random_graph)
      virus_simulation = VirusSimulation(
        random_graph, relationship_scores=relationship_scores,
        modify_graph=False, weighted=True, plot_name="Random",
        plot_shape=networkx.spring_layout)
      print(SimulationResult(*virus_simulation.run_virus()))

    if FACEBOOK_GRAPH:
      print("Real-world Data - Facebook graph:")
      relationship_scores = generate_relationship_scores(facebook_graph)
      virus_simulation = VirusSimulation(
        facebook_graph, relationship_scores=relationship_scores,
        modify_graph=False, weighted=True, plot_name="Facebook",
        plot_shape=networkx.spring_layout)
      print(SimulationResult(*virus_simulation.run_virus()))
