#!/usr/bin/env python3

'''
title: SwipProteomics
description: Leidenalg Clustering of the Enhanced Protein Covariation Network
author: Tyler W A Bradshaw

This script does the work of clustering a network provided as an adjacency
matrix using the Leiden algorithm and the Python Leidenalg library. 

Leidenalg is the work of Vincent Traag:
http://dx.doi.org/10.1038/s41598-019-41695-z

'''


## ---- Input

root = "~/projects/SwipProteomics"

# There is only one piece of input data, an input NxN adjacency
# matrix saved as a csv file in  root/rdata/.
adjm_file = 'ne_adjm.csv'

# Optimization methods:
optimization_method = recursive_method = 'Surprise'

# If module_size > max_size, then cluster recursively.
recursive = False

# Prefix out output partition, saved as .csv.
output_name = 'swip'

# Parameters for multiresolution methods:
rmin = 1 # Min resolution for multi-resolution methods.
rmax = 1 # Max resolution for multi-resolution methods.
nsteps = 1 # Number of steps to take between rmin and rmax.
max_size = 100 # Maximum allowable size of a module.

## Optimization parameters
# Not the number of recursive iterations, but the number
# of optimization iterations.
n_iterations = -1

## Output
# Saved in root/rdata/
# [output_name]_partitions.csv
# * a partition of the network saved in root/rdata

## NOTE: some params may not be used if not required
## by optimization_method

import numpy as np
from igraph import Graph
from pandas import DataFrame


## ---- Prepare the workspace

import os
import sys
import glob
from sys import stderr
from os.path import dirname

import numpy as np
from numpy import linspace
from importlib import import_module
from progressbar import ProgressBar
from leidenalg import Optimiser, find_partition

from igraph import Graph
from pandas import read_csv, DataFrame

# project Directories:
rdatdir = os.path.join(root,"rdata")
funcdir = os.path.join(root,"Py")


## ---- Functions

def graph_from_adjm(adjm,weighted=True,signed=True):
    if not signed: adjm = abs(adjm)
    edges = adjm.stack().reset_index()
    edges.columns = ['nodeA','nodeB','weight']
    edges = edges[edges.weight != 0]
    edge_tuples = list(zip(edges.nodeA,edges.nodeB,edges.weight))
    if weighted: g = Graph.TupleList(edge_tuples,weights=True)
    if not weighted: g = Graph.TupleList(edge_tuples,weights=False)
    return g
#EOF


## ---- Leidenalg qualitity metrics

# Leidenalg supports the following optimization methods:
methods = {
        # Modularity
        "Modularity": {'partition_type' : 'ModularityVertexPartition',
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations,
            'multi_resolution' : False },
        # Surprise
        "Surprise": {'partition_type' : 'SurpriseVertexPartition',
            'weights' : True, 'signed' : False,
            'resolution_parameter' : None, 'n_iterations' : n_iterations,
            'multi_resolution' : False },
        # RBConfiguration
        "RBConfiguration": {'partition_type' : 'RBConfigurationVertexPartition',
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations,
            'multi_resolution' : True },
        # RBER
        "RBER": {'partition_type' : 'RBERVertexPartition',
            'weights' : True, 'signed' : False,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations,
            'multi_resolution' : True },
        # CPM
        "CPM": {'partition_type' : 'CPMVertexPartition',
            'weights' : True, 'signed' : True,
            'resolution_parameter' : {'start':rmin,'stop':rmax,'num':nsteps},
            'n_iterations' : n_iterations,
            'multi_resolution' : True },
        # Significance
        # FIXME: Significance method doesn't seem to be working.
        "Significance":
        {'partition_type' : 'SignificanceVertexPartition',
            'weights': None, 'signed' : False,
            'resolution_parameter' : None,
            'n_iterations' : n_iterations,
            'multi_resolution' : False }
        }

# Get method specific parameters for clustering.
parameters = methods.get(optimization_method)
method = parameters.get('partition_type')

# Status report.
print("Performing Leidenalg clustering utilizing the {}".format(method),
        "method to find optimal partition(s).", file=stderr)

if recursive:
    msg = "Recursively splitting modules larger than {} nodes with '{}'."
    print(msg.format(max_size,recursive_method))
#EIS

if parameters.get('multi_resolution') is True:
    msg = "Analyzing graph at {} resolutions."
    R = len(linspace(**parameters.get('resolution_parameter')))
    print(msg.format(R), file=stderr)
#EIS



## ---- Load input adjacency matrix and create an igraph object

# Load graph adjacency matrix.
myfile = os.path.join(rdatdir,adjm_file)
adjm = read_csv(myfile, header = 0, index_col = 0)

# Create igraph graph object and add to parameters dictionary.
# Note, this can take several minutes.
if parameters.get('weights') is True:
    # Create a weighted graph.
    g = graph_from_adjm(adjm,weighted=True,signed=parameters.pop('signed'))
    parameters['weights'] = 'weight'
    parameters['graph'] = g
else:
    # Create an unweighted graph.
    g = graph_from_adjm(adjm,weighted=False,signed=parameters.pop('signed'))
    parameters['graph'] = g
#EIS

# the input graph
print("Input graph: {}".format(g.summary()),file=stderr)

# NOTE: UNW is UNdirected and Weighted graph!


## ---- Community detection with the Leiden algorithm

# FIXME: really we need two functions:
# * single-resolution (+/- recursive)
# * multi-resolution methods

# Update partition type parameter.
# Dynamically load the partition_type class.
# This is the method to be used to optimize the clustering.
parameters['partition_type'] = getattr(import_module('leidenalg'), method)

# Remove any None type parameters.
out = [key for key in parameters if parameters.get(key) is None]
for key in out: del parameters[key]


## ---- Perform Leidenalg module detection

## MULTIRESOLUTION METHODS

if parameters.pop('multi_resolution') is True:
  # collect clustering params
  profile = list()
  res_range = linspace(**parameters.get('resolution_parameter'))
  pbar = ProgressBar()
  #n_iter = parameters.pop(n_iterations)
  # for res in pbar(linespace(**p.pop('resolution_range'))):
  for resolution in pbar(res_range):
        # optimize partition
        parameters['resolution_parameter'] = resolution
        partition = find_partition(**parameters)
        optimiser = Optimiser()
        diff = optimiser.optimise_partition(partition,parameters['n_iterations'])
        # if recursive, split modules that are larger than min size
        if recursive:
            # update params
            print("Initial partition: {}".format(partition.summary()),
                    file=stderr)
            initial_membership = partition.membership
            new_params = parameters.copy()
            new_method = methods.get(recursive_method).get('partition_type')
            if type(new_method) is str:
              new_params['partition_type'] = getattr(import_module('leidenalg'), new_method)
            else:
              new_params['partition_type'] = new_method
            if 'resolution_parameter' in new_params.keys():
              new_params.pop('resolution_parameter')
            # get modules to be split
            subgraphs = partition.subgraphs()
            too_big = [subg.vcount() > max_size for subg in subgraphs]
            n_big = sum(too_big)
            # recursively split modules that are too big
            while any(too_big):
                idx = [i for i, too_big in enumerate(too_big) if too_big]
                new_graph = subgraphs.pop(idx[0])
                # mask negative edges
                edge_weights = new_graph.es['weight']
                new_weights = [e if e > 0 else 0 for e in edge_weights]
                new_graph.es['weight'] = new_weights
                new_params['graph'] = new_graph
                # find optimal partition of subgraph
                part = find_partition(**new_params)
                optimiser = Optimiser()
                diff = optimiser.optimise_partition(part,n_iterations=-1)
                subgraphs.extend(part.subgraphs())
                too_big = [subg.vcount() > max_size for subg in subgraphs]
            #EOL to split modules
            # Collect subgraph membership as a single partition.
            nodes = [subg.vs['name'] for subg in subgraphs]
            parts = [dict(zip(n,[i]*len(n))) for i, n in enumerate(nodes)]
            new_part = {k: v for d in parts for k, v in d.items()}
            # Set membership as partition after recursive splitting
            membership = [new_part.get(node) for node in partition.graph.vs['name']]
            partition.set_membership(membership)
        # EIS if recursive
        profile.append(partition)
    # EOL through resolutions
else:
    ## SINGLE RESOLUTION METHODS
    # Single resolution methods: first iteration if recursive
    profile = list()
    partition = find_partition(**parameters)
    optimiser = Optimiser()
    diff = optimiser.optimise_partition(partition,n_iterations=-1)
    profile.append(partition)
    if recursive:
        # update clustering params
        initial_membership = partition.membership
        new_params = parameters.copy()
        new_method = methods.get(recursive_method).get('partition_type')
        if type(new_method) is str:
            new_params['partition_type'] = getattr(import_module('leidenalg'), new_method)
        else:
            new_params['partition_type'] = new_method
        if 'resolution_parameter' in new_params.keys():
            new_params.pop('resolution_parameter')
        subgraphs = partition.subgraphs()
        too_big = [subg.vcount() > max_size for subg in subgraphs]
        n_big = sum(too_big)
        while any(too_big):
            # Perform clustering for any subgraphs that are too big.
            idx = [i for i, too_big in enumerate(too_big) if too_big]
            new_graph = subgraphs.pop(idx[0])
            # mask negative edges
            edge_weights = new_graph.es['weight']
            new_weights = [e if e > 0 else 0 for e in edge_weights]
            new_graph.es['weight'] = new_weights
            new_params['graph'] = new_graph
            # find optimal partition of subgraph
            part = find_partition(**new_params)
            optimiser = Optimiser()
            diff = optimiser.optimise_partition(part,n_iterations=-1)
            subgraphs.extend(part.subgraphs())
            too_big = [subg.vcount() > max_size for subg in subgraphs]
        #EOL to split modules
        # Collect subgraph membership as a single partition.
        nodes = [subg.vs['name'] for subg in subgraphs]
        parts = [dict(zip(n,[i]*len(n))) for i, n in enumerate(nodes)]
        new_part = {k: v for d in parts for k, v in d.items()}
        # Set membership of initial graph.
        membership = [new_part.get(node) for node in partition.graph.vs['name']]
        partition.set_membership(membership)
        profile[0] = partition
    # EIS if recursive
#EIS if multi_resolution else single resolution

print("Final partition:  {}".format(partition.summary()), file=stderr)


## ---- Save Leidenalg clustering results

# collect matrix in which each row is a partition (nrow = nRes)
df = DataFrame(columns = profile[0].graph.vs['name'])
for i in range(len(profile)):
    df.loc[i] = profile[i].membership
#EOL

# save the data as csv
myfile = os.path.join(rdatdir, output_name + "_partition.csv")
df.to_csv(myfile)
print("saved: {}".format(myfile), file=stderr)
