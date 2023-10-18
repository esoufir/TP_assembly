#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open(fastq_file, "r+") as f_in:
        for line in f_in:
            # line is already the first line
            yield next(f_in).strip()
            # Skipping the +
            next(f_in)
            # Skipping the quality
            next(f_in)


def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    # Sliding window : 
    for i in range(len(read)-kmer_size+1):
        yield read[i:i+kmer_size]
    


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    dict_kmer = dict()
    for seq in read_fastq(fastq_file):
        for km in cut_kmer(seq, kmer_size):
            dict_kmer[km]=0
        for km in cut_kmer(seq, kmer_size):
            dict_kmer[km]+=1
    return dict_kmer


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graph = nx.DiGraph()
    for k, v in kmer_dict.items():
        graph.add_edge(k[:len(k)-1], k[1:len(k)], weight = v)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    
    for path in path_list:
        # Tous les nœuds d’un chemin sont supprimés.
        if delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)

        # Tous les nœuds d’un chemin sont supprimés sauf le premier et le dernier.
        elif not delete_entry_node and not delete_sink_node:
            graph.remove_nodes_from(path[1:len(path)-1]) # -1 ?
        # Tous les nœuds d’un chemin sont supprimés sauf le dernier nœud.
        elif delete_entry_node:
            graph.remove_nodes_from(path[:len(path)-1])
        elif delete_sink_node:
            graph.remove_nodes_from(path[1:])
    return graph


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    sdv_weights = statistics.stdev(weight_avg_list)

    indice_max_path = 0
    max_val = 0
    if sdv_weights>0:
        for i in range(len(weight_avg_list)):
            if weight_avg_list[i]>max_val:
                indice_max_path=i
                max_val = weight_avg_list[i]
    else:
        sdv_length = statistics.stdev(path_length)
        if sdv_length > 0:
            #Chemin le plus long : 
            for i in range(len(path_length)):
                if path_length[i]>max_val:
                    indice_max_path=i
                    max_val = path_length[i]
        else : 
            indice_max_path = random.randint(0, len(path_length))
    print("INDICE MAX", indice_max_path)
    path_to_remove = path_list[:indice_max_path] + path_list[indice_max_path+1:]
    
    print("TO RM",path_list[indice_max_path])
    # Filtering the graph:
    graph_cleaned = remove_paths(graph, path_to_remove, delete_entry_node, delete_sink_node)
    print(graph_cleaned.edges())
    return graph_cleaned
            
    
def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    # Path between the two nodes :
    path_length = []
    path_weights = []
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    for path in path_list:
        path_length.append(len(path))
        # Poids des arretes : 
        path_weights.append(path_average_weight(graph, path))
    # Cleaning :
    
    graph = select_best_path(graph, path_list, path_length, path_weights, delete_entry_node=False, delete_sink_node=False)
    return graph

def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    pass

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    list_entry_nodes = []
    nodes = list(graph.nodes)
    for node in nodes:
        if len(list(graph.predecessors(node))) == 0:
            list_entry_nodes.append(node)
    return list_entry_nodes

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    list_exit_nodes = []
    nodes = list(graph.nodes)
    for node in nodes:
        if len(list(graph.successors(node)))== 0:
            list_exit_nodes.append(node)
    return list_exit_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            # If there is a pathway between the two nodes = contig
            if nx.has_path(graph, starting_node, ending_node):
                for path in nx.all_simple_paths(graph, starting_node, ending_node): 
                    contig = path[0]
                    for node in path[1:]:
                        contig+=node[-1]
                    contigs.append((contig,len(contig)))
    return contigs

def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, "w") as f_out:
        for i in range(len(contigs_list)):
            seq,longueur = contigs_list[i]
            print(seq, longueur)
            # Writing the header:
            f_out.write(f">contig_{i} len={longueur}\n")
            f_out.write(f"{textwrap.fill(seq, width=80)}\n")



def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    dict_km = build_kmer_dict(args.fastq_file, 7)
    graph = (build_graph(dict_km))
    starting_nodes = get_starting_nodes(graph)

    ending_nodes = get_sink_nodes(graph)
    contigs = get_contigs(graph, starting_nodes, ending_nodes)
    # Saving the contigs :
    save_contigs(contigs, "test_output.fasta")
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()
