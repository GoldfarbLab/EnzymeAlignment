#!/usr/bin/env python
import sys
import os
import argparse


########################################################################################################################
# Command line argument parsing
########################################################################################################################

description = "This is the entry point to the program. It will execute the requested tasks."

# initialize the parser
parser = argparse.ArgumentParser(description=description)


# generic arguments
parser.add_argument("task",
                    help="execute the specified task",
                    choices=["pipeline", "align-global", "align-local", "MSA-opt",  "HMM-make", "HMM-search",
                             "thread", "MSA", "annotate"])

parser.add_argument("-n", "--num-jobs",
                    help="number of jobs if running on a cluster",
                    type=int)

parser.add_argument("-c", "--compute-environment",
                    help="choose execution environment",
                    choices=["local", "LSF", "PBS"],
                    default="local")

parser.add_argument("-d", "--database",
                    help="path to sequence database")

parser.add_argument("-D", "--database-name",
                    help="nice name for sequence database")

parser.add_argument("-e", "--enzyme",
                    help="path to enzyme data")

parser.add_argument("-E", "--enzyme-name",
                    help="nice name for enzyme")


# output arguments
parser.add_argument("-o", "--output-path")


# global pairwise alignment arguments
parser.add_argument("--align-global-identity-threshold",
                    help="initial percent identity threshold for global pairwise alignment",
                    type=float,
                    default=0.2)


# local pairwise alignment arguments


# MSA optimization arguments


# HMM creation arguments


# HMM search arguments


# threading arguments


# annotation arguments


args = parser.parse_args()

if args.task:
    print(args.task)

else:
    print('fail')

