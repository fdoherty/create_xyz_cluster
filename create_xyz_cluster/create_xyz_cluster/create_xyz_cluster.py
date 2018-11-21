#!/usr/bin/env python3
"""
create_xyz_cluster.py
Script used for creating a text file with atom coordinates:
 - Reads file output from genetic algorithm for supported nanoclusters
 - Locates minimum energy structure from list of all structures
 - Uses relative atom positions to create xyz list of cluster geometry
 - Writes to .txt file with proper syntax for Atomic Simulation Environment
"""

from __future__ import print_function
import sys
from argparse import ArgumentParser
import numpy as np
import os

SUCCESS = 0
IO_ERROR = 1

DEFAULT_INPUT_FILE_NAME = 'minima.xyz'


def warning(*objs):
    """Writes a message to stderr."""
    print("WARNING: ", *objs, file=sys.stderr)


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = ArgumentParser(description='Reads in minima.xyz file from genetic algorithm output and creates '
                                        'a separate .txt file with coordinates of most stable configuration')
    parser.add_argument("-m", "--minima_file", help="The location (directory and file name) of the minima.xyz file",
                        default=DEFAULT_INPUT_FILE_NAME)
    args = None
    try:
        args = parser.parse_args(argv)
        with open(args.minima_file, "r") as f:
            args.minima_data = f.read().splitlines()
        args.minima_data = [i.split() for i in args.minima_data]
        # builds data into list of lists to handle variable row lengths
    except IOError as e:
        warning("Problems reading file:", e)
        parser.print_help()
        return args, IO_ERROR

    return args, SUCCESS


def minima_stats(listoflists):
    """
    Finds the minimum energy structure from the input data and returns the
    cartesian coordinates of the cluster geometry along with the cluster size.

    Parameters
    ----------
    listoflists : List of lists containing the delimited data from minima.xyz

    Returns
    -------
    npsize : Size (number of atoms) of the metal cluster
    output : Full list of relative cartesian coordinates for the cluster
    position : Cartesian coordinates of the first atom in the cluster
    """
    # Load first column to find energy data and number of total atoms
    first = []
    for i in listoflists:
        try:
            first.append(float(i[0]))
        except:
            first.append(i[0])

    # Find index of minimum energy
    index_min = np.argmin(first)

    # Find number of total atoms
    natoms = int(first[0])

    # Create list of elements within structure
    elem = [i[0] for i in listoflists[index_min + 1: index_min + natoms + 1]]

    # Build x,y,z coordinates for each atom in structure
    x = [float(i[1]) for i in listoflists[index_min + 1: index_min + natoms + 1]]
    y = [float(i[2]) for i in listoflists[index_min + 1: index_min + natoms + 1]]
    z = [float(i[3]) for i in listoflists[index_min + 1: index_min + natoms + 1]]

    # Count number of catalyst metal atoms in structure
    npsize = sum(i == 'Rh' for i in elem)

    # Build full list of relative cartesian coordinates for the cluster
    output = [(0, 0, 0)]
    output.extend((x[-npsize + i] - x[-npsize],
                   y[-npsize + i] - y[-npsize],
                   z[-npsize + i] - z[-npsize]) for i in range(1, npsize))

    # Store cartesian coordinates of the first atom in the cluster, relative to
    # the metal oxide slab
    position = [x[-npsize] - 5.98528059,
                y[-npsize] - 3.78923,
                z[-npsize] - max(z[:-npsize])]

    return npsize, output, position


def main(argv=None):
    args, ret = parse_cmdline(argv)
    if ret != SUCCESS:
        return ret
    npsize, output, position = minima_stats(args.minima_data)

    # Store cluster size and position in syntax required by Atomic Simulation
    # Environment
    string = """cluster = Atoms('Rh%d', %s)
    add_adsorbate(slab, cluster, %f, position=(%f,%f))""" % (npsize, output,
                                                             position[2],
                                                             position[0],
                                                             position[1])

    # Write new .txt file with Atomic Simulation Environment input
    f = open("coords.txt", "w+")
    f.write(string)
    f.close()
    return SUCCESS


if __name__ == "__main__":
    status = main()
    sys.exit(status)
