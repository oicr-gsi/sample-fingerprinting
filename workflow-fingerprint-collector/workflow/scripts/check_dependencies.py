#!/usr/bin/env python
from argparse import ArgumentParser
import os
#from evac.varcallers import run_caller, list_callers

# Main

def main(script_dir):
    parser = ArgumentParser()
    parser.add_argument(
        '-b', '--bam',
        default='-', metavar="PATH",
        help="Location of the BAM file")
    parser.add_argument(
        '-o', '--output',
        default="-", metavar="PATH",
        help="Directory to output.")
    parser.add_argument(
        '-c', '--caller',
        choices=list_callers(), default='mpileup',
        help="The caller pipeline to use")
