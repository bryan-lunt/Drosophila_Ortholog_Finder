#!/usr/bin/env python

import Bio.SeqIO as SIO
import Bio.pairwise2 as pw2

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--gapopen",default=-5.0,type=float)
parser.add_argument("--gapextend",default=-0.5,type=float)
parser.add_argument("IN1",type=str)
parser.add_argument("IN2",type=str)
args, other = parser.parse_known_args()

in_1 = list(SIO.parse(args.IN1,"fasta"))
in_2 = list(SIO.parse(args.IN2,"fasta"))

GAPOPEN = args.gapopen
if GAPOPEN > 0:
	GAPOPEN = -GAPOPEN

GAPEXTEND = args.gapextend
if GAPEXTEND > 0:
	GAPEXTEND = -GAPEXTEND

aligns = list()

for left, right in zip(in_1,in_2):
	forwardalign = pw2.align.globalxs(left.seq.upper(),right.seq.upper(),GAPOPEN,GAPEXTEND)[0]
	revcompalign = pw2.align.globalxs(left.seq.upper(),right.seq.reverse_complement().upper(),GAPOPEN,GAPEXTEND)[0]	
	
	aligns.append((forwardalign, revcompalign))

