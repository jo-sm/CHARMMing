#!/usr/bin/env python

import glob
import lib.Etc as Etc

def mean(iterable):
    return sum(iterable)/float(len(iterable))

def stdDev(iterable):
    avg = mean(iterable)
    devSq = ( (x - avg)**2 for x in iterable )
    return (sum(devSq)/float(len(iterable)))**.5

def min(iterable):
    iterable = map(float,iterable)
    min = mean(iterable)
    for taco in iterable:
        if taco < min: min = taco
    return min

def max(iterable):
    iterable = map(float,iterable)
    max = mean(iterable)
    for taco in iterable:
        if taco > max: max = taco
    return max

def load_nat(file):
    return [ float(line) for line in open(file) ]

def load_all():
    all_nat = glob.glob('*.nat')
    all_nat.sort()
    return Etc.flatten([ load_nat(nat) for nat in all_nat ])




