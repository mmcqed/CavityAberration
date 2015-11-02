#!/bin/bash

rm /farmshare/user_data/papag/StandardErrorFiles/*
rm /farmshare/user_data/papag/StandardOutputFiles/*

rm /afs/ir.stanford.edu/users/p/a/papag/Documents/LevLab/CavityAberrationSimulation/CavDataDir/*
rm /afs/ir.stanford.edu/users/p/a/papag/Documents/LevLab/CavityAberrationSimulation/CavityModeFunctions/*

python CavityModeMaker.py
qsub -t 1-121 qsubArray
