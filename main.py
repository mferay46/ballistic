# -*- coding: utf-8 -*-
"""
Main
@author: FERAY.M
"""
from model import AnalyticSolution
from model import NumericalSolution

from pop import start_popup
from params import *

task, mode = start_popup()
def main():
    if task == 0:
        model = AnalyticSolution(params)
        model.run()
        model.plot_all()
    if task == 1:
        model = NumericalSolution(params)
        model.run_all(mode) 
    if task == 2:
        NumericalSolution.plot_trajectories(params)   
    if task == 3:
        NumericalSolution.isocontours(params)
main()