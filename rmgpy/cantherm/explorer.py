#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import os
import mpmath as mp
import numpy as np

from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.pdep import PDepNetwork
from rmgpy.molecule.molecule import Molecule

class ExplorerJob(object):
    def __init__(self, source, pdepjob, explore_tol, energy_tol, flux_tol, ):
        self.source = source
        self.explore_tol = explore_tol
        self.energy_tol = energy_tol
        self.flux_tol = flux_tol
        
        self.pdepjob = pdepjob
        self.network = pdepjob.network
        self.network.__class__ = PDepNetwork
        self.network.source = source
        self.network.explored = []
        self.network.index = -1

    def copy(self):
        """
        Return a copy of the pressure dependence job.
        """
        return ExplorerJob(
            source=self.source,
            pdepjob=self.pdepjob,
            explore_tol=self.explore_tol,
            energy_tol=self.energy_tol,
            flux_tol=self.flux_tol
        )
        
    def execute(self):
        
        reactionModel = CoreEdgeReactionModel()
        
        bathgases = []
        for g in self.network.bathGas.keys():
            if not g in self.source:
                g.reactive = False
                bathgases.append(g)
        
        for spc in bathgases+self.source: #add species to model
            spc,isNew = reactionModel.makeNewSpecies(spc)
            reactionModel.enlarge(spc,reactEdge=False,unimolecularReact=np.ones(len(reactionModel.core.species)),
                      bimolecularReact=np.zeros((len(reactionModel.core.species),len(reactionModel.core.species))))
        
        #react initial species
        reactionModel.enlarge(reactEdge=True,unimolecularReact=np.ones(len(reactionModel.core.species)),
                      bimolecularReact=np.zeros((len(reactionModel.core.species),len(reactionModel.core.species))))
        
        #find the network we're interested in, merge it with ours and put it back in
        for i,nwk in enumerate(reactionModel.networkList):
            if nwk.source == self.source:
                network = nwk
                ind = i
                break
        else:
            raise ValueError('did not generate a network with the requested source')
        
        index = network.index
        
        network = self.network.merge(network)
        network.index = index
        
        reactionModel.networkList[ind] = network
        
        #get the molecular formula for the network
        
        mmol = Molecule()
        for spc in self.source:
            mmol.merge(spc.molecule[0])
            
        form = mmol.getFormula()
        
        #determine T and P combinations
        
        if self.pdepjob.Tlist:
            Tlist = self.pdepjob.Tlist
        else:
            Tlist = np.linspace(self.pdepjob.Tmin,self.pdepjob.Tmax,self.pdepjob.Tcount)
            
        if self.pdepjob.Plist:
            Plist = self.pdepjob.Plist
        else:
            Plist = np.linspace(self.pdepjob.Pmin,self.pdepjob.Pmax,self.pdepjob.Pcount)
            
        #generate the network
        
        explore_tol = self.explore_tol
        
        incomplete = True
        
        while incomplete:
            incomplete = False
            for T in Tlist:
                for P in Plist:
                    if network.getLeakCoefficient(T=T,P=P) > explore_tol:
                        spc = network.getMaximumLeakSpecies(T=T,P=P)
                        flags = np.array([s.molecule[0].getFormula()==form for s in reactionModel.core.species])
                        reactionModel.enlarge((network,spc),reactEdge=False,unimolecularReact=flags,
                                          bimolecularReact=np.zeros((len(reactionModel.core.species),len(reactionModel.core.species))))
                    
                        flags = np.array([s.molecule[0].getFormula()==form for s in reactionModel.core.species])
                        reactionModel.enlarge(reactEdge=True,unimolecularReact=flags,
                                          bimolecularReact=np.zeros((len(reactionModel.core.species),len(reactionModel.core.species))))
                        incomplete = True
        
        
        #clean up output files
        
        path = os.path.join(reactionModel.pressureDependence.outputFile,'pdep')
        for name in os.listdir(path):
            if name.endswith('.py') and '_' in name:
                if int(name.split('_')[-1].split('.')[0]) != len(network.isomers):
                    os.remove(os.path.join(path,name))
                else:
                    os.rename(os.path.join(path,name),os.path.join(path,'network_full.py'))
                    
        
        
        
        
        
        
        
        