#!/usr/env/bin python
#
# ReFlex3D: Refined Flexible Alignment of Molecules using Shape and Electrostatics
# by Thomas Christian Schmidt, David A Cosgrove and Jonas Bostrom
# https://pubs.acs.org/doi/10.1021/acs.jcim.7b00618
# 30/03/2018

from __future__ import print_function
from openeye import *
from math import *

class MultiOptMol(object):
    def __init__(self,mol):
        self.mol=mol.CreateCopy()
        print("MM mol initialized with ",self.mol.NumAtoms()," atoms")
        #self.epsilon    =   80.
        #self.saltconc   =   .08
        self.alpha=1.
        self.beta =1.

        ### initialize Szybki
        self.szybkiopts=OESzybkiOptions()
        self.szybkiopts.SetRunType(OERunType_CartesiansOpt)
        self.szybki=OESzybki(self.szybkiopts)
        self.szybkiresults=OESzybkiResults()
        ### compute reference energy
        tmol=self.mol.CreateCopy()
        self.szybki(tmol,self.szybkiresults)
        self.refenergy=self.szybkiresults.GetTotalEnergy()
        print("reference energy set to ",self.refenergy)

        ### configure Shape Overlap Comparison
        self.overlap=OEOverlap()                                                # mandatory
        self.overlapresults=OEOverlapResults()                                  # mandatory
        self.overlap.SetMethod(OEOverlapMethod_Exact)                           # optional
        self.overlap.SetUseHydrogens(True)                                      # optional
        self.coloroverlapresults=OEColorResults()                               # mandatory
        self.coloroverlap=OEColorOverlap()                                      # mandatory
        self.coloroverlap.SetColorForceField(OEColorFFType_ImplicitMillsDean)   # optional

    def GetCoords(self):
        cvec=[]
        for atom in self.mol.GetAtoms():
            coords=self.mol.GetCoords(atom)
            for i in range(3):
                cvec.append(coords[i])
        return cvec

    def SetRefMol(self,mol):
        self.refmol=mol.CreateCopy()
        print("reference molecule set")
        self.overlap.SetRefMol(self.refmol)                                     # mandatory
        self.coloroverlap.SetRefMol(self.refmol)                                # mandatory

    def MMscore(self,cvec):
        try:
            self.refmol
        except:
            exit('no similarity without reference molecule')

        ### create local copy of mol and assign coordinates
        for i,atom in enumerate(self.mol.GetAtoms()):
            coords=(cvec[3*i],cvec[3*i+1],cvec[3*i+2])
            self.mol.SetCoords(atom,coords)


        ### compute relative energy
        self.szybki.SetRunType(OERunType_SinglePoint)
        self.szybki(self.mol,self.szybkiresults)
        self.energy=self.szybkiresults.GetTotalEnergy()-self.refenergy
#       print("relative energy is ",self.energy)

        ### compute shape similarity
        self.overlap.Overlap(self.mol,self.overlapresults)
        shapescore= -self.overlapresults.GetTversky(self.alpha,self.beta)
        #print("shape overlap: ",shapescore)

        ### compute color similarity
        self.coloroverlap.ColorScore(self.mol,self.coloroverlapresults)
        colorscore= -self.coloroverlapresults.GetTversky(self.alpha,self.beta)
        #print("color overlap: ",colorscore)

        score=0.
        score += self.energy
        score += shapescore * 250 * exp(-.01 * self.energy)
        score += colorscore * 250 * exp(-.01 * self.energy)
        #print("score is ",score)

        return score

