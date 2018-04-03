#!/usr/env/bin python
#
# ReFlex3D: Refined Flexible Alignment of Molecules using Shape and Electrostatics
# by Thomas Christian Schmidt, David A Cosgrove and Jonas Bostrom
# https://pubs.acs.org/doi/10.1021/acs.jcim.7b00618
# 30/03/2018
#

from __future__ import print_function
from numpy import *
from scipy import optimize
import sys
from openeye import *
from scipy.optimize import minimize
from scimodule import *

# # # # # # # # # # # # # # # ################ # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # ###   MAIN   ### # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # ################ # # # # # # # # # # # # # # # #
def main(argv=[__name__]):

    if len(sys.argv) < 4:
        print("Usage\t 'python ./reflex3d.py refMol.oeb dbMol.oeb aligned.oeb'")
        print("Example\t 'python ./reflex3d.py 1fkgLigand.oeb 1fkhLigand.oeb 1fkhLigandAligned.oeb'")
        sys.exit(1)

    OEThrow.SetLevel(99)

    ifs=oemolistream(sys.argv[1])
    refmol=OEMol()
    OEReadMolecule(ifs,refmol)

    ifs=oemolistream(sys.argv[2])
    mol=OEMol()
    molecules=[]
    while OEReadMolecule(ifs,mol):
        molecules.append(mol.CreateCopy())

    ofs=oemolostream(sys.argv[3])
    print('all done reading input')



    for mol in molecules:
        print(mol.GetTitle())

        for conf in mol.GetConfs():
            print(conf.GetIdx()+1,' of ',mol.NumConfs())

            ### initialize new object wrapped around OEMol()
            OEAddExplicitHydrogens(conf)
            mopt=MultiOptMol(OEGraphMol(conf))
            cvec=mopt.GetCoords()
            mopt.SetRefMol(refmol)

            ### loop over iterations - check for convergence
            res=minimize(mopt.MMscore,cvec,method='L-BFGS-B',jac=None,\
                        tol=1.e-3,\
                        callback=None,\
                        options={\
                                 'disp': True,\
                                 'iprint':  -1,\
                                 'gtol':    1e-04,\
                                 'eps':     1e-03,\
                                 'maxiter': 15000,\
                                 'ftol':    2.220446049250313e-09,\
                                 'maxcor':  10,\
                                 'maxfun':  15000})


            print("Optimization converged?:  ",res.success)
            #print(res)


            print("success: ",res.success)
            print("status: ",res.status)
            print("message: ",res.message)
            print("fun: ",res.fun)
            print("nfev: ",res.nfev)
            print("nit: ",res.nit)

            ### write optimized coordinates back to mol
            for i,atom in enumerate(conf.GetAtoms()):
                    coords=(res.x[3*i],res.x[3*i+1],res.x[3*i+2])
                    conf.SetCoords(atom,coords)
        print('###')
        print('molecule ',mol.GetTitle(),' done')
        print('###')
        OEWriteMolecule(ofs,mol)
        ofs.flush()

    #os.unlink('iterate.dat')



    print('\n all done\n')
# # # # # # # # # # # # # ####################### # # # # # # # # # # # # # #
# # # # # # # # # # # # # ###   END OF MAIN   ### # # # # # # # # # # # # # #
# # # # # # # # # # # # # ####################### # # # # # # # # # # # # # #
if __name__ == "__main__":
    sys.exit(main(sys.argv))
