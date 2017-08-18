import matplotlib as mpl
import nuSQUIDSpy as nsq
import matplotlib.pyplot as plt
import sys

import nuSQUIDSTools
import numpy as np
mpl.rc('font', family='serif', size=20)

units = nsq.Const()
interactions = False

E_min = 1.*units.GeV
E_max = 1.e3*units.GeV
E_nodes = 100 #nodes for constructor
Ebins = np.logspace(np.log10(E_min),np.log10(E_max),E_nodes)
Enum = 1000 #bins for print out

czmin=-1.
czmax=1.
Ncz = 120 #nodes for constructor
czbins = np.linspace(czmin,czmax,Ncz)
cnum = 200 #bins for print out
neutrino_flavors = 3

#NSI params
epsilon_mutau = 0.010
#epsilone_mumu = 0.
#epsilon_tautau = 0.0
#epsilon_emu = 0.
#epsilon_ee = 0.0
#epsilon_etau = 0.

nus_atm = nsq.nuSQUIDSNSIAtm(czbins,Ebins,neutrino_flavors,nsq.NeutrinoType.neutrino,interactions)

nus_atm.Set_MixingAngle(0,1,0.585732); #33.56 degrees (from nu-fit.org)
nus_atm.Set_MixingAngle(0,2,0.147655); #8.46 degrees (from nu-fit.org)
nus_atm.Set_MixingAngle(1,2,0.726057); #41.6 degrees (from nu-fit.org)
nus_atm.Set_SquareMassDifference(1,7.5e-05); #7.50 x 10^-5 (from nu-fit.org)
nus_atm.Set_SquareMassDifference(2,0.002524+7.5e-05); #2.524 x 10^-3 (from nu-fit.org)

for icz in range(Ncz):
    nus_atm.GetnuSQuIDS(icz).Set_epsilon_mutau(epsilon_mutau)
    #nus_atm.GetnuSQuIDS(icz).Set_epsilon_etau(epsilon_tautau)
    #nus_atm.GetnuSQuIDS(icz).Set_epsilon_memu(epsilon_mumu)
    #nus_atm.GetnuSQuIDS(icz).Set_epsilon_etau(epsilon_etau)
    #nus_atm.GetnuSQuIDS(icz).Set_epsilon_emu(epsilon_emu)
    #nus_atm.GetnuSQuIDS(icz).Set_epsilon_ee(epsilon_eee)

nus_atm.Set_rel_error(1.0e-6)
nus_atm.Set_abs_error(1.0e-6)

#Set initial state
inistate = np.zeros((nus_atm.GetNumCos(),nus_atm.GetNumE(),neutrino_flavors))
for ci in range(0,nus_atm.GetNumCos()):
        for ei in range(0,nus_atm.GetNumE()):
            for flv in range(0,neutrino_flavors):
                if flv ==1:
                    inistate[ci][ei][flv] = 1
                else:
                    inistate[ci][ei][flv] = 0

nus_atm.Set_initial_state(inistate,nsq.Basis.flavor)

nus_atm.EvolveState()

if (len(sys.argv) > 1):
    filename = ("%s.txt"%sys.argv[1])
else:
    filename = "fluxes_flavor_python_nsi.txt"

cz_range = np.linspace(czmin,czmax,cnum)
E_range = np.logspace(np.log10(E_min),np.log10(E_max),Enum)

f = open(filename,"wb")
f.write("# log10(E) cos(zenith) E flux_e flux_mu flux_tau\n")
for cz in cz_range:
    for ei in E_range:
        energy = np.log10(ei/units.GeV)
        nu_e = nus_atm.EvalFlavor(0,cz,ei)
        nu_mu = nus_atm.EvalFlavor(1,cz,ei)
        nu_tau = nus_atm.EvalFlavor(2,cz,ei)
        f.write("%f %f %f %f %f %f\n"%(energy,cz,ei/units.GeV,nu_e,nu_mu,nu_tau))
    f.write("\n")
f.close()
