# TargetSim

The instructions below are for **seaquestgpvm01.fnal.gov** where the **seaquest-offline** libraries are compiled at **/exp/seaquest/app/users/yuhw/install**.

## setup environment
```
source /exp/seaquest/app/users/yuhw/setup.sh
```

## check out analysis repository
```
git clone https://github.com/E1039-Collaboration/seaquest-analysis.git
```

## analysis module
'TruthEval' contains the analysis module.


## How to run
```
cd <path-to-seaquest-analysis/TargetSim>
root -l Fun4All_G4_E1039_R2.C
```

- 'Fun4All_G4_E1039_R2.C' is the main macro to run the simulation
- 'ana.C' is the ploting macro
