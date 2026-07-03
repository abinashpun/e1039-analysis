## Usage

### Setup

On seaquestgpvm01
```
source /exp/seaquest/appa/users/yuhw/setup.sh
```

### checkout analysis repository

```
git checkout git@github.com:HaiwangYu/seaquest-analysis.git
```

### Run simulation macro
```
cd <path-to/seaquest-analysis/Pythia8Demo>
root -l -q Fun4Pythia8Demo.C //Pythia only or
root -l -q Fun4Sim.C //Geant4 simulation
```
The output DSTRead.root file contanins a tree with TClonesArray branches.
TClonesArray may constains:
- [PHG4VtxPoint](https://haiwangyu.github.io/seaquest-offline-doc/d6/d81/classPHG4VtxPoint.html) for interaction vertex
- [PHG4Particle](https://haiwangyu.github.io/seaquest-offline-doc/de/dc9/classPHG4Particle.html) for truth particles
- [PHG4Hit](https://haiwangyu.github.io/seaquest-offline-doc/d3/d9e/classPHG4Hit.html) for Geant4 Hit

You could use methods after loading the dinamic library.
The [ana.C](https://github.com/HaiwangYu/seaquest-analysis/blob/master/Pythia8Demo/ana.C) is a simple example looping over PHG4Particle TClonesArray.
Other TClonesArray branches can be accessed similarly.


Refer: https://haiwangyu.github.io/seaquest-offline-doc/index.html

### Run example analysis macro
```
root -l -q ana.C
```


