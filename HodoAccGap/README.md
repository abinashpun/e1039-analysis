

### clone the analysis repository
```
git clone https://github.com/E1039-Collaboration/seaquest-analysis.git
```

### run the Fun4All macro to perform the simulation
```
cd <path-to-hodo-acc-gap>
root -l Fun4HodoAccGap.C
```
Several root files are generated.
The 'trk_eval.root' contains the information from the evaluator [`TrkEval`](https://e1039-collaboration.github.io/seaquest-offline-doc/de/d2f/classTrkEval.html).
You may open that file and check the contents.

```
$root -l trk_eval.root 
root [0] 
Attaching file trk_eval.root as _file0...
root [1] .ls
TFile**		trk_eval.root	
 TFile*		trk_eval.root	
  KEY: TTree	T;1	TrkEval
```

### analysis

The 'ana.C' is an example macro for analysis.
You may test it using the example data.

```
ln -s /exp/seaquest/app/users/yuhw/seaquest-analysis/HodoAccGap/300cm/*.root .

root -l
root [0] .L ana.C 
root [1] drawRelAcc()
root [2] drawAccPhi()
```

