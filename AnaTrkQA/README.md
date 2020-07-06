## AnaTrkQA
Work on progress anslysis moudle for quality assurance of tracks.

## Building

The analysis code is located in `e1039-analysis/AnaTrkQA/src`.
The following commands set up the shell environment and compile the code;
```
cd e1039-analysis/AnaTrkQA
source setup.sh
cd build
cmake-ana-TrkQA
make install
```
### Running the analysis macro

Load libana_trkqa.so in Fun4Sim macro in SimChainDev area and register the AnaTrkQA after the recontruction process as
```
AnaTrkQA* trackQA = new AnaTrkQA();
trackQA->set_hit_container_choice("Vector");
trackQA->set_out_name("trackQA.root");
se->registerSubsystem(trackQA);
```
### Plotting macors
 Plotting macro `plot_trackQA.C` is situated inside macros directory. It loads `trackQA.root` produced by Fun4Sim macro and plots various histograms and graphs for quality assurance of the tracks.
