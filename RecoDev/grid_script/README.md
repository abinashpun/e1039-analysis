```
gridrun_data.sh
gridsub_data.sh
```
 These are scripts to run real data in grid with the macro, `RecoE1039Data.C`. Since `/data2` area is not acessible during job submission, this script assumes you copied the dst data files from `/data2` area to `/pnfs/e906/persistent/users/$USER/DstRun` area. Please refer to the [wiki page](https://github.com/E1039-Collaboration/e1039-wiki/wiki/data-file-on-grid) about setting the file location. To run list of the runs from `run_list.txt` in grid, use the following command;
```
./gridsub_data.sh ouput_dir 1 run_list.txt 0 single > test_single-grid.log
```
Here, the "0" options lets you run all the available events in the run. For the longer runs which can't be finished within the 24 hours limit, you need to split the data in some spill number interval. The `splitDST.C`macro
splits the longer DST runs in specified number of spills. For now, the default number of spills is set to 10. 

```
root -b -q splitDst.C\(run_number\)
```

The splitted dst files are outputed to `/pnfs/e906/persistent/users/$USER/DstRun/run_number/` (You might see some error msg regarding the persmission to write in persistent area, please ignore the message for now). After splitting the data for particular run, use the following command to run all the splitted files from the runs in run_list.txt;

```
./gridsub_data.sh ouput_dir 1 run_list.txt 0 split > test_split-grid.log
```

