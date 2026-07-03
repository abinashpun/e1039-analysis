## Fun4AllVectEventInputManager
The `Fun4AllVectEventInputManager` is an input file manager in the Fun4All framework, designed to manage and register various types of input data from a simple vector-based ROOT file, such as triggers, hits, and event information, into the framework's event processing system.

## Fun4AllVectEventOutputManager
The `Fun4AllVectEventOutputManager` is an output file manager in the Fun4All framework, designed to store triggers, hits, and event information in a vector-based file format.

## Goal
The Fun4All framework does not support random indexing in the Fun4All event processing server. Therefore, translating the data into a convenient format enables us to perform pre- and post-processing using the input/output managers.

# Updated Variable Names and Types

## Event-Level Variables
| Variable Name      | Type               | Description                          |
|--------------------|--------------------|--------------------------------------|
| `eventID`          | `int`              | Unique identifier for the event      |
| `runID`            | `int`              | Identifier for the current run       |
| `spillID`          | `int`              | Identifier for the spill in the run  |
| `fpgaTriggers`     | `int[5]`           | Array of FPGA trigger counts         |
| `nimTriggers`      | `int[5]`           | Array of NIM trigger counts          |
| `intensity`        | `int[33]`          | Array for QIE RF intensities         |

## Hit-Level Variables
| Variable Name          | Type                     | Description                                  |
|------------------------|--------------------------|----------------------------------------------|
| `detectorID`           | `std::vector<int>`       | Detector IDs for all hits                    |
| `elementID`            | `std::vector<int>`       | Element IDs associated with each hit         |
| `driftDistance`        | `std::vector<double>`    | Drift distances for each hit                 |
| `tdcTime`              | `std::vector<double>`    | TDC timing values for each hit               |
| `hitInTime`            | `std::vector<bool>`      | Flags indicating if hits are within time     |

## Trigger-Level Hit Variables
| Variable Name              | Type                     | Description                                  |
|----------------------------|--------------------------|----------------------------------------------|
| `trigDetectorID`           | `std::vector<int>`       | Detector IDs for triggered hits              |
| `trigElementID`            | `std::vector<int>`       | Element IDs for triggered hits               |
| `trigDriftDistance`        | `std::vector<double>`    | Drift distances for triggered hits           |
| `trigTdcTime`              | `std::vector<double>`    | TDC timing values for triggered hits         |
| `trigHitInTime`            | `std::vector<bool>`      | Flags indicating if triggered hits are in time |


``` Compilation
source Fun4AllVectFileManager ../setup.sh
cmake-this
make-this
```

### Uses

1. **Convert the DST file to a vector file**  
   If you'd like to convert the DST file to a vector format for easier event processing, use the following commands:

   ```bash
   cd Convert
   root -b ConvertToVect.C
   ```
2. **Run the Pre- or Post-Process Data in the Fun4All Server**

   After completing the event processing, if you would like to add the triggers, hits, and event information into the Fun4All framework's event processing server, follow these steps:

   Go to the `work_reco` directory, where you will find two subdirectories:

   - **(i) Reco_DST_in**:  
     Run `Fun4Sim.C` and make sure to use the correct input file in DST format. (This directory was created to verify that both input file systems yield the same results in the output.)

   - **(ii) Vect_DST_in**:  
     Run `Fun4Sim.C` and make sure to use the correct input file in vector format, which you obtained from the Convert directory. You can definitely perform preprocessing as required for your analysis before registering the data into the Fun4All framework and running reconstruction and vertexing.

