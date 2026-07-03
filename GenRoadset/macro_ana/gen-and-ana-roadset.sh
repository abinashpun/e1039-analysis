#!/bin/bash

root -b GenerateRoad.C
root -b ApplyRoadset2Signal.C
root -b 'ApplyRoadset2BG.C(0)'
root -b 'ApplyRoadset2BG.C(45000)'
root -b 'ApplyRoadset2BG.C(63000)'
root -b 'ApplyRoadset2BG.C(80500)'
root -b DrawRS.C
