#!/bin/bash
## Script to make lists of input signal files and input BG files.

function MakeOneList {
    local -r DIR_IN=$1
    local -r  FN_IN=$2
    local -r FN_OUT=$3
    echo "DIR_IN = $DIR_IN"
    echo "FN_IN  = $FN_IN"
    echo "FN_OUT = $FN_OUT"
    find $DIR_IN -name $FN_IN | sort >$FN_OUT
    echo "N of inputs found = $(cat $FN_OUT | wc -l)"
    echo
}

#echo "## Signal List: Reverse KMag polarity"
#MakeOneList ../macro_gen_signal/data/main_reverse signal_tree.root list_signal_reverseKMAG.txt

#echo "## Signal List: Normal KMag polarity"
#MakeOneList ../macro_gen_signal/data/main_v2 signal_tree.root list_signal.txt

#echo "## Signal List: Reverse KMag polarity, H1X gap"
#MakeOneList ../macro_gen_signal/data/main_reverse_H1Xgap signal_tree.root list_signal_reverseKMAG.txt

#echo "## Signal List: Normal KMag polarity, H1X gap"
#MakeOneList ../macro_gen_signal/data/main_v2_H1Xgap signal_tree.root list_signal.txt

#echo "## Signal List: No KMag field, H1X gap"
#MakeOneList ../macro_gen_signal/data/main_zero_H1Xgap signal_tree.root list_signal_zeroKMAG.txt

echo "## Signal List: KMag = 50% Norm, H1X gap"
MakeOneList ../macro_gen_signal/data/main_normal050 signal_tree.root list_signal_KMagNormal050.txt

#echo "## BG List: Reverse KMag polarity"
#MakeOneList ../macro_full_bg/data/main_01_nim3pot_reverseMag bg_data.root list_bg_fullsimRun06_reverseKMAG.txt

#echo "## BG List: Normal KMag polarity"
#MakeOneList ../macro_full_bg/data/main_01_nim3potv2 bg_data.root list_bg_fullsimRun06.txt

#echo "## BG List: Reverse KMag polarity, H1X gap"
#MakeOneList ../macro_full_bg/data/main_01_nim3pot_reverseMag_H1Xgap bg_data.root list_bg_fullsimRun06_reverseKMAG.txt

#echo "## BG List: Normal KMag polarity, H1X gap"
#MakeOneList ../macro_full_bg/data/main_01_nim3potv2_H1Xgap bg_data.root list_bg_fullsimRun06.txt

#echo "## BG List: Normal KMag polarity, H1X gap"
#MakeOneList ../macro_full_bg/data/main_01_nim3pot_normal050 bg_data.root list_bg_KMagNormal050.txt

