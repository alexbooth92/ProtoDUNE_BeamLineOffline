#!/bin/sh

g++ -std=c++11 -o Module_MakeEventLevelTree.exe Module_MakeEventLevelTree.C BeamLineClasses/class_AcquisitionXBPF.C BeamLineClasses/class_AcquisitionXBTF.C BeamLineClasses/class_Detector.C BeamLineClasses/class_BeamLine.C BeamLineClasses/class_DetectorCoincidences.C BeamLineClasses/class_Spill.C BeamLineClasses/unpack.C BeamLineClasses/plotting.C BeamLineClasses/Analyse.C `root-config --cflags --glibs`
