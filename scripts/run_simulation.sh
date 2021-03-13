#!/bin/bash

npsim --runType batch --numberOfEvents 100000 \
      --compactFile ./topside.xml \
      --inputFiles  ../datasets/data/emcal_electron_0GeVto30GeV_100kEvt.hepmc \
      --outputFile  ./sim_electron_0GeVto30GeV_100k_input.root
