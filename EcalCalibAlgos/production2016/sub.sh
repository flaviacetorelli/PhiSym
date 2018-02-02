#!/bin/bash

# for era in B C D E;
# do

#     ./submit.sh -d /AlCaPhiSym/Run2017${era}-v1/RAW -g 92X_dataRun2_Prompt_v9 -m Run2017${era}_NovReReco_v2 -j Cert_294927-305185_13TeV_PromptReco_Collisions17_JSON.txt -l multifit --storage T2_CH_CERN

# done

for era in C;
do
    ./submit.sh -d /AlCaPhiSym/dpg_ecal-crab_PHISYM-CMSSW_9_2_12-multifit-92X_dataRun2_Prompt_v9-Run2017${era}_NovReReco_v2-5020d2c0332ae0a06a741c48c9f15d93/USER -g 92X_dataRun2_Prompt_v9 --template run_MERGER_template.py --merger --storage T2_CH_CERN -m Run2017${era}_NovReReco_v2-merged_v1 --iovmap=../ntuples/Run2017${era}_NovReReco_v1/IOVMap2017${era}.root -l multifit

done

