#!/bin/bash

DATASET=/AlCaPhiSym/Run2016G-v1/RAW
INSTANCE=global
if [ ! -z "${2}" ]; then    
    DATASET=${2}
fi
if [ ! -z "${3}" ]; then    
    INSTANCE=${3}
fi
    
RUNLIST=`../scripts/das_client.py --query "lumi run dataset=${DATASET} instance=prod/${INSTANCE} | sort run.run_number" --limit 0 --secondary-key lumi.number`

echo -ne '{\n '$RUNLIST'\n}' | sed 's:]] :]],\n :g' | sed 's:\[\[:\: \[\[:g' | sed 's: \::" \::g' | sed 's:^[^{}]: ":g' > ${1}

exit 0
