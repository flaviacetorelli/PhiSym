#!/bin/bash

### MAIN ###
TEMP=`getopt -o d:g:j:m:l:p:h --long dataset:,globaltag:,iovmap:,lumijson:,mod:,localreco:,proxy:,storage:,template:,merger,dryrun,help -n 'submit.sh' -- "$@"`

if [ $? != 0 ] ; then echo "Options are wrong..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

# Defaults
merger=0
dryrun=0
lumijson=''
mod=''
reco='weights'
proxy='/tmp/x509up_u81383'
storage='T3_IT_MIB'
template='run_RECO_template.py'
user_tmpl=0
help=0

while true; do
    case "$1" in
        -d | --dataset) dataset="$2"; shift 2;;
        -g | --globaltag) globaltag="$2"; shift 2;;        
        -j | --lumijson) lumijson="$2"; shift 2;;
        -m | --mod) mod="$2"; shift 2;;
        -l | --localreco) reco="$2"; shift 2;;        
        -p | --proxy) proxy="$2"; shift 2;;
        -h | --help ) help=1; shift;;
        --iovmap ) cp "$2" IOVmap.root; shift 2;;
        --storage ) storage="$2"; shift 2;;
        --template ) template="$2"; user_tmpl=1; shift 2;;
        --merger ) merger=1; shift;;
        --dryrun ) dryrun=1; shift;;
        -- ) shift; break ;;
        * ) break ;;
    esac
done

if [ $help == 1 ]; then
    echo "Options: "
    echo "-d | --dataset"
    echo "-g | --globaltag"
    echo "-j | --lumijson"
    echo "-m | --mod"
    echo "-l | --localreco"
    echo "-p | --proxy"
    echo "-h | --help"
    echo "--iovmap"
    echo "--storage"
    echo "--template"
    echo "--merger"
    echo "--dryrun"
    exit 0
fi

echo "Submission log:"
echo "---------------"
echo "Dataset:       ${dataset}"
echo "GlobalTag:     ${globaltag}"
echo "Lumi JSON:     ${lumijson}"
echo "CRAB template: ${template}"
echo "Local reco:    ${reco}"
echo "Stage out:     ${storage}"
echo "---------------"

if [ $dryrun != 1 ]; then
    if [ $merger == 0 ]; then
        set -x
        crab submit --proxy ${proxy} ${template} General.requestName="PHISYM-${CMSSW_VERSION}-${reco}-${globaltag}-${mod}" JobType.psetName="PhiSymProducer_${reco}_cfg.py" Site.storageSite="${storage}" Data.inputDataset="${dataset}" Data.lumiMask="${lumijson}"
        set +x
    else
        if [ $user_tmpl == 0 ]; then
           template='run_MERGER_template.py'
        fi
        set -x
        crab submit --proxy ${proxy} ${template} General.requestName="PHISYM-${CMSSW_VERSION}-${reco}-${globaltag}-${mod}" Site.storageSite="${storage}" Data.inputDataset="${dataset}" Data.lumiMask="${lumijson}"
        set +x
    fi
fi

# 80X_dataRun2_Prompt_v4
# /AlCaPhiSym/Commissioning2016-v1/RAW
