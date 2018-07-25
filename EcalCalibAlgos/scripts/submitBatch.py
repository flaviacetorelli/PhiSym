#!/usr/bin/python

import sys
import os
import commands
from commands import getstatusoutput
from commands import getoutput
import datetime
import argparse

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def lxbatchSubmitJob (eosdir, iovmap, iov, basedir, outdir, queue, dryrun):
    jobname = outdir+'/PhiSymCalibration_'+queue+'_'+str(iov)+'.sh'
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('export X509_USER_PROXY=/afs/cern.ch/user/s/spigazzi/x509up_u68758 \n\n')
    f.write ('cp '+iovmap+' IOVMap.root \n')
    f.write ('cp '+basedir+'/src/PhiSym/EcalCalibAlgos/test/PhiSymCalibration_cfg.py . \n')
    f.write ('cd '+basedir+' \n')
    f.write ('eval `scramv1 runtime -sh` \n')
    f.write ('cd - \n\n')
    f.write ('PhiSymCalibration PhiSymCalibration_cfg.py eosdirs='+eosdir+' iovmaps=IOVMap.root niovs=1 firstiov='+str(iov)+'\n\n')
    f.write ('rm IOVMap.root \n')
    f.write ('cp *.root '+outdir+'\n')
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    if not dryrun:
        getstatusoutput ('cd '+outdir+'; bsub -q ' + queue + ' ' + '-u $USER@cern.ch ' + jobname + '; cd -')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'Submit PhiSymCalibration jobs: 1 job = 1 IOV')
    parser.add_argument('-d', '--eosdir' , default='', help='dataset to be processed (EOS path to ntuples)')
    parser.add_argument('-m', '--iovmap' , default = '', help='ROOT file with IOV definition')
    parser.add_argument('-q', '--queue' , default = '8nm', help='batch queue (1nh)')
    parser.add_argument('-o', '--outdir' , default = '$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/', help='storage path')
    parser.add_argument('--dryrun' , action="store_true", default=False, help='do not submit the jobs, just create them')
    
    args = parser.parse_args ()
    ## expand output
    stage_out_dir = os.path.abspath(os.path.expandvars(args.outdir))
    print("Output will be copied to: "+stage_out_dir)
    getstatusoutput('mkdir -p '+stage_out_dir)
    ## get list of iovs.
    iovfile = ROOT.TFile.Open(args.iovmap, "READ")
    iovtree = iovfile.Get("iov_map")
    n_jobs = iovtree.GetEntries()

    ## create jobs
    basedir = os.path.expandvars("$CMSSW_BASE")
    print("Submitting "+str(n_jobs)+" jobs to queue "+args.queue+"...")
    for iov in range(0, n_jobs):
        lxbatchSubmitJob(args.eosdir, os.path.abspath(args.iovmap), iov, basedir, stage_out_dir, args.queue, args.dryrun) 
