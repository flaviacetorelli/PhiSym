#!/usr/bin/python

import sys
import os
import commands
from commands import getstatusoutput
from commands import getoutput
import datetime
import argparse
import subprocess

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def lxbatchSubmitJob(eosdir, iovmap, iov, basedir, outdir, queue, proxy, dryrun):
    jobname = outdir+'/PhiSymCalibration_'+queue+'_'+str(iov)+'.sh'
    logname = jobname.replace('.sh', '.log')
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('export X509_USER_PROXY='+proxy+' \n\n')
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
        getstatusoutput('cd '+outdir+'; bsub -q ' + queue + ' -oo ' + logname + ' -u "" ' + jobname + '; cd -')

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def resubmitJobs(outdir, queue):
    """
    Check log files in directory outdir and resubmit failed jobs
    """

    getFailedCmd = subprocess.Popen(
        ['grep -r "Exited with exit code" '+outdir+' | sed -e s/:.*$//g'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    logs_list, err = getFailedCmd.communicate()
    logs_list = logs_list.split('\n')
    logs_list.pop()
    if len(logs_list) == 0:
        print("No job to resubmit, calibration complete!")
    for log in logs_list:
        jobname = log.replace('.log', '.sh')
        print log
        getstatusoutput ('cd '+outdir+'; bsub -q ' + queue + ' -oo ' + log + ' -u "" ' + jobname + '; cd -')
        
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'Submit PhiSymCalibration jobs: 1 job = 1 IOV')
    parser.add_argument('-d', '--eosdir' , default='', help='dataset to be processed (EOS path to ntuples)')
    parser.add_argument('-m', '--iovmap' , default = '', help='ROOT file with IOV definition')
    parser.add_argument('-q', '--queue' , default = '8nm', help='batch queue (1nh)')
    parser.add_argument('-o', '--outdir' , default = '$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/', help='storage path')
    parser.add_argument('--resub' , action="store_true", default=False, help='resubmit failed jobs, requires -o')
    parser.add_argument('--dryrun' , action="store_true", default=False, help='do not submit the jobs, just create them')
    
    
    args = parser.parse_args ()

    ## expand output
    stage_out_dir = os.path.abspath(os.path.expandvars(args.outdir))
    print("Output will be copied to: "+stage_out_dir)
    getstatusoutput('mkdir -p '+stage_out_dir)

    ## get GRID proxy
    # check validity
    stat,out = commands.getstatusoutput("voms-proxy-info -e --valid 2:00")
    if stat:
        raise Exception,"voms proxy not found or validity less than 2 hours:\n%s.\n\033[1;34m Please run: voms-proxy-init -voms cms\033[0;10m" % out
    # get proxy name
    stat,proxy = commands.getstatusoutput("voms-proxy-info -p")
    if stat:
        raise Exception,"Unable to voms proxy:\n%s" % proxy
    commands.getstatusoutput("cp %s %s/" % (proxy, stage_out_dir))
    proxy = stage_out_dir+"/"+proxy.split("/")[2]

    ## resubmit failed jobs
    if args.resub:
        resubmitJobs(stage_out_dir, args.queue)
    else:
        ## get list of iovs.
        iovfile = ROOT.TFile.Open(args.iovmap, "READ")
        iovtree = iovfile.Get("iov_map")
        n_jobs = iovtree.GetEntries()
        
        ## create jobs
        basedir = os.path.expandvars("$CMSSW_BASE")
        print("Submitting "+str(n_jobs)+" jobs to queue "+args.queue+"...")
        for iov in range(0, n_jobs):
            lxbatchSubmitJob(args.eosdir, os.path.abspath(args.iovmap), iov, basedir, stage_out_dir, args.queue, proxy, args.dryrun) 
