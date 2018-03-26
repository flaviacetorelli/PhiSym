#!/usr/bin/python

import os
import re
import json
import subprocess
import collections

from optparse import OptionParser

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('', '--debug', dest='debug', action='store_true')
    parser.add_option('-d', dest='dataset', type='string', default='', help='merge LS info for this dataset')
    parser.add_option('-p', '--pu-file', dest='pu_file', type='string', default='',
                      help='Pre-generated PU info from brilcalc file')    
    parser.add_option('-o', '--output', dest='output', type='string', default='',
                      help='specify a custom output file, by default the merged json will be put dataset directory if it is on EOS')
    (options, args) = parser.parse_args()

    ###---Get list of LS info files
    files = []
    print "Getting files from DAS for dataset "+options.dataset
    cmd_DAS_query = subprocess.Popen(
        ["dasgoclient --query='file dataset="+options.dataset+" instance=prod/phys03' --limit 0 | grep '/store/'"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    file_list_string, err_file_list = cmd_DAS_query.communicate()
    files.extend([re.sub('phisym.*lumis', 'phisym_lumi_info_json', '/eos/cms'+ifile) for ifile in file_list_string.split("\n")])
    files.pop()
    if options.debug:
        for ifile in files:
            print(ifile)
    eosdir = os.path.abspath(os.path.dirname(files[0])+"/../")
            
    ###---Read PU data file or get data from brilcalc if file is not specified
    pu_data = {}
    pu_file_str = options.pu_file
    if pu_file_str == "":
        cmd_runs = subprocess.Popen(
            ["dasgoclient --query='run dataset="+options.dataset+" instance=prod/phys03'  --limit 0"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        runs_string, err_runs = cmd_runs.communicate()
        runs_string = runs_string.replace("\n", ",")[:-1]
        if options.debug:
            print("List of runs read from: "+options.dataset, runs_string)
        cmd_pu = subprocess.Popen(
            ['${CMSSW_BASE}/src/PhiSym/EcalCalibAlgos/scripts/get_pu_info.sh '+runs_string+' > /tmp/$USER/pu_list.json'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        pu_string, err_pu = cmd_pu.communicate()
        pu_file_str = '/tmp/$USER/pu_list.json'
        if options.debug:
		print 'Pile up file is ', pu_file_str
    with open(os.path.expandvars(pu_file_str)) as pu_file:
        pu_data = json.load(pu_file)
            
    ###---Read LS info files and combine LS and PU info
    data = {}
    for i, ifile in enumerate(files):
        with open(ifile) as json_file:
            print "Reading files... (", i, "/", len(files), ")\r",
            new_data = json.load(json_file)
            for time, info in new_data.items():
                run_num = str(info[0])
                lumi_num = str(info[1])
                if run_num in pu_data.keys() and lumi_num in pu_data[run_num].keys():
                    if float(pu_data[run_num][lumi_num]) > 0:
                        info[3] = info[3]*float(pu_data[run_num][lumi_num])
                        data[time] = info

    merged_json = eosdir+'/ls_info.json' if options.output == "" else options.output
    with open(merged_json, "w") as json_file:
        json.dump(data, json_file, sort_keys=True, indent=4, separators=(',', ' : '))

    print "#########Summary##################", ""
    print "\033[1;34m LS summary \033[0;10m json file stored in: \033[1;34m"+merged_json+"\033[0;10m"
    print "\033[1;34m PU summary \033[0;10m file: \033[1;34m"+os.path.expandvars(pu_file_str)+"\033[0;10m"
    if options.pu_file == "":
        print "\033[1;34m NOTE \033[0;10m: You can copy the PU file to speed up further reprocessing passing it to this script with the -p option"
