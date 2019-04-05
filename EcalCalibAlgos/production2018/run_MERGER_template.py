from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName     = 'PHISYM'
config.General.transferLogs    = True
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName      = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName        = 'PhiSymMerger_cfg.py'
config.JobType.inputFiles      = ['IOVmap.root']
config.JobType.priority        = 20

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
config.Data.inputDataset   = ''

#config.Data.useParent = True
config.Data.inputDBS           = 'phys03'
config.Data.splitting          = 'LumiBased'
config.Data.lumiMask           = ''
config.Data.unitsPerJob        = 100
config.Data.totalUnits         = -1
config.Data.publication        = False
#config.Data.ignoreLocality     = True

# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase      =  '/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite        = 'T2_CH_CERN'
#config.Site.blacklist          = ['T2_US_Nebraska']
config.Site.whitelist          = ['T1_IT_CNAF', 'T2_CH_CERN', 'T2_IT_Bari']
