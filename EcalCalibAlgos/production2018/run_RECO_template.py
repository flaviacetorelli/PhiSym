from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName       = 'PHISYM'
config.General.transferLogs      = True
config.General.transferOutputs   = True

config.section_('JobType')
config.JobType.pluginName        = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName          = 'PhiSymProducer_cfg.py'
#config.JobType.inputFiles        = ['alphas_eflow2012.db']
config.JobType.outputFiles        = ['phisym_lumi_info_json.root']
config.JobType.priority          = 30

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
config.Data.inputDataset         = ''
config.Data.inputDBS             = 'global'
config.Data.splitting            = 'LumiBased'
config.Data.lumiMask             = ''
config.Data.unitsPerJob          = 15
config.Data.totalUnits           = -1
config.Data.publication          = True
config.Data.ignoreLocality       = True
#config.Data.isbchecksum          = 'aa6d0f694fbde2c60e6338a92e82f36e'

# This string is used to construct the output dataset name
config.Data.publishWithGroupName = True
config.Data.outLFNDirBase        = '/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite          = 'T2_CH_CERN'
#config.Site.blacklist            = ['']
config.Site.whitelist            = ['T1_IT_CNAF', 'T2_CH_CERN', 'T1_US_FNAL']
