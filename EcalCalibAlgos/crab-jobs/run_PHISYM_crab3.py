from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName     = 'PHISYM_test_production_Run2015A_v3'
config.General.transferLogs    = True
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName      = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName        = 'PhiSymProducer_cfg.py'
config.JobType.priority        = 20

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
config.Data.inputDataset       = '/AlCaPhiSym/Run2015A-v1/RAW'
#config.Data.useParent = True
config.Data.inputDBS           = 'global'
config.Data.splitting          = 'LumiBased'
config.Data.lumiMask           = 'Run2015A_v0_filtered.json'
config.Data.unitsPerJob        = 750
config.Data.totalUnits         = -1
config.Data.publication        = True
#config.Data.ignoreLocality     = True

# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase      =  '/store/user/spigazzi/'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite        = 'T3_IT_MIB'

