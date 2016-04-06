from WMCore.Configuration import Configuration

#recoType = "weights"
recoType = "multifit"

config = Configuration()

config.section_('General')
config.General.requestName       = 'PHISYM-CMSSW_7415-'+recoType+'-74X_dataRun2_Prompt_v4-Run2015B_v2'
config.General.transferLogs      = True
config.General.transferOutputs   = True

config.section_('JobType')
config.JobType.pluginName        = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName          = 'PhiSymProducer_'+recoType+'_cfg.py'
config.JobType.inputFiles      = ['ecaltemplates_popcon_weekly_best.db', 'ecalcovariances_popcon_weekly_best.db']
config.JobType.priority          = 30

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
config.Data.inputDataset         = '/AlCaPhiSym/Run2015B-v1/RAW'
config.Data.inputDBS             = 'global'
config.Data.splitting            = 'LumiBased'
config.Data.lumiMask             = 'Run2015ABCD_DCSONLY.json'
config.Data.unitsPerJob          = 10
config.Data.totalUnits           = -1
config.Data.publication          = True
#config.Data.isbchecksum        = 'aa6d0f694fbde2c60e6338a92e82f36e'

# This string is used to construct the output dataset name
# config.Data.publishWithGroupName = True
# config.Data.outLFNDirBase        = '/store/group/dpg_ecal/alca_ecalcalib/phiSymmetry/'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite          = 'T3_IT_MIB'

