from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'CMSSW_7_1_25_patch5_GENntuple_ggFHgg_MINLOHJJ_20July2017' 
config.General.workArea = 'crab'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/GenPhotonTuplizer_71X_GEN.py'
config.JobType.outputFiles = ['genPhotonNtuple.root']

config.section_("Data")

config.Data.inputDataset = '/GluGluHToGG_M125_13TeV_powheg2_minloHJJ_pythia8/RunIIWinter15GenOnly-MCRUN2_71_V1-v2/GEN'
#config.Data.lumiMask = 'data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.Data.inputDBS = 'global' #change this according to the DBS instance (usually 'global') of the target dataset
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 500
config.Data.publication = False
config.Data.ignoreLocality = True #enable AAA

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/user/zhicaiz/GenPhotonNtuples/MINLOHJJ_20July2017/'
