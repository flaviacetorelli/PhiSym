import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.applyCorrections = cms.bool(False)
process.variables = cms.vstring("IC")
process.ioFilesOpt = cms.PSet(    
    inputFiles = cms.vstring(

        # 2012 ABCD
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_193834_194050.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_194051_194120.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_194150_194210.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_194223_194428.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_194429_194479.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_194480_194644.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_194691_194897.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_194912_195115.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_195147_195390.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_195396_195552.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_195633_195868.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_195913_195950.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_195963_196239.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_196249_196438.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_196452_196531.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_198202_198249.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_198268_198522.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_198941_199008.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_199011_199336.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_199356_199436.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_199698_199754.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_199804_199877.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_199960_200091.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_200152_200466.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_200473_200601.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_200990_201115.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_201159_201176.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_201278_201625.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_201657_201824.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_202012_202075.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_202084_202209.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_202237_202477.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_202478_202478.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012ABCD_timeEvo/summed_202504_202504.root"

        # 2015CD history weights
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_254790-346_254790-611.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_254790-612_254790-784.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_254790-90_254790-345.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256677-366_256729-326.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256729-1100_256729-1510.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256729-1511_256801-84.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256729-327_256729-708.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256729-709_256729-1099.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256801-85_256843-158.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256843-1300_256868-200.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256843-159_256843-528.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256843-529_256843-1299.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256868-202_256926-35.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_256926-36_256941-294.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257394-41_257397-56.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257397-57_257399-271.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257400-1354_257487-220.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257400-1_257400-357.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257400-358_257400-735.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257400-736_257400-1353.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257487-1069_257490-254.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257487-221_257487-492.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257487-493_257487-773.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257487-774_257487-1068.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257490-255_257490-585.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257613-1119_257645-374.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257613-14_257613-541.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257613-542_257613-1118.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257645-375_257645-931.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257645-932_257682-342.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257682-343_257751-267.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257751-268_257816-98.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257816-99_257822-33.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257822-1363_257969-147.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257822-34_257822-669.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257822-670_257822-1362.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_257969-148_257969-634.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258129-30_258158-344.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258158-1572_258159-483.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258158-345_258158-925.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258158-926_258158-1571.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258159-484_258177-579.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258177-1260_258211-122.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258177-580_258177-1259.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258211-123_258287-193.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258425-3_258434-144.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258434-145_258440-320.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258440-321_258443-290.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258443-291_258448-731.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258655-60_258703-237.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258703-238_258706-344.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258706-345_258712-255.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258712-256_258742-361.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_258742-362_258750-197.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_259626-83_259685-132.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_259685-133_259686-252.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_259686-253_259721-233.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_259721-234_259721-408.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_259809-53_259821-104.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_259821-105_259862-433.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_weights/summed_259862-434_259891-108.root"

        # 2015CD history multifit        
        # # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_254790-346_254790-611.root",
        # # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_254790-612_254790-784.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_254790-90_254790-345.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256677-366_256729-326.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256729-1100_256729-1510.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256729-1511_256801-84.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256729-327_256729-708.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256729-709_256729-1099.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256801-85_256843-158.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256843-1300_256868-200.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256843-159_256843-528.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256843-529_256843-1299.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256868-202_256926-35.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_256926-36_256941-294.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257394-41_257397-56.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257397-57_257399-271.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257400-1354_257487-220.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257400-1_257400-357.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257400-358_257400-735.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257400-736_257400-1353.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257487-1069_257490-254.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257487-221_257487-492.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257487-493_257487-773.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257487-774_257487-1068.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257490-255_257490-585.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257613-1119_257645-374.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257613-14_257613-541.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257613-542_257613-1118.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257645-375_257645-931.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257645-932_257682-342.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257682-343_257751-267.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257751-268_257816-98.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257816-99_257822-33.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257822-1363_257969-147.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257822-34_257822-669.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257822-670_257822-1362.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_257969-148_257969-634.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258129-30_258158-344.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258158-1572_258159-483.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258158-345_258158-925.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258158-926_258158-1571.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258159-484_258177-579.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258177-1260_258211-122.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258177-580_258177-1259.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258211-123_258287-193.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258425-3_258434-144.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258434-145_258440-320.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258440-321_258443-290.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258443-291_258448-731.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258655-60_258703-237.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258703-238_258706-344.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258706-345_258712-255.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258712-256_258742-361.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_258742-362_258750-197.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_259626-83_259685-132.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_259685-133_259686-252.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_259686-253_259721-233.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_259721-234_259721-408.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_259809-53_259821-104.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_259821-105_259862-433.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/history_multifit/summed_259862-434_259891-108.root"

        ## weights vs multifit
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/summed_251562_251562.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415_m/summed_251562_251562.root"

        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_7415/summed_254790_254790.root",
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015C_7415_m/summed_254790_254790.root"

        ## TRANSPORT
        
        # # 2012C Boff
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/summed_200781_200798.root",
        
        # # 2015A Boff
        # "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/summed_250896_250896.root",
        
        # # 2012D Bon
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2012D_newThr/summed_208538_208686.root",
        
        # # 2015B Bon
        #"$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015B_7415/summed_251562_251562.root",
        "$CMSSW_BASE/src/PhiSym/EcalCalibAlgos/ntuples/2015D_7415/summed_256843_256843.root"
    ),

    correctionsFiles = cms.vstring(
        # #---Boff
        # "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2012C_v1/corrections_200781_200798.txt",
        # "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_12_patch4/src/PhiSym/EcalCalibAlgos/ntuples/2015A_v2/corrections_250896_250896.txt"

        #---Bon
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_15/src/PhiSym/EcalCalibAlgos/ntuples/2012D_newThr/corrections_208538_208686.txt",
        # "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_15/src/PhiSym/EcalCalibAlgos/corrections_208538_208686.txt",
        # "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_15/src/PhiSym/EcalCalibAlgos/ntuples/2015B_newGT_v5/corrections_251562_251562.txt"
        "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_15/src/PhiSym/EcalCalibAlgos/ntuples/2015D_7415/corrections_256843_256843.txt"
        # "/afs/cern.ch/user/s/spigazzi/work/ECAL/CMSSW_7_4_15/src/PhiSym/EcalCalibAlgos/corrections_256843_256843.txt"        
    )
)

