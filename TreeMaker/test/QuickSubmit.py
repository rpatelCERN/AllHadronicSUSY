import sys
import os
#flist=open("%s" %sys.argv[1], 'r')
#print "%s" %sys.argv[1]
count=0

#for line in flist:
#	print line
#	rline=line.rstrip('\n')
for i in range(0, 110):
#for line in flist:
		#inc=(i-3)*50+1000
		#count=count+1
		skip=(i*1000)
		f=open("bsub%d.sh" %i, 'w')
		f.write("set SCRAM_ARCH=slc6_amd64_gcc481\n")
		f.write("cd /afs/cern.ch/user/r/rpatel/CMSSW_7_2_3_patch1/src\n")
		f.write("eval `scramv1 runtime -sh`\n")
		f.write("cd AllHadronicSUSY/TreeMaker/test\n")
		f.write("cmsRun runMakeTreeFromMiniAOD_cfg.py  %d skip=%d \n" %(i,skip))
#		f.write("rm -f core.*")
		#f.write("cmsRun step2_DIGI_L1_L1TrackTrigger_DIGI2RAW_RECO_PU.py %s %d %d" %(rline, count,inc))
		f.close()
	#	os.system("cmsRun test_L1PixelTrack.py  %d %d \n" %(i,skip))
		os.system("chmod 744 bsub%d.sh" %i)
		os.system("bsub -R \"pool>30000\" -q 1nh -J job%d < bsub%d.sh" %(i,i))
