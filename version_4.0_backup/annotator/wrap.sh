#### Step-1
#for i in `seq 1 22` 'X' 'Y';
#do
#	nohup python pph-step1.py -i ../conf/config-4.1.json -c $i &
#done 

#### Step-2
#nohup python pph-step2.py -i ../conf/config-4.1.json &


#### Step-3
for i in `seq 1 5` 'X' 'Y';
do
       nohup python pph-step3.py -i ../conf/config-4.1.json -c $i &
done 


#### Step-4
#for i in `seq 1 22` 'X' 'Y';
#do
#	inFile="outdir2/pph-features-step3.chr$i"".txt"
#	outFile="outdir2/pph-predictions.chr$i"".txt"
#	logFile="outdir2/pph-step4.chr$i"".log"
#	/softwares/polyphen-2.2.2/bin/run_weka.pl $inFile 1> $outFile  2> $logFile
#done 


#### Step-5
#nohup python pph-step5.py -i ../conf/config-4.1.json &


#### Step-6
#for i in `seq 1 22` 'X' 'Y';
#do
#       nohup python netnglyc-step1.py -i ../conf/config-4.1.json -c $i &
#done 


#### Step-7
#nohup python netnglyc-step2.py -i ../conf/config-4.1.json &




