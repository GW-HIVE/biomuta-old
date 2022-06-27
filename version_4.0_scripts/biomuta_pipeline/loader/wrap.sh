#### Step-1
#Load data to database tables following this order:

#nohup python do-loader.py -i ../conf/config-4.1.json &

#nohup python mut-loader.py -i ../conf/config-4.1.json &

#for i in 'civic'  'clinvar'  'cosmic' 'icgc'  'tcga';
#do
#	python freq-loader.py -i ../conf/config-4.1.json -s $i 
#done 

#nohup python eff-loader.py -i ../conf/config-4.1.json &
#nohup python protein-loader.py -i ../conf/config-4.1.json &
nohup python pph-loader.py -i ../conf/config-4.1.json &
#nohup python netnglyc-loader.py -i ../conf/config-4.1.json &
#nohup python uniprotann-loader.py -i ../conf/config-4.1.json &
#nohup python stat-loader.py -i ../conf/config-4.1.json &


#delete-records.py
#remove-freq-records.py
#dump-freq-stat.py


