#!/bin/bash
commFile=$1 #community profile
N=$2
cancerType=$3

begin=$(date +%S)

mkdir $cancerType.entrez

python sym2entrez_v2.py sym2entrez_dir/$cancerType.* ../$cancerType.top$N.community_detection_mlc $cancerType.entrez/$cancerType.top$N.entrez
python myGOenrichment.py $cancerType.entrez/$cancerType.top$N.entrez -trait2genes GOBPname2genes.human.BP.txt -pcut 0.05 -topK 1 > ../$cancerType.top$N.GO

end=$(date +%S)
echo $((($end-$begin)/3600)):$((($end-$begin)/60)):$(($end-$begin))	
