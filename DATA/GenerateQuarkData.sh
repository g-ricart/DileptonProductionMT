cat ECurveQCD.txt | awk '{print $1}' > wT.list

rm QCurveQCD.txt

for wT in `cat wT.list`
do

cat EvolutionParameters_ee.txt | awk -v wVal=${wT} 'BEGIN{wLow=0.0; wHigh=0.0; QLow=0.0; QHigh=0.0;} {wLow=wHigh; QLow=QHigh; wHigh=$1; QHigh=$3/$2; if($1>wVal){print wVal,((wHigh-wVal)*QHigh+(wVal-wLow)*QLow)/(wHigh-wLow); exit;}}' >> QCurveQCD.txt

done

rm wT.list
