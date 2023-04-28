abc=$1
t=$2

PED=$(head -n 1 ./datalist)
echo $PED
for line in $(grep "PM" ./datalist)
do
	PM=${PM},$line
done
PM=${PM#*,}
echo $PED
echo $PM
echo $abc
echo $t
./main_wudr /home/run/SPMT_zhangsh/SPMT_Readout_zhangsh/data/OnsiteQuickResults/0_bin2root/rootFile/ $PED $PM $abc $t
