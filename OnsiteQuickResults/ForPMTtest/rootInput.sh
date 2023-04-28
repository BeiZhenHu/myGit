gcu=$1
date=$2
path=/home/run/SPMT_zhangsh/SPMT_Readout_zhangsh/data/OnsiteQuickResults/0_bin2root/rootFile/

rm datalist
for line in $(ls $path*$1*PED*$2*)
do
	line=${line##*/}
	echo ${line%.*} >> datalist
done
for line in $(ls $path*$1*PM*$2*)
do
	line=${line##*/}
	echo ${line%.*} >> datalist
done
