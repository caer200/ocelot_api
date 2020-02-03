for i in *.cif
do
sed '/_shelx_res_file/q' $i >> $i
done
