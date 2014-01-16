LIST=$1
OUTPUT=$2

rm -f $OUTPUT
echo "{}" > $OUTPUT
for l in `cat $LIST`
do
./combine_JSON.py -a $OUTPUT -b $l -o tmp -r or
mv tmp $OUTPUT
done
