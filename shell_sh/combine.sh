mkdir data
for i in {1..24}
do
mv work"$i"/output"$i".dat data/.
done
cd data/
cat output*.dat > all.dat
