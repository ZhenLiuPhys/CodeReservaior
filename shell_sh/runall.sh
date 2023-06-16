for i in {1..24}
do
echo "$i"
cd work"$i"
./main mssm3.par "$i" output"$i".dat &
cd ..
done
