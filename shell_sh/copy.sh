for i in {1..24}
do
echo "$1"
mkdir work"$i"
cp -r work/* work"$i"/.
cd work"$i"/
make main=main.F
cd ..
done
