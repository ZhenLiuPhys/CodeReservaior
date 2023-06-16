j=51
k=9
for i in $(cat widthtable2.dat); do
echo $i
j=`expr $j + 1`
k=`expr $k + 1`
echo $i
cd run_$j
root -l unweighted_events.root <<EOF
LHEF->MakeClass("s8")
EOF
cp ../s8.C .
root -l s8.C <<EOF
s8 t
t.Loop()
EOF
root -l output.root <<EOF
cout << "result" << '\t' << `echo "$k * 100"` << '\t' << `echo "$k * 100 - $i * 3 / 2"` << '\t' << `echo $i` << '\t' << mytree->GetEntries() << '\t' << mytree->GetEntries("TMath::Abs(eta1)<=2.5&&TMath::Abs(eta2)<=2.5") << '\t' << mytree->GetEntries("TMath::Abs(eta1)<=2.5&&TMath::Abs(eta2)<=2.5&&TMath::Abs(deta)<=1.3") << '\t' << mytree->GetEntries("TMath::Abs(eta1)<=2.5&&TMath::Abs(eta2)<=2.5&&TMath::Abs(deta)<=1.3&&minv>`echo "$k * 100 - $i * 3 / 2"`") << '\n';
EOF
cd ..
done
