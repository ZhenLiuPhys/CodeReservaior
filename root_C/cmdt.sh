for i in {72..80}
do
cp ../Events/run_$i/tag_1_pythia_lhe_events.root .
root -b -l Sig.C <<EOF
Sig t
t.Loop()
EOF
cp Sigoutput.root Sigoutput$i.root
done
