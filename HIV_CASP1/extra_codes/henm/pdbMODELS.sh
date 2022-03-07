cpdb=$1

echo 'running awk to convert each frame to a seperate models'
echo 'using input file: ' $cpdb

awk 'BEGIN{n=1}; /CRYST/ {print "MODEL "n; n++}; ! /CRYST/ && ! /END/ {print $0}; /END/ {print "ENDMDL\nMODEL "n; n++}' $cpdb > temp.pdb

echo 'removing tailing MODEL keyword and overwriting'
sed '$d' < temp.pdb > $cpdb

rm temp.pdb
