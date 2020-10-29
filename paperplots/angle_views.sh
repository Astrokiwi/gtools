for j in 15 30 45 60 75 90
do
    python sph_anim.py 3032 $1 --rad 20. --phi $j --savemap --view face --snap0 200 --maxsnapf 200 --plot view8mic,view12mic,view850mic --L 100 --suffix _theta=$j &
done
