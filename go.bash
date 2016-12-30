#for((i=20;i<=120;i+=20)) do
#    echo "$i"
#    python simulate_mvp_tree_auto_range.py -n 100000 -d "$i" -b u -a 10 > result/mvp_100000_"$i"_u_10.out
#done
#for((i=20;i<=120;i+=20)) do
#    echo "$i"
#    python simulate_mvp_tree_auto_range.py -n 100000 -d "$i" -b 10 -a 10 > result/mvp_100000_"$i"_10_10.out
#done
for((i=40;i<=120;i+=20)) do
    echo "$i"
    python simulate_mvp_tree_auto_range.py -n 100000 -d "$i"  -b 1 -a 10 > result/mvp_100000_"$i"_1_10.out
done
