#for((i=20;i<=120;i+=20)) do
#    echo "$i"
#    python simulate_mvp_tree_auto_range.py -n 100000 -d "$i" -b u -a 10 > result/mvp_100000_"$i"_u_10.out
#done
#for((i=20;i<=120;i+=20)) do
#    echo "$i"
#    python simulate_mvp_tree_auto_range.py -n 100000 -d "$i" -b 10 -a 10 > result/mvp_100000_"$i"_10_10.out
#done
for((i=20;i<=120;i+=20)) do
    echo "$i"
    python simulate_kd_tree_knn.py -n 100000 -a 10 -b u -v ppgeh -m 15 -k 1 -d "$i" > knn_result/kd_100000_"$i"_u_10_1.out
    python simulate_kd_tree_knn.py -n 100000 -a 10 -b u -v ppgeh -m 15 -k 5 -d "$i" > knn_result/kd_100000_"$i"_u_10_5.out
    python simulate_kd_tree_knn.py -n 100000 -a 10 -b u -v ppgeh -m 15 -k 10 -d "$i" > knn_result/kd_100000_"$i"_u_10_10.out
done
