echo "20d"
python simulate_nd_tree_knn.py -n 100000 -d 20 -b u -a 10 > result/nd_knn_100000_20_u.out
echo "40d"
python simulate_nd_tree_knn.py -n 100000 -d 40 -b u -a 10 > result/nd_knn_100000_40_u.out
echo "60d"
python simulate_nd_tree_knn.py -n 100000 -d 60 -b u -a 10 > result/nd_knn_100000_60_u.out
echo "80d"
python simulate_nd_tree_knn.py -n 100000 -d 80 -b u -a 10 > result/nd_knn_100000_80_u.out
echo "100d"
python simulate_nd_tree_knn.py -n 100000 -d 100 -b u -a 10 > result/nd_knn_100000_100_u.out
echo "120d"
python simulate_nd_tree_knn.py -n 100000 -d 120 -b u -a 10 > result/nd_knn_100000_120_u.out
