# NDDS

## utils.py
여러 스크립트에서 공통적으로 사용하는 모듈들을 모아놨습니다.

- getOptions() : arguments를 options로 바꿔 반환해줍니다. arguments의 옵션을 추가하고 싶으면 여기에 추가하면 됩니다.
  - -n : number of data 
  - -m : number of vantage point
  - -d : number of dimension
  - -a : number of alphabet
  - -b : distribution [ u : uniform, c1 : skewed, c10 : 10-cluster ]
  - -v : type of vantage point [ greedy, random, AA, AB, etc ... ]
  - -r : query range
- createDirectory(directoryFileName) : 말그대로 디렉토리를 생성해줍니다.
- saveGraph(imageFileName,xp,yp) : x에 대한 y값 그래프를 그려서 png image로 저장해줍니다.
- saveGraphWithHighValue(imageFileName,xp,yp,highValue) : x에 대한 y값 그래프와 최고점의 값을 그려서 png image로 저장해줍니다.
- saveGraphUsingPointWithCC(imageFileName,xp,yp,cc,dim) : x에 대한 y값 그래프와 correlation coefficient를 그려서 png image로 저장해줍니다.
- saveGraphUsing3DSurfaceWithCC(imageFileName,xp,yp,zp,cc,dim) : x, y, z에 대한 그래프와 correlation coefficient를 그려서 png image로 저장해줍니다.
- getDataInFile(filename) : gdp perl script의 데이터 파일 형식이 파싱을 필요로 하기 때문에 모든 데이터는 이와 같은 형식을 갖고 있습니다. 이러한 데이터 파일을 파싱해서 2차원 리스트로 반환해줍니다.
- hammingDistance(a,b) : 두 포인트의 해밍 거리를 구해줍니다.
- writeDataToFile(filename,datas) : 2차원 리스트를 파일에 출력합니다.
- readDataFromFile(filename) : 파일로부터 2차원 리스트를 로드합니다.
- calculateCorrelationCoefficient(vp1,vp2,datas) : 두 빈티지 포인트의 Correlation Coefficient를 구해줍니다.
- getDataFileName(options) : argument options을 이용하여 데이터 파일 이름을 반환해줍니다.
- getQueryFileName(options) : argument options을 이용하여 Query 파일 이름을 반환해줍니다.
- getVPFileName(options) : argument options을 이용하여 빈티지 포인트 파일 이름을 반환해줍니다.
- getCDSDataFileName(options) : argument options을 이용하여 cds 데이터 파일 이름을 반환해줍니다.
- getCDSQueryFileName(options) : argument options을 이용하여 cds Query 파일 이름을 반환해줍니다.
- getImageFileName(options,tag) : argument options을 이용하여 이미지 파일 이름을 반환해줍니다.
- getNDTDataFileName(options) : argument options을 이용하여 nd_tree 버전 데이터 파일 이름을 반환해줍니다.
- getNDTQueryFileName(options) : argument optios을 이용하여 nd_tree 버전 Query 파일 이름을 반환해줍니다.
- getRQResultFileName(options) : argument options를 이용하여 range query result filename을 반환해줍니다.
- getFigurePairName(options,id1,id2) : argument options를 이용하여 figure pair filename을 반환해줍니다.
- executeCommand(command) : commands 라이브러리를 이용하여 command를 수행합니다.


## generate_data.py
데이터를 생성해주는 스크립트입니다.

- data directory에 있는 gdp.pl perl script를 이용합니다.

ex) python generate_data.py -n 100000 -d 25 -a 4 -b u

## generate_vantage_points.py
빈티지 포인트를 생성하는 스크립트입니다.

- 사용자가 원하는 데로 수정해서 사용해야합니다.
- generateVantagePoints(options) 를 수정해서 빈티지 포인트를 뽑을 수있습니다.
- generateGreedyVantagePoints(options) 로 2차년도 과제의 greedy 알고리즘을 이용하여 빈티지 포인트를 생성합니다. |A|가 4일 때, 10차원을 초과하면 굉장히 늦습니다.

ex) python generate_vantage_points.py -d 10 -m 10 -a 4 

## generate_query.py
query를 생성하는 스크립트입니다.

- data directory에 존재하는 모든 데이터에 대해서 100개의 랜덤 query points를 생성합니다.

ex) python generate_query.py

## calculate_cc.py
Correlation Coefficient를 계산해주는 스크립트입니다.

- 옵션을 주어 데이터 파일과 빈티지 포인트 파일에 있는 것으로 Correlation Coefficient를 구해줍니다.

ex) python calculate_cc.py -n 1000 -d 10 -a 2 -b u -m 2 -v greedy

## convert_ndds_to_cds.py
빈티지 포인트를 이용하여 NDDS 데이터 포인트와 쿼리 포인트를 cds 로 변환해줍니다.

- cds_data, cds_query directory에 저장됩니다.

ex) python convert_ndds_to_cds.py -n 100000 -d 10 -m 10 -a 4 -b u -v greedy

## draw_deep_graph.py
빈티지 포인트들을 이용하여 데이터 포인트가 얼마나 퍼지는지 그래프를 그려주는 스크립트입니다. 

- 모든 상태를 보여줄 수 없으니, 각 빈티지 포인트에 의해 매핑된 좌표 중 가장 많이 나타난 좌표에 매핑된 포인트를 기준으로 다음 빈티지 포인트로 매핑해줍니다.
- 얼마만큼 줄어들었느냐가 퍼짐의 척도가 됩니다.

ex) python draw_deep_graph.py -n 100000 -d 10 -m 10 -a 4 -b u -v greedy

## convert_data_for_ndt.py
NDDS 데이터 파일과 Query 파일 형식을 nd_tree의 파일 형식에 맞게 변경해주는 스크립트입니다.

- nd_tree의 데이터 파일과 range query, box query 파일의 형식은 정수형으로 이루어져 있습니다.
  - A = 0, B = 1, C = 2, D = 3, ...
  - 아마 ND-tree가 그룹을 비트 마스트를 이용하여 쉽게 하기 위해서 인듯 합니다.
  - 이를 이용하면 알파벳 개수에 제한 없이 그루핑을 할 수 있는 것 같습니다.
- data/*.txt 와 query/*.txt의 데이터 형식을 변형하여 ndt_data 와 ndt_query 디렉토리로 옮겨줍니다.

ex) python convert_data_for_ndt.py

## simulate_nd_tree.py
nd_tree 디렉토리의 NDT 파일을 실행해주는 스크립트입니다.

- NDT 프로그램은 특이하게 config.h 파일 안에 차원이나 데이터의 개수 등등을 const int로 정의하였습니다.
- 이를 수정하기보다, 기존의 config.h 파일을 original_config.h로 옮겨두고, NDT를 실행할 때마다 옵션 값을 수정하여 config.h를 새로 생성하도록 만들었습니다.
- NDT 프로그램은 데이터의 파일 형식도 sourceData%d+0.txt, rangequeryAll.txt, boxqueryAll.txt 로 고정되어 있으므로 ndt_data, ndt_query 디렉토리 하위에 있는 파일을 nd_tree디렉토리로 옮기며 이름을 위와 같이 변경해주도록 구현했습니다.

- 이 스크립트의 수행 flow는 다음과 같습니다.
  - options을 argument로 받습니다.
  - original_config.h 파일을 options에 맞게 수정하여 config.h 파일을 생성합니다.
  - ndt_data 디렉토리 하위에 있는 options에 맞는 데이터 파일을 nd_tree/sourceData%d+0.txt 로 옮깁니다.
  - ndt_query 디렉토리 하위에 있는 options에 맞는 Query 파일을 nd_tree/rangequeryAll.txt, nd_tree/boxqueryAll.txt 로 옮깁니다.
  - nd_tree 디렉토리의 make를 수행합니다.
  - NDT 프로그램을 실행합니다.

ex) python simulate_nd_tree.py -n 100000 -d 10 -b u -a 4

## simulate_kd_tree.py
kd_tree 디렉토리의 kdtree 파일을 실행해주는 스크립트입니다.

- kdtree 프로그램은 boost_1_55_0 라이브러리의 설치를 필요로 합니다.
- http://www.boost.org/users/history/version_1_55_0.html 에서 받을 수 있고, /usr/local/ 에 설치하면됩니다.

- 이 스크립트의 수행 flow는 다음과 같습니다.
  - cds_data, cds_query 디렉토리 안의 데이터와 query포인트를 이용합니다.
  - simulate_nd_tree.py 스크립트와 유사하게 cds_data, cds_query 디렉토리 안의 데이터와 query 파일을 kd_tree/data.txt, kd_tree/query.txt 로 카피합니다.
  - make를 수행하고, cd kd_tree ; ./kdtree -load_file data.txt -orig_dim %d -dim %d -rqfile query.txt -range %d -count %d 를 수행합니다.
  - range는 -r로 줍니다.

ex) python simulate_kd_tree.py -n 100000 -d 10 -m 10 -r 1 -v greedy

## simulate_kd_tree_auto_range.py
simulate_kd_tree.py 의 queryRange를 1부터 dimension-1까지 실행해주는 스크립트입니다.

ex) python simulate_kd_tree_auto_range.py -n 100000 -d 10 -m 10 -v greedy

## generate_rq_result.py
data file과 query file로 linear한 range query 결과를 생성해주는 스크립트입니다.

- 당연히 true negative는 존재하지 않습니다.

ex) python generate_rq_result.py -n 100000 -d 10 -a 4 -b u

## draw_pair_graph.py
빈티지 포인트 집합안에서 모든 pair의 공간변환 후 그래프를 그려줍니다.

- 2d 그래프로 그리면 얼만큼 뭉쳐있는지를 알 수 없어서 3d로 그립니다.
- z 축은 x,y 에 매핑된 데이터 포인트의 수를 나타냅니다.

ex) python draw_pair_graph.py -n 100000 -d 10 -m 10 -a 4 -b u -v random

## calculate_linear_time.cpp
선형 Query time을 재는 코드입니다.

- data directory안에 있는 모든 데이터 파일과 query directory안에 있는 모든 query 파일에 대한 결과를 출력합니다.
- 결과는 result/result_*.txt에 저장됩니다.

ex) g++ calculate_linear_time.cpp , ./a.out

## generate_vantage_point_using_many.py

- python generate_vantage_point_using_many.py -n 100000 -m 5 -d 40 -a 10 -b 1 -v minor

## calculate_all_pair_distance.py
모든 페어 거리를 보여주는 코드입니다.

- python calculate_all_pair_distance.py -n 100000 -m 5 -d 40 -a 10 -b 1 -v minor
