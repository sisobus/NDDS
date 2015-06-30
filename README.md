# NDDS

## utils.py
여러 스크립트에서 공통적으로 사용하는 모듈들을 모아놨습니다.

- getOptions() : arguments를 options로 바꿔 반환해줍니다. arguments의 옵션을 추가하고 싶으면 여기에 추가하면 됩니다.
- createDirectory(directoryFileName) : 말그대로 디렉토리를 생성해줍니다.
- saveGraph(imageFileName,xp,yp) : x에 대한 y값 그래프를 그려서 png image로 저장해줍니다.
- getDataInFile(filename) : gdp perl script의 데이터 파일 형식이 파싱을 필요로 하기 때문에 모든 데이터는 이와 같은 형식을 갖고 있습니다. 이러한 데이터 파일을 파싱해서 2차원 리스트로 반환해줍니다.
- hammingDistance(a,b) : 두 포인트의 해밍 거리를 구해줍니다.
- writeDataToFile(filename,datas) : 2차원 리스트를 파일에 출력합니다.
- readDataFromFile(filename) : 파일로부터 2차원 리스트를 로드합니다.
- calculateCorrelationCoefficient(vp1,vp2,datas) : 두 빈티지 포인트의 Correlation Coefficient를 구해줍니다.


## generate_data.py
데이터를 생성해주는 스크립트입니다.

- -n : number of data 
- -m : number of vantage point
- -d : number of dimension
- -a : number of alphabet
- -b : distribution [ u : uniform, c1 : skewed, c10 : 10-cluster ]
- -v : type of vantage point [ greedy, random, AA, AB, etc ... ]
- data directory에 있는 gdp.pl perl script를 이용합니다.

ex) python generate_data.py -n 100000 -d 25 -a 4 -b u

## generate_vantage_point.py
빈티지 포인트를 생성하는 스크립트입니다.

- 사용자가 원하는 데로 수정해서 사용해야합니다.
- generateVantagePoints(options) 를 수정해서 빈티지 포인트를 뽑을 수있습니다.
- generateGreedyVantagePoints(options) 로 2차년도 과제의 greedy 알고리즘을 이용하여 빈티지 포인트를 생성합니다. |A|가 4일 때, 10차원을 초과하면 굉장히 늦습니다.

## generate_query.py
query를 생성하는 스크립트입니다.

- data directory에 존재하는 모든 데이터에 대해서 100개의 랜덤 query points를 생성합니다.

## calculate_cc.py
Correlation Coefficient를 계산해주는 스크립트입니다.

- 옵션을 주어 데이터 파일과 빈티지 포인트 디렉토리에 있는 것으로 Correlation Coefficient를 구해줍니다.

## convert_ndds_to_cds.py
빈티지 포인트를 이용하여 NDDS 데이터 포인트와 쿼리 포인트를 cds 로 변환해줍니다.

- cds_data, cds_query directory에 저장됩니다.

ex) python convert_ndds_to_cds.py -n 100000 -d 10 -m 10 -a 4 -b u -v greedy
