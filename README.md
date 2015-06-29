# NDDS

## utils.py
여러 스크립트에서 공통적으로 사용하는 모듈들을 모아놨습니다.

## generate_data.py
데이터를 생성해주는 스크립트입니다.
-n : number of data 
-m : number of vantage point
-d : number of dimension
-a : number of alphabet
-b : distribution [ u : uniform, c1 : skewed, c10 : 10-cluster ]

data directory에 있는 gdp.pl perl script를 이용합니다.
ex) python generate_data.py -n 100000 -d 25 -a 4 -b u

## generate_vantage_point.py
빈티지 포인트를 생성하는 스크립트입니다.

사용자가 원하는 데로 수정해서 사용해야합니다.
generateVantagePoints(options) 를 수정해서 빈티지 포인트를 뽑을 수있습니다.
generateGreedyVantagePoints(options) 로 2차년도 과제의 greedy 알고리즘을 이용하여 빈티지 포인트를 생성합니다. 10차원을 초과하면 굉장히 늦습니다.

## generate_query.py
query를 생성하는 스크립트입니다.

data directory에 존재하는 모든 데이터에 대해서 100개의 랜덤 query points를 생성합니다.

## calculate_cc.py
Correlation Coefficient를 계산해주는 스크립트입니다.

옵션을 주어 데이터 파일과 빈티지 포인트 디렉토리에 있는 것으로 구해줍니다.
