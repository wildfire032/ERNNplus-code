# ERNNplus-code
g++ -O3 -g -std=c++17 main.cpp -o run

./run [method] [dataset] [num of facility vertex] [num of round] [number budget] [total weight budget]

method:

DBEISTAR : DBEI* (vldb2024 baseline)

SPTgreedy : our purposed method SPTG

SPT   : our purposed method SPT

SPT1 : our purposed method SPT-Fast

Example：

./run SPT CT 1000 50 100 10000

./run DBEISTAR CT 1000 50 200 50000
