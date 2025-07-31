# ERNNplus-code
g++ -O3 -g -std=c++17 main.cpp -o run

./run [method] [dataset] [num of facility vertex] [num of round] [number budget] [total weight budget]

method:
DBEISTAR vldb2024 (baseline)
SPTgreedy our purposed method SPTG
SPT   our purposed method SPT
SPT1 our purposed method SPT-Fast

Example：

./run SPT CT 100 1 4 5000

./run DBEISTAR CT 100 1 4 5000
