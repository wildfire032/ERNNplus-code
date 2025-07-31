#include "graph.h"

int returnNumOfVertices(string filename) {
    Graph g(filename.c_str());
    return g.nVertices;
}

vector<int> returndeglist(string filename) {
    Graph g(filename.c_str());
    return g.degreelist;
}

void testDynamic(string filename, vector<int> poi) {
    Graph g(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();
}

pair<int, double> Exact(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.exact(targetID, budget);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout << "Exact: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    return make_pair(gain, duration.count() * 1.0 / 1000);
}

pair<int, double> Basic(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.greedy(targetID, budget, Prune::NoPrune);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout << "Basic: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    return make_pair(gain, duration.count() * 1.0 / 1000);
}

pair<int, double> DBEI(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.greedy(targetID, budget, Prune::Prune3);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout << "DBEI: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    return make_pair(gain, duration.count() * 1.0 / 1000); 
}

pair<int, double> DBEIPLUS(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.greedy(targetID, budget, Prune::Prune1 | Prune::Prune3);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout << "DBEIPLUS: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    return make_pair(gain, duration.count() * 1.0 / 1000); 
}

pair<int, double> DBEISTAR(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    Graph orig(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.greedy(targetID, budget,
                        Prune::Prune1 | Prune::Prune2 | Prune::Prune3);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout << "DBEISTAR: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
        //cout << "DBEIPLUS: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    //cout<<check(g,orig,budget.first,budget.second)<<"\n";
    //assert(check(g,orig,budget.first,budget.second));
    return make_pair(gain, duration.count() * 1.0 / 1000);
}
pair<int, double> potential(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.convert_potential(targetID, budget, 0);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout << "DBEIPLUS: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    return make_pair(gain, duration.count() * 1.0 / 1000); 
}
pair<int, double> SPT(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    Graph orig(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.shortestpathtree_method0(targetID, budget, 0);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout << "DBEIPLUS: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    //cout<<check(g,orig,budget.first,budget.second)<<"\n";
    //assert(check(g,orig,budget.first,budget.second));
    return make_pair(gain, duration.count() * 1.0 / 1000); 
}
pair<int, double> SPT1(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    Graph orig(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.shortestpathtree_method(targetID, budget, 0);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout<<check(g,orig,budget.first,budget.second)<<"\n";
    //cout << "DBEIPLUS: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    assert(check(g,orig,budget.first,budget.second));
    return make_pair(gain, duration.count() * 1.0 / 1000); 
}
pair<int, double> SPTgreedy(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    Graph orig(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.shortestpathtree_methodgreedy(targetID, budget, 0);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout<<check(g,orig,budget.first,budget.second)<<"\n";
    //cout << "DBEIPLUS: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    //assert(check(g,orig,budget.first,budget.second));
    return make_pair(gain, duration.count() * 1.0 / 1000); 
}
pair<int, double> Weight(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    g.setPOI(poi);

    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.maxWeight(targetID, budget);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout << "Weight: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    return make_pair(gain, duration.count() * 1.0 / 1000);
}

pair<int, double> Neighbor(string filename, vector<int> poi, int targetID, pair<int,int> budget) {
    Graph g(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();

    auto start = chrono::high_resolution_clock::now();
    int gain = g.chooseNbr(targetID, budget);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    //cout << "Neighbor: Gain = " << gain
    //     << ", Time = " << duration.count() * 1.0 / 1000 << endl;
    return make_pair(gain, duration.count() * 1.0 / 1000);
}

int main(int argc, char *argv[]) {
    string method = argv[1];
    string strArg = argv[2];
    string filename = "data/" + strArg + ".tmp"; //"data/NY.tmp";
    int n = returnNumOfVertices(filename);

    int numOfPOI = std::atoi(argv[3]);
    int numberOfTest = std::atoi(argv[4]);
    int budget_number = std::atoi(argv[5]);
    long long budget_weight =std::atoi(argv[6]);
    pair<int,long long> budget ={budget_number,budget_weight};
    vector<int> poi;
    poi.clear();
    srand(1234); // Fixed seed for consistency
    for (int i = 0; i < numOfPOI; i++) {
        int randomNum = rand() % n;
        poi.push_back(randomNum);
    }
    // vector<int> degreelist = returndeglist(filename);
    // vector<int> tmp;
    // for(int i=0;i<n;++i)
    // {
    //     for(int j=1;j<=degreelist[i];++j) tmp.push_back(i);
    // }
    //     for (int i = 0; i < numOfPOI; i++) {
    //     int randomNum = tmp[rand() % tmp.size()];
    //     poi.push_back(randomNum);
    // }
    
    vector<int> vid(numOfPOI);
    iota(vid.begin(), vid.end(), 0);

    default_random_engine engine(1234);
    shuffle(vid.begin(), vid.end(), engine);
    double TotalGain = 0;
    double TotalTime = 0;
    for (int i = 0; i < numberOfTest; ++i) {
        int id = vid[i];
        int targetID = poi[id];

        //cout << "Round : " << i+1 << ", targetID = " << targetID+1 << endl;
        if(method == "Weight"){
            pair<int, double> ans = Weight(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else if(method == "Neighbor"){
            pair<int, double> ans = Neighbor(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else if(method == "DBEISTAR"){
            pair<int, double> ans = DBEISTAR(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else if(method == "DBEIPLUS"){
            pair<int, double> ans = DBEIPLUS(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else if(method == "DBEI"){
            pair<int, double> ans = DBEI(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else if(method == "Basic"){
            pair<int, double> ans = Basic(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else if(method == "Exact"){
            pair<int, double> ans = Exact(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else if(method == "SPTgreedy"){
            pair<int, double> ans = SPTgreedy(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else if(method == "SPT"){
            pair<int, double> ans = SPT(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else if(method == "SPT1"){
            pair<int, double> ans = SPT1(filename, poi, targetID, budget);
            //cout << ans.first << " " << ans.second << endl;
            if(ans.first <= 0 or ans.first >= 1e9) ans.first = 0;
            TotalGain += ans.first;
            TotalTime += ans.second;
        }
        else{
            cout << "Parameter Wrong" << endl;
            return 0;
        }
    }
    cout << "Average Gain = " << 1.0 * TotalGain / numberOfTest << ", "
         << "Average Time = " << 1.0 * TotalTime / numberOfTest << endl;
    return 0;
}
