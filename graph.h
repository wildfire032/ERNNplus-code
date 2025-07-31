#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <thread>
#include <vector>
#include <cassert>

//#define DEBUG

using namespace std;

const long long INF = 1e18;
const long long intINF = 1e9;
typedef pair<int, int> pii;

enum Prune { NoPrune = 0, Prune1 = 1, Prune2 = 2, Prune3 = 4 };

class Vertex {
  public:
    vector<pair<int, int>> adj;
    int id;
    int isPOI;
    int deg;
    long long dist;
    int closest_node;

    bool operator<(const Vertex &rhs) const { return dist > rhs.dist; }
  public:
    int ispot;
    int father_node;//father in shortest path from node target
    long long target_dis;
    int target_edge;// zuiduanlu shang zuishaobianshu


};

class Graph {
  public:
    vector<Vertex> graph;
    int nVertices;
    int maxedgeweight=0;
    map<pair<int,int>,int> edge;
    vector<int> degreelist;

  public:
    Graph(string path);
    ~Graph();

  public:
    void dijkstra();
    int updateWeight(int u, int v, int delta);
    void dynamicDijkstra(int u, int v, int delta);
    void recomputeDijkstra(int u, int v, int delta);
    void dijkstra_onetarget(int target);

  public:
    int exact(int targetID, pair<int,int> budget);
    int greedy(int targetID, pair<int,long long> budget, int pruneOP);
    int convert_potential(int targetID, pair<int,int> budget, int pruneOP);

  public:
    vector<int> boundList;
    void commputeBound(int targetID);

  public:
    int maxWeight(int targetID, pair<int,int> budget);
    int chooseNbr(int targetID, pair<int,int> budget);

  public:
    void setPOI(vector<int> poi);

  public:
  int shortestpathtree_methodgreedy(int targetID, pair<int,long long> budget, int pruneOP);// caclulate true increment
  int shortestpathtree_method0(int targetID, pair<int,long long> budget, int pruneOP);//use lca to calculate the convert value (by definition)
    int shortestpathtree_method(int targetID, pair<int,long long> budget, int pruneOP);// calculate by jump one by one
    int shortestpathtree_method2(int targetID, pair<int,long long> budget, int pruneOP);//binary lifting speed up
    //bool check(Graph now,Graph ori,int budgetn,int budgetw);
    // Graph Origin;
};

// Implementation
bool check(Graph now,Graph ori,int budgetn,long long budgetw){
    int n=ori.nVertices;
    int nown=0;long long noww=0;
    for(int i=0;i<n;++i)
    {
        for(int j=0;j<now.graph[i].adj.size();j++)
        {
            auto u=now.graph[i].adj[j],v=ori.graph[i].adj[j];
            assert(u.first==v.first);
            if(u.first>i){
                if(u.second!=v.second){
                    nown++;noww+=abs(u.second-v.second);
                }
            }
        }
    }
    //cout<<nown<<" "<<noww<<"\n";
    if(nown<=budgetn&&noww<=budgetw) return true;
    return false;
}
Graph::Graph(string filename) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Could not open file " << filename << endl;
        return;
    }

    int n, m;
    infile >> n >> m;

    nVertices = n;
    graph.resize(n);
    int u, v, w;
    set<pair<int, int>> uvPairs;

    for (int i = 0; i < m; i++) {
        infile >> u >> v >> w; // read edge
        u--, v--;              // starting from one
        maxedgeweight=max(maxedgeweight,w);
        auto [it, inserted] = uvPairs.insert(make_pair(u, v));
        if (inserted) {
            graph[u].adj.push_back(make_pair(v, w));
        }
        auto [it2, inserted2] = uvPairs.insert(make_pair(v, u));
        if (inserted2) {
            graph[v].adj.push_back(make_pair(u, w));
        }
        edge[{u,v}]=w;
        edge[{v,u}]=w;
    }
    //cout<<maxedgeweight<<" \n";
    infile.close();

    for (int i = 0; i < n; i++) {
        graph[i].id = i;
        graph[i].isPOI = 0;
        graph[i].dist = INF;
        graph[i].closest_node = -1;
        graph[i].deg = graph[i].adj.size();
        degreelist.push_back(graph[i].deg);
    }
    // cout << n << " " << m << endl;
}

Graph::~Graph() {}

void Graph::dijkstra() {
    for (int i = 0; i < nVertices; i++) {
        graph[i].dist = INF;
        graph[i].closest_node = -1;
    }

    vector<int> P; //facilities
    for (int i = 0; i < nVertices; i++) {
        if (graph[i].isPOI == 1)
            P.push_back(i);
    }

    priority_queue<Vertex> Q;
    vector<bool> visited(nVertices, false);

    for (int u : P) {
        graph[u].dist = 0;
        graph[u].closest_node = u;
        Q.push(graph[u]);
    }

    while (!Q.empty()) {
        Vertex current = Q.top();
        Q.pop();
        int v = current.id;
        long long d = current.dist;
        int u = current.closest_node;

        if (visited[v])
            continue;
        visited[v] = true;

        for (const auto &edge : graph[v].adj) {
            int w = edge.first;
            int dist_v_w = edge.second;

            if (!visited[w] && d + dist_v_w < graph[w].dist) {
                graph[w].dist = d + dist_v_w;
                graph[w].closest_node = u;
                Q.push(graph[w]);
            }
        }
    }
}
void Graph::dijkstra_onetarget(int target) { // to be done
    for (int i = 0; i < nVertices; i++) {
        graph[i].target_dis = INF;
        graph[i].target_edge = intINF;
        graph[i].father_node = -1;
       // graph[i].closest_node = -1;
    }

    vector<int> P;
    P.push_back(target);
    // for (int i = 0; i < nVertices; i++) {
    //     if (graph[i].isPOI == 1 && i!=target)
    //         P.push_back(i);
    // }
    auto cmp=[&](Vertex a,Vertex b){return a.target_dis>b.target_dis;};
    priority_queue<Vertex,vector<Vertex>,decltype(cmp)> Q(cmp);
    vector<bool> visited(nVertices, false);

    for (int u : P) {
        graph[u].target_dis = 0;
        graph[u].target_edge = 0;
        //graph[u].closest_node = u;
        //graph[u].father_node=u;
        Q.push(graph[u]);
    }

    while (!Q.empty()) {
        
        Vertex current = Q.top();
        Q.pop();
        int v = current.id;
        long long d = current.target_dis;
        //int u = current.closest_node;

        if (visited[v])
            continue;
        visited[v] = true;
        for (const auto &edge : graph[v].adj) {
            int w = edge.first;
            if(graph[w].isPOI==1) continue;
            int dist_v_w = edge.second;

            if (!visited[w] && d + dist_v_w < graph[w].target_dis) {
                graph[w].target_dis = d + dist_v_w;
                graph[w].target_edge=graph[v].target_edge+1;
                graph[w].father_node = v;
                //graph[w].closest_node = u;
                Q.push(graph[w]);
            }
        }
    }
}
void Graph::setPOI(vector<int> poi) {
    for (int i : poi)
        graph[i].isPOI = 1;
}

int Graph::updateWeight(int u, int v, int delta) {
    int ans = 0;
    for (auto &edge : graph[u].adj) {
        if (edge.first == v) {
            if (delta == -1)
                ans=edge.second,edge.second = 0;
            else
                edge.second = max(edge.second - delta, 0);
            //ans = edge.second;
            break;
        }
    }
    for (auto &edge : graph[v].adj) {
        if (edge.first == u) {
            if (delta == -1)
                edge.second = 0;
            else
                edge.second = max(edge.second - delta, 0);
            break;
        }
    }
    return ans;
}

void Graph::dynamicDijkstra(int u, int v, int delta) {
    int w = updateWeight(u, v, delta);

    priority_queue<Vertex> Q;
    vector<bool> visited(nVertices, false);

    if (graph[u].dist + w < graph[v].dist) {
        graph[v].dist = graph[u].dist + w;
        graph[v].closest_node = graph[u].closest_node;
        Q.push(graph[v]);
    }

    if (graph[v].dist + w < graph[u].dist) {
        graph[u].dist = graph[v].dist + w;
        graph[u].closest_node = graph[v].closest_node;
        Q.push(graph[u]);
    }

    while (!Q.empty()) {
        Vertex current = Q.top();
        Q.pop();
        int v = current.id;
        long long d = current.dist;
        int u = current.closest_node;

        if (visited[v])
            continue;
        visited[v] = true;

        for (const auto &edge : graph[v].adj) {
            int w = edge.first;
            int dist_v_w = edge.second;

            if (!visited[w] && d + dist_v_w < graph[w].dist) {
                graph[w].dist = d + dist_v_w;
                graph[w].closest_node = u;
                Q.push(graph[w]);
            }
        }
    }
}

void Graph::recomputeDijkstra(int u, int v, int delta) {
    int w = updateWeight(u, v, delta);
    dijkstra();
}

int Graph::exact(int targetID, pair<int,int> budget) {
    int nBefore = 0;
    int gain = -intINF;
    int cost = 0;
    int budgetn=budget.first;
    int budgetw=budget.second;
    vector<tuple<int, int, int>> el;
    for (int i = 0; i < nVertices; ++i) {
        for (auto j : graph[i].adj) {
            if (i < j.first) {
                el.push_back(make_tuple(i, j.first, j.second));
            }
        }
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif

    vector<pair<int, int>> tempAns;
    int m = el.size();
    vector<int> a(m);
    iota(a.begin(), a.end(), 0);
    queue<pair<vector<int>, int>> subsetsQueue;
    subsetsQueue.push({{}, 0});

    while (!subsetsQueue.empty()) {
        auto currentSubset = subsetsQueue.front().first;
        int startIndex = subsetsQueue.front().second;
        subsetsQueue.pop();

        vector<Vertex> gtemp = graph;
#ifdef DEBUG
        cout << "Choose ";
#endif
        int sumw=0;
        if (currentSubset.size() <= budgetn) {
            sumw=0;
            for (int element : currentSubset) {
                int u = get<0>(el[element]);
                int v = get<1>(el[element]);
                int w = get<2>(el[element]);
                sumw+=w;
                dynamicDijkstra(u, v, w);
#ifdef DEBUG
                cout << "(" << u << ", " << v << "), ";
#endif
            }
        }

        int nafter = 0;
        for (int i = 0; i < nVertices; ++i) {
            if (graph[i].closest_node == targetID)
                nafter++;
        }
#ifdef DEBUG
        cout << "After = " << nafter << ", gain = " << nafter - nBefore << endl;
#endif
        int gainAfter = nafter - nBefore;
        if (gain < gainAfter &&sumw<=budgetw) {
            gain = gainAfter;
            tempAns.clear();
            for (int element : currentSubset) {
                int u = get<0>(el[element]);
                int v = get<1>(el[element]);
                tempAns.push_back(make_pair(u, v));
            }
        }
        graph = gtemp;
        cost++;

        if (currentSubset.size() == budgetn) {
            continue;
        }

        for (int i = startIndex; i < m; i++) {
            vector<int> newSubset = currentSubset;
            newSubset.push_back(a[i]);
            subsetsQueue.push({newSubset, i + 1});
        }
    }
#ifdef DEBUG
    cout << "maximum gain is " << gain << endl;
    cout << "select: ";
    for (auto ans : tempAns) {
        cout << "(" << ans.first << ", " << ans.second << "), ";
    }
    cout << endl;
    cout << "cost = " << cost << endl;
#endif
    return gain;
}

void Graph::commputeBound(int targetID) {
    boundList.clear();
    vector<int> a;
    a.clear();
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node != targetID and graph[i].dist != 0)
            a.push_back(graph[i].dist);
    }
    // for(auto e : a) cout << e << endl;

    vector<int> b;
    map<int, int> c;
    map<int, int> d;
    map<int, int> e;

    // Step 1: Generate vector b
    for (int i : a) {
        b.push_back(i - 1);
    }

    // Step 2: Generate map c
    for (int i : b) {
        c[i]++;
    }

    // Step 3: Generate map d (accumulate values of c in reverse order)
    int accumulated = 0;
    for (auto it = c.rbegin(); it != c.rend(); ++it) {
        accumulated += it->second;
        d[it->first] = accumulated;
    }

    // Step 4: Generate map e (fill in gaps between distances)
    int lastKey = d.rbegin()->first;
    int lastValue = d[lastKey];
    for (int i = lastKey; i >= 0; --i) {
        if (d.find(i) != d.end()) {
            e[i] = d[i];
            lastValue = d[i];
        } else {
            e[i] = lastValue;
        }
    }

    boundList.clear();
    // Step 5: Generate vector boundList (only preserve the values in each
    // position)
    for (int i = 0; i < e.size(); ++i) {
        boundList.push_back(e[i]);
    }
}
int Graph::convert_potential(int targetID, pair<int,int> budget, int pruneOP) {//daobile
    int nBefore = 0;
    int nAfter = 0;
    int gain =0;   // each round
    int totalGain = 0; // whole round
    int cost = 0;

    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
    int roundID = 0;
    int budget_n =budget.first;
    int budget_w =budget.second;
    if(budget_w<=0) budget_w=maxedgeweight*budget_n;
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif
    while (budget_n > 0 && budget_w>0) {
    // 1 choose m pot vertex randomly
    int numofpot=0,nownum=0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].isPOI) numofpot++;
    }
    vector<int> pot;
    pot.clear();
    srand(1234^roundID); // Fixed seed for consistency
    for (int i = 0; i < nVertices; i++) {
        int randomNum = rand() % nVertices;
        if(graph[randomNum].isPOI!=1&&graph[randomNum].ispot!=1) graph[randomNum].ispot=1,nownum++,pot.push_back(randomNum);
        if(nownum>=numofpot) break;
    }
    // 2 dij from target to compute each pot vertex 
    // 2.5 dij from all poi and pot to caculate its potiential
    // 3 sort vertex to choose
    // 4 convert the path between
    // 5 recaculate
    for (int i = 0; i < nVertices; i++) {
        graph[i].ispot=0;
    }
    roundID++;
    }
    return totalGain;
}

// int Graph::shortestpathtree_method(int targetID, pair<int,long long> budget, int pruneOP) {
//     int nBefore = 0;
//     int nAfter = 0;
//     int gain =0;   // each round
//     int totalGain = 0; // whole round
//     int cost = 0;

//     for (int i = 0; i < nVertices; ++i) {
//         if (graph[i].closest_node == targetID)
//             nBefore++;
//     }
//     int roundID = 0;
//     int budget_n =budget.first;
//     long long budget_w =budget.second;
//     if(budget_w<=0) budget_w=maxedgeweight*budget_n;
// #ifdef DEBUG
//     cout << "Before = " << nBefore << endl;
// #endif
//     // cout<<"Facility node: ";
//     // for(int i=0;i<nVertices;++i)
//     // {
//     //     if(graph[i].isPOI==1) cout<<i+1<<" ";
//     // }
//     // cout<<"\n";
//     //1 build the spt (root=target node=only user)
//     dijkstra_onetarget(targetID);
//     // for(int i=0;i<nVertices;++i)
//     // {
//     //     cout<<i+1<<" "<<graph[i].dist<<" "<<graph[i].target_dis<<" "<<graph[i].father_node+1<<"\n";
//     // }
    
//     //1.5 calculate the disbound, convertvalue and value for each edge
//     std::vector<std::vector<int>> spt_son(nVertices+1);
//     for(int i=0;i<nVertices;++i){
//         //cout<<i<<" "<<graph[i].father_node<<"\n";
//         if(graph[i].father_node!=-1) spt_son[graph[i].father_node].push_back(i);

//     } 
//     // cout<<"Son: ";
//     // for(int i=0;i<nVertices;++i)
//     // {
//     //     cout<<i<<": ";
//     //     for(auto u:spt_son[i]) cout<<u<<" ";
//     //     cout<<"\n";
//     // }
//     int lim=1;
//     while((1<<lim)<=nVertices) lim++;
//     std::vector<std::vector<int>> fa(nVertices, std::vector<int>(lim+1, -1));
//     std::vector<std::vector<long long>> faweight(nVertices, std::vector<long long>(lim+1, INF));
//     std::vector<int> dis_bound(nVertices+1,0),convert_value(nVertices+1,0),high_fa(nVertices+1,0);
//     for(int i=0;i<nVertices;++i) dis_bound[i]=graph[i].dist;
//     for(int i=0;i<nVertices;++i) fa[i][0]=i,faweight[i][0]=0;
//     //for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?edge[{i,fa[i][1]}]:INF; //test shujuji you chongbian!
//     for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?graph[i].target_dis-graph[fa[i][1]].target_dis:INF;
//     for(int j=2;j<=lim;++j)
//     {
//         for(int i=0;i<nVertices;++i)
//         {
//             if(fa[i][j-1]!=-1)
//             {
//                 fa[i][j]=fa[fa[i][j-1]][j-1];
//                 if(fa[i][j]!=-1) faweight[i][j]=faweight[i][j-1]+faweight[fa[i][j-1]][j-1];
//                 else faweight[i][j]=INF;
//             }
//             else
//             {
//                 fa[i][j]=-1;faweight[i][j]=INF;
//             }
//         }

//     }
//     // for(int i=0;i<nVertices;++i)
//     // {
//     //     for(int j=1;j<=lim;++j) cout<<i<<"'s "<<(1<<(j-1))<<"th ansector :"<<fa[i][j]<<" "<<faweight[i][j]<<"\n";
//     // }
//     for(int i=0;i<nVertices;++i)
//     {
//         int tmp=dis_bound[i];
//         if(graph[i].isPOI==1)continue;
//         if(graph[i].closest_node==targetID) continue;
        
//     }
//     for(int i=0;i<nVertices;i++)
//     {
//         int tmpd=dis_bound[i],ti=i,tmw=0;
//         while(ti!=-1)
//         {
//             assert(tmpd>=0);
//             int j=0;
//             while(j<=lim&&faweight[ti][j+1]<=tmpd) j++;
//             tmpd-=faweight[ti][j];
//             tmw+=faweight[ti][j];
//             ti=fa[ti][j];
//             if(j==0)
//             {
//                 //if(tmw!=graph[i].target_dis-graph[ti].target_dis) cout<<i<<"\n";
//                 assert(tmw==graph[i].target_dis-graph[ti].target_dis);
//                 high_fa[i]=ti;break;
//             }
//         }
//     }
//     for(int i=0;i<nVertices;i++)
//     {
//         if(graph[i].isPOI==1)continue;
//         if(graph[i].closest_node==targetID) continue;
//         //cout<<i<<" "<<dis_bound[i]<<" "<<high_fa[i]<<"\n";
//         // if(graph[i].target_dis-graph[high_fa[i]].target_dis>dis_bound[i])
//         // {
//         //     cout<<i<<" "<<high_fa[i]<<" "<<graph[i].target_dis-graph[high_fa[i]].target_dis<<" "<<dis_bound[i]<<"\n";
//         // }
//         assert(graph[i].target_dis-graph[high_fa[i]].target_dis<=dis_bound[i]);
//         assert(fa[high_fa[i]][1]==-1||graph[i].target_dis-graph[fa[high_fa[i]][1]].target_dis>dis_bound[i]);
//         convert_value[high_fa[i]]++;
//     }
//     auto dfs = [&](auto&& self, int x) -> void {
//         for(auto u:spt_son[x])
//         {
//             convert_value[u]+=convert_value[x];
//             self(self,u);
//         }
//         return;
//     };
//     dfs(dfs,targetID);
//     // for(int i=0;i<nVertices;++i)
//     // {
//     //     if(fa[i][1]==-1) continue;
//     //     cout<<i<<" "<<convert_value[i]-convert_value[fa[i][1]]<<" "<<graph[i].target_edge<<"\n";
//     // };


//     //2 convert neighbor edge with max v/w
//     //std::vector<std::tuple<int, int, int>> a(nVertices);
//     auto cmp =[&](std::tuple<int, int, int>& t1, std::tuple<int, int, int>& t2){
//         if(pruneOP==1) return 1ll*std::get<1>(t1) * 1ll*std::get<2>(t2) < 1ll * std::get<1>(t2)*1ll * std::get<2>(t1);
//         return 1ll*std::get<1>(t1) < 1ll * std::get<1>(t2);
//     };
//     std::priority_queue<std::tuple<int, int, int>, std::vector<std::tuple<int, int, int>>, decltype(cmp)> pq(cmp);
//     for(auto u:spt_son[targetID]) pq.push(make_tuple(u, convert_value[u] - convert_value[graph[u].father_node], graph[u].target_dis - graph[graph[u].father_node].target_dis));
//     int nown=0,noww=0,exp_gain=0;
//     vector<int> ansnode;
//     while (!pq.empty() && budget_n > nown && budget_w > noww) {
//         auto t=pq.top();
//         pq.pop();
//         //cout<<"pq :"<<get<0>(t)<<" "<<get<1>(t)<<" "<<get<2>(t)<<" \n";
//         if(get<2>(t)+noww <= budget_w)
//         {
//             exp_gain+=get<1>(t);
//             ansnode.push_back(get<0>(t));
//             nown++;
//             noww+=get<2>(t);
//             for(auto u:spt_son[get<0>(t)]) pq.push(make_tuple(u, convert_value[u] - convert_value[graph[u].father_node], graph[u].target_dis - graph[graph[u].father_node].target_dis));
//         }
//     }
//     //cout<<budget_n<<" "<<budget_w<<" \n";
//     //cout<<nown<<" "<<noww<<" "<<exp_gain<< "\n";
//     //cout<<"!!:";
//     for(auto u:ansnode)
//     {
//         //cout<<u<<" ";
//         updateWeight(u,graph[u].father_node,-1);
//     }
//     //cout<<"\n";
//     dijkstra();
//     int nNow=0;
//     for (int i = 0; i < nVertices; ++i) {
//         if (graph[i].closest_node == targetID)
//             nNow++;
//     }
//     totalGain=nNow-nBefore;
//     // if(totalGain<exp_gain) // to be done
//     // {
//     //     cout<<"?? "<<totalGain<<" "<<exp_gain<<"\n";
//     // }
//     // assert(totalGain>=exp_gain);
//     return totalGain;
// }


// calculate lca to calculate convert value
int Graph::shortestpathtree_methodgreedy(int targetID, pair<int,long long> budget, int pruneOP) { //slow
    int nBefore = 0;
    int nAfter = 0;
    int gain =0;   // each round
    int totalGain = 0; // whole round
    int cost = 0;
    
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
    int roundID = 0;
    int budget_n =budget.first;
    long long budget_w =budget.second;
    if(budget_w<=0) budget_w=maxedgeweight*budget_n;
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif
    int nCurrent=nBefore;
    long long nown=0,noww=0;
    auto find_best_in_tree=[&]()
    {
        
        dijkstra_onetarget(targetID);
        std::vector<std::vector<int>> spt_son(nVertices+1);
        for(int i=0;i<nVertices;++i){
            if(graph[i].father_node!=-1) spt_son[graph[i].father_node].push_back(i);
        } 

        int lim=1;
        while((1<<lim)<=nVertices) lim++;
        std::vector<std::vector<int>> fa(nVertices, std::vector<int>(lim+1, -1));
        std::vector<std::vector<long long>> faweight(nVertices, std::vector<long long>(lim+1, INF));
        std::vector<long long> dis_bound(nVertices+1,0),convert_value(nVertices+1,0),high_fa(nVertices+1,0);
        for(int i=0;i<nVertices;++i) dis_bound[i]=graph[i].dist;
        for(int i=0;i<nVertices;++i) fa[i][0]=i,faweight[i][0]=0;
        //for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?edge[{i,fa[i][1]}]:INF; //test shujuji you chongbian!
        for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?graph[i].target_dis-graph[fa[i][1]].target_dis:INF;
        for(int i=0;i<nVertices;++i)
        {
            long long tmp=dis_bound[i];
            if(graph[i].isPOI==1)continue;
            if(graph[i].closest_node==targetID) continue;
            
        }
        int bestid=-1,bestvalue=0,bestedge=0;
        vector<Vertex> gtemp = graph;
        for(int i=0;i<nVertices;i++)
        {

            if(i==targetID)continue;
            if(nown+graph[i].target_edge<=budget_n&&noww+graph[i].target_dis<=budget_w){
                //cout<<i<<" "<<graph[i].target_edge<<" "<<graph[i].target_dis<<endl;

                int tmpid=i;
                while (tmpid!=-1&&tmpid!=targetID)
                {
                    int tm=updateWeight(tmpid,graph[tmpid].father_node,-1);
                    //cout<<tm<<"?\n";
                    //if(tm) nown++;
                    tmpid=graph[tmpid].father_node;
                }
                dijkstra();
                int nnow=0;
                for (int i = 0; i < nVertices; ++i) {
                    if (graph[i].closest_node == targetID)
                    nnow++;
                }
                //if((nnow-nCurrent)*bestedge>bestvalue*(graph[i].target_edge))
                if(nnow-nCurrent>bestvalue)
                {
                    bestid=i; bestvalue=nnow-nCurrent;bestedge=graph[i].target_edge;
                }
                graph=gtemp;
                //bestid=i; bestvalue=convert_value[i];bestedge=graph[i].target_edge;
            }
            

        }
        int tmpgain=0,now=nCurrent,tmpid=bestid;
        //cout<<"????:"<<bestid<<endl;
        if(tmpid!=-1)
        {
            //nown=nown+graph[tmpid].target_edge;
            noww=noww+graph[tmpid].target_dis;
            while (tmpid!=-1&&tmpid!=targetID)
            {
                int tm=updateWeight(tmpid,graph[tmpid].father_node,-1);
                //cout<<tm<<"?\n";
                if(tm) nown++;
                tmpid=graph[tmpid].father_node;
            }
            
        }
        //cout<<"????:"<<bestid<<endl;
        
        dijkstra();
        int nNow=0;
        for (int i = 0; i < nVertices; ++i) {
            if (graph[i].closest_node == targetID)
            nNow++;
        }
        //if(bestid!=-1) assert(nNow-nCurrent>=convert_value[bestid]); //to be done
        if(nNow==nCurrent) return -1;
        nCurrent=nNow;
        //cout<<"??:"<<bestid<<endl;
        return bestid;

    };
    while(nown<=budget_n&&noww<=budget_w)
    {
        int tm=find_best_in_tree();
        if(tm==-1) break;
    }
    //cout<<nown<<" "<<noww<<"\n";
    dijkstra();
    int nNow=0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nNow++;
    }
    totalGain=nNow-nBefore;

    return totalGain;
}

int Graph::shortestpathtree_method0(int targetID, pair<int,long long> budget, int pruneOP) { //slow
    int nBefore = 0;
    int nAfter = 0;
    int gain =0;   // each round
    int totalGain = 0; // whole round
    int cost = 0;

    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
    int roundID = 0;
    int budget_n =budget.first;
    long long budget_w =budget.second;
    if(budget_w<=0) budget_w=maxedgeweight*budget_n;
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif
    int nCurrent=nBefore;
    long long nown=0,noww=0;
    auto find_best_in_tree=[&]()
    {
        dijkstra_onetarget(targetID);
        std::vector<std::vector<int>> spt_son(nVertices+1);
        for(int i=0;i<nVertices;++i){
            if(graph[i].father_node!=-1) spt_son[graph[i].father_node].push_back(i);
        } 

        int lim=1;
        while((1<<lim)<=nVertices) lim++;
        std::vector<std::vector<int>> fa(nVertices, std::vector<int>(lim+1, -1));
        std::vector<std::vector<long long>> faweight(nVertices, std::vector<long long>(lim+1, INF));
        std::vector<long long> dis_bound(nVertices+1,0),convert_value(nVertices+1,0),high_fa(nVertices+1,0);
        for(int i=0;i<nVertices;++i) dis_bound[i]=graph[i].dist;
        for(int i=0;i<nVertices;++i) fa[i][0]=i,faweight[i][0]=0;
        //for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?edge[{i,fa[i][1]}]:INF; //test shujuji you chongbian!
        for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?graph[i].target_dis-graph[fa[i][1]].target_dis:INF;
        // for(int j=2;j<=lim;++j)
        // {
        //     for(int i=0;i<nVertices;++i)
        //     {
        //         if(fa[i][j-1]!=-1)
        //         {
        //             fa[i][j]=fa[fa[i][j-1]][j-1];
        //             if(fa[i][j]!=-1) faweight[i][j]=faweight[i][j-1]+faweight[fa[i][j-1]][j-1];
        //             else faweight[i][j]=INF;
        //         }
        //         else
        //         {
        //             fa[i][j]=-1;faweight[i][j]=INF;
        //         }
        //     }

        // }
        for(int i=0;i<nVertices;++i)
        {
            long long tmp=dis_bound[i];
            if(graph[i].isPOI==1)continue;
            if(graph[i].closest_node==targetID) continue;
            
        }
        for(int i=0;i<nVertices;i++)
        {
            long long tmpd=dis_bound[i],ti=i,tmw=0,flag=0;
            if(graph[i].isPOI==1&&i!=targetID)continue;
            long long tmcv=0;
            //if(graph[i].closest_node==targetID) continue;
            auto dfs1 = [&](auto&& self, int x) -> void {
                for(auto u:spt_son[x])
                {
                    //convert_value[u]+=convert_value[x];
                    if(graph[u].closest_node!=targetID)
                        if(graph[u].target_dis-graph[i].target_dis<=dis_bound[u]&&(graph[i].father_node==-1||graph[u].target_dis-graph[graph[i].father_node].target_dis>dis_bound[u])) tmcv++;
                    self(self,u);
                }
                return;
            };
            dfs1(dfs1,i);
            convert_value[i]=tmcv;
        }

        // for(int i=0;i<nVertices;i++)
        // {
        //     if(graph[i].isPOI==1)continue;
        //     if(graph[i].closest_node==targetID) continue;

        //     assert(graph[i].target_dis-graph[high_fa[i]].target_dis<=dis_bound[i]);
        //     assert(fa[high_fa[i]][1]==-1||graph[i].target_dis-graph[fa[high_fa[i]][1]].target_dis>dis_bound[i]);
        //     convert_value[high_fa[i]]++;
        // }

        auto dfs = [&](auto&& self, int x) -> void {
            for(auto u:spt_son[x])
            {
                convert_value[u]+=convert_value[x];
                self(self,u);
            }
            return;
        };
        dfs(dfs,targetID);
        int bestid=-1,bestvalue=0,bestedge=0;
        for(int i=0;i<nVertices;i++)
        {

            if(graph[i].isPOI==1)continue;
            if(nown+graph[i].target_edge<=budget_n&&noww+graph[i].target_dis<=budget_w&&(convert_value[i]>=bestvalue||convert_value[i]==bestvalue&&graph[i].target_edge<=bestedge)){
                bestid=i; bestvalue=convert_value[i];bestedge=graph[i].target_edge;
            }
        }
        int tmpgain=0,now=nCurrent,tmpid=bestid;
        if(tmpid!=-1)
        {
            //nown=nown+graph[tmpid].target_edge;
            noww=noww+graph[tmpid].target_dis;
            while (tmpid!=-1&&tmpid!=targetID)
            {
                int tm=updateWeight(tmpid,graph[tmpid].father_node,-1);
                //cout<<tm<<"?\n";
                if(tm) nown++;
                tmpid=graph[tmpid].father_node;
            }

        }
        //cout<<"????:"<<bestid<<endl;
        
        dijkstra();
        int nNow=0;
        for (int i = 0; i < nVertices; ++i) {
            if (graph[i].closest_node == targetID)
            nNow++;
        }
        //if(bestid!=-1) assert(nNow-nCurrent>=convert_value[bestid]); //to be done
        if(nNow==nCurrent) return -1;
        nCurrent=nNow;
        //cout<<"??:"<<bestid<<endl;
        return bestid;

    };
    while(nown<=budget_n&&noww<=budget_w)
    {
        int tm=find_best_in_tree();
        if(tm==-1) break;
    }
    //cout<<nown<<" "<<noww<<"\n";
    dijkstra();
    int nNow=0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nNow++;
    }
    totalGain=nNow-nBefore;

    return totalGain;
}


//current best version
int Graph::shortestpathtree_method(int targetID, pair<int,long long> budget, int pruneOP) { //slow
    int nBefore = 0;
    int nAfter = 0;
    int gain =0;   // each round
    int totalGain = 0; // whole round
    int cost = 0;

    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
    int roundID = 0;
    int budget_n =budget.first;
    long long budget_w =budget.second;
    if(budget_w<=0) budget_w=maxedgeweight*budget_n;
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif
    int nCurrent=nBefore;
    long long nown=0,noww=0;
    auto find_best_in_tree=[&]()
    {
        dijkstra_onetarget(targetID);
        std::vector<std::vector<int>> spt_son(nVertices+1);
        for(int i=0;i<nVertices;++i){
            if(graph[i].father_node!=-1) spt_son[graph[i].father_node].push_back(i);
        } 

        int lim=1;
        while((1<<lim)<=nVertices) lim++;
        std::vector<std::vector<int>> fa(nVertices, std::vector<int>(lim+1, -1));
        std::vector<std::vector<long long>> faweight(nVertices, std::vector<long long>(lim+1, INF));
        std::vector<long long> dis_bound(nVertices+1,0),convert_value(nVertices+1,0),high_fa(nVertices+1,0);
        for(int i=0;i<nVertices;++i) dis_bound[i]=graph[i].dist;
        for(int i=0;i<nVertices;++i) fa[i][0]=i,faweight[i][0]=0;
        //for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?edge[{i,fa[i][1]}]:INF; //test shujuji you chongbian!
        for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?graph[i].target_dis-graph[fa[i][1]].target_dis:INF;

        
        for(int i=0;i<nVertices;++i)
        {
            long long tmp=dis_bound[i];
            if(graph[i].isPOI==1)continue;
            if(graph[i].closest_node==targetID) continue;
            
        }
        for(int i=0;i<nVertices;i++)
        {
            long long tmpd=dis_bound[i],ti=i,tmw=0,flag=0;
            if(graph[i].isPOI==1)continue;
            if(graph[i].closest_node==targetID) continue;
            while(ti!=-1)
            {
                assert(tmpd>=0);
                int j=0;
                while(tmw+faweight[ti][1]<=tmpd) tmw+=faweight[ti][1],ti=fa[ti][1];
                //if(tmw!=graph[i].target_dis-graph[ti].target_dis) cout<<i<<"\n";
                assert(tmw==graph[i].target_dis-graph[ti].target_dis);
                high_fa[i]=ti;break;
            }
 
        }
        for(int i=0;i<nVertices;i++)
        {
            if(graph[i].isPOI==1)continue;
            if(graph[i].closest_node==targetID) continue;
            assert(graph[i].target_dis-graph[high_fa[i]].target_dis<=dis_bound[i]);
            assert(fa[high_fa[i]][1]==-1||graph[i].target_dis-graph[fa[high_fa[i]][1]].target_dis>dis_bound[i]);
            convert_value[high_fa[i]]++;
        }
        auto dfs = [&](auto&& self, int x) -> void {
            for(auto u:spt_son[x])
            {
                convert_value[u]+=convert_value[x];
                self(self,u);
            }
            return;
        };
        dfs(dfs,targetID);
        int bestid=-1,bestvalue=0,bestedge=0;
        for(int i=0;i<nVertices;i++)
        {

            if(graph[i].isPOI==1)continue;
            if(nown+graph[i].target_edge<=budget_n&&noww+graph[i].target_dis<=budget_w&&(convert_value[i]>=bestvalue||convert_value[i]==bestvalue&&graph[i].target_edge<=bestedge)){
                bestid=i; bestvalue=convert_value[i];bestedge=graph[i].target_edge;
            }
            
            // if(nown+graph[i].target_edge<=budget_n&&noww+graph[i].target_dis<=budget_w&&((long long)convert_value[i]*bestedge>=(long long)bestvalue*graph[i].target_edge)){
            //     bestid=i; bestvalue=convert_value[i];bestedge=graph[i].target_edge;
            // }
        }
        int tmpgain=0,now=nCurrent,tmpid=bestid;
        if(tmpid!=-1)
        {
            //nown=nown+graph[tmpid].target_edge;
            noww=noww+graph[tmpid].target_dis;
            while (tmpid!=-1&&tmpid!=targetID)
            {
                int tm=updateWeight(tmpid,graph[tmpid].father_node,-1);
                //cout<<tm<<"?\n";
                if(tm) nown++;
                tmpid=graph[tmpid].father_node;
            }

        }

        
        dijkstra();
        int nNow=0;
        for (int i = 0; i < nVertices; ++i) {
            if (graph[i].closest_node == targetID)
            nNow++;
        }
        //if(bestid!=-1) assert(nNow-nCurrent>=convert_value[bestid]); //to be done
        if(nNow==nCurrent) return -1;
        nCurrent=nNow;
        //cout<<"??:"<<bestid<<endl;
        return bestid;

    };
    while(nown<=budget_n&&noww<=budget_w)
    {
        int tm=find_best_in_tree();
        if(tm==-1) break;
    }
    //cout<<nown<<" "<<noww<<"\n";
    dijkstra();
    int nNow=0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nNow++;
    }
    totalGain=nNow-nBefore;

    return totalGain;
}



// discard version
int Graph::shortestpathtree_method2(int targetID, pair<int,long long> budget, int pruneOP) {  //using binary lifting to speed up
    int nBefore = 0;
    int nAfter = 0;
    int gain =0;   // each round
    int totalGain = 0; // whole round
    int cost = 0;

    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
    int roundID = 0;
    int budget_n =budget.first;
    long long budget_w =budget.second;
    if(budget_w<=0) budget_w=maxedgeweight*budget_n;
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif
    int nCurrent=nBefore;
    long long nown=0,noww=0;
    auto find_best_in_tree=[&]()
    {
        dijkstra_onetarget(targetID);
        std::vector<std::vector<int>> spt_son(nVertices+1);
        for(int i=0;i<nVertices;++i){
            if(graph[i].father_node!=-1) spt_son[graph[i].father_node].push_back(i);
        } 

        int lim=1;
        while((1<<lim)<=nVertices) lim++;
        std::vector<std::vector<int>> fa(nVertices, std::vector<int>(lim+1, -1));
        std::vector<std::vector<long long>> faweight(nVertices, std::vector<long long>(lim+1, INF));
        std::vector<long long> dis_bound(nVertices+1,0),convert_value(nVertices+1,0),high_fa(nVertices+1,0);
        for(int i=0;i<nVertices;++i) dis_bound[i]=graph[i].dist;
        for(int i=0;i<nVertices;++i) fa[i][0]=i,faweight[i][0]=0;
        //for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?edge[{i,fa[i][1]}]:INF; //test shujuji you chongbian!
        for(int i=0;i<nVertices;++i) fa[i][1]=graph[i].father_node,faweight[i][1]=(fa[i][1]!=-1)?graph[i].target_dis-graph[fa[i][1]].target_dis:INF;
        for(int j=2;j<=lim;++j)
        {
            for(int i=0;i<nVertices;++i)
            {
                if(fa[i][j-1]!=-1)
                {
                    fa[i][j]=fa[fa[i][j-1]][j-1];
                    if(fa[i][j]!=-1) faweight[i][j]=faweight[i][j-1]+faweight[fa[i][j-1]][j-1];
                    else faweight[i][j]=INF;
                }
                else
                {
                    fa[i][j]=-1;faweight[i][j]=INF;
                }
            }

        }
        // for(int i=0;i<nVertices;++i)
        // {
        //     long long tmp=dis_bound[i];
        //     if(graph[i].isPOI==1)continue;
        //     if(graph[i].closest_node==targetID) continue;
            
        // }
        for(int i=0;i<nVertices;i++)
        {
            if(graph[i].isPOI==1)continue;
            if(graph[i].closest_node==targetID) continue;
            long long tmpd=dis_bound[i],ti=i,tmw=0,nowlim=lim;
            while(1)
            {
                assert(tmpd>=0);
                int j=nowlim;
                while(j>0 && faweight[ti][j]>tmpd) j--;
                tmpd-=faweight[ti][j];
                tmw+=faweight[ti][j];
                ti=fa[ti][j];
                nowlim=j-1;
                if(j==0)
                {
                    //if(tmw!=graph[i].target_dis-graph[ti].target_dis) cout<<i<<"\n";
                    assert(tmw==graph[i].target_dis-graph[ti].target_dis);
                    high_fa[i]=ti;break;
                }
                
            }
            // while(ti!=-1)
            // {
            //     assert(tmpd>=0);
            //     int j=0;
            //     while(j<=lim && faweight[ti][j+1]<=tmpd) j++;
            //     tmpd-=faweight[ti][j];
            //     tmw+=faweight[ti][j];
            //     ti=fa[ti][j];
            //     if(j==0)
            //     {
            //         //if(tmw!=graph[i].target_dis-graph[ti].target_dis) cout<<i<<"\n";
            //         assert(tmw==graph[i].target_dis-graph[ti].target_dis);
            //         high_fa[i]=ti;break;
            //     }
                
            // }
        }
        for(int i=0;i<nVertices;i++)
        {
            if(graph[i].isPOI==1)continue;
            if(graph[i].closest_node==targetID) continue;

            assert(graph[i].target_dis-graph[high_fa[i]].target_dis<=dis_bound[i]);
            assert(fa[high_fa[i]][1]==-1||graph[i].target_dis-graph[fa[high_fa[i]][1]].target_dis>dis_bound[i]);
            convert_value[high_fa[i]]++;
        }
        auto dfs = [&](auto&& self, int x) -> void {
            for(auto u:spt_son[x])
            {
                convert_value[u]+=convert_value[x];
                self(self,u);
            }
            return;
        };
        dfs(dfs,targetID);
        int bestid=-1,bestvalue=0,bestedge=0;
        for(int i=0;i<nVertices;i++)
        {

            if(graph[i].isPOI==1)continue;
            if(nown+graph[i].target_edge<=budget_n&&noww+graph[i].target_dis<=budget_w&&(convert_value[i]>=bestvalue||convert_value[i]==bestvalue&&graph[i].target_edge<=bestedge)){
                bestid=i; bestvalue=convert_value[i];bestedge=graph[i].target_edge;
            }
        }
        int tmpgain=0,now=nCurrent,tmpid=bestid;
        if(tmpid!=-1)
        {
            //nown=nown+graph[tmpid].target_edge;
            noww=noww+graph[tmpid].target_dis;
            while (tmpid!=-1&&tmpid!=targetID)
            {
                int tm=updateWeight(tmpid,graph[tmpid].father_node,-1);
                //cout<<tm<<"?\n";
                if(tm) nown++;
                tmpid=graph[tmpid].father_node;
            }

        }
        //cout<<"????:"<<bestid<<endl;
        
        dijkstra();
        int nNow=0;
        for (int i = 0; i < nVertices; ++i) {
            if (graph[i].closest_node == targetID)
            nNow++;
        }
        //if(bestid!=-1) assert(nNow-nCurrent>=convert_value[bestid]); //to be done
        if(nNow==nCurrent) return -1;
        nCurrent=nNow;
        //cout<<"??:"<<bestid<<endl;
        return bestid;

    };
    while(nown<=budget_n&&noww<=budget_w)
    {
        int tm=find_best_in_tree();
        if(tm==-1) break;
    }
    //cout<<nown<<" "<<noww<<"\n";
    dijkstra();
    int nNow=0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nNow++;
    }
    totalGain=nNow-nBefore;

    return totalGain;
}



int Graph::greedy(int targetID, pair<int,long long> budget, int pruneOP) {
    int nBefore = 0;
    int nAfter = 0;
    int gain = -intINF;   // each round
    int totalGain = 0; // whole round
    int cost = 0;

    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif
    commputeBound(targetID);

    int roundID = 0;
    int budget_n =budget.first;
    long long budget_w =budget.second;
    if(budget_w<=0) budget_w=maxedgeweight*budget_n;
    //cout<<budget_w<<endl;
    while (budget_n > 0 && budget_w>0) {
        budget_n--;
        tuple<int, int, int> tempAnsEpoh(0, 0, 0);
        vector<Vertex> gtemp = graph;
        vector<int> dist(nVertices, numeric_limits<int>::max());
        vector<bool> visited(nVertices, false);
        vector<int> max_edge_on_path(nVertices, 0);

        priority_queue<pii, vector<pii>, greater<pii>> pq;
        pq.push({0, targetID});
        dist[targetID] = 0;
        set<pair<int, int>> uvPairs;
        uvPairs.clear();

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            if (visited[u])
                continue;
            visited[u] = true;

            if (pruneOP & Prune::Prune3) {
                if (boundList[dist[u]] <= gain)
                    break;
            }

            for (const auto &edge : graph[u].adj) {
                int v = edge.first;
                int w = edge.second;
                if(w>budget_w) continue;
                auto [it, inserted] = uvPairs.insert(make_pair(u, v));
                auto [it2, inserted2] = uvPairs.insert(make_pair(v, u));
                if (inserted == false or inserted2 == false)
                    continue;

                if (pruneOP & Prune::Prune1) {

                    // prune if both u and v not contain a targetID
                    if (graph[u].closest_node != targetID and
                        graph[v].closest_node != targetID)
                        continue;
                    // prune if v not contains a targetID and is far
                    if (graph[v].closest_node != targetID and
                        graph[v].dist < graph[u].dist)
                        continue;
                    if (graph[u].closest_node != targetID and
                        graph[u].dist < graph[v].dist)
                        continue;
                }

                if (pruneOP & Prune::Prune2) {
                    // prune if the weight is smaller
                    if (w <= max_edge_on_path[u])
                        continue;
                }

                cost++;
                dynamicDijkstra(u, v, w);
#ifdef DEBUG
                cout << "(" << u << ", " << v << "), ";
#endif
                int cnt = 0;
                for (int j = 0; j < nVertices; ++j) {
                    if (graph[j].closest_node == targetID)
                        cnt++;
                }
                //if(cnt>10000)  cout << "bao:(" << u << ", " << v << "), target="<<targetID<<" "<<cost<<endl; //you turan tebieduo de qingkuang
                int gainAfter = cnt - nBefore; // gap
                // cout << "After = " << cnt << ", Before = " << nBefore
                //      << ", gain = " << gainAfter << endl;
                
#ifdef DEBUG
                cout << "After = " << cnt << ", Before = " << nBefore
                     << ", gain = " << gainAfter << endl;
#endif
                if (gainAfter > gain) {
                    gain = gainAfter;
                    tempAnsEpoh = make_tuple(u, v, w);
                    nAfter = cnt;
                }
                graph = gtemp;
            }

            for (const auto &edge : graph[u].adj) {
                int v = edge.first;
                int w = edge.second;

                if (dist[u] + w < dist[v]) {
                    dist[v] = dist[u] + w;
                    max_edge_on_path[v] = max(max_edge_on_path[u], w);
                    pq.push({dist[v], v});
                } else if (dist[u] + w == dist[v]) {
                    max_edge_on_path[v] =
                        max(max_edge_on_path[v], max(max_edge_on_path[u], w));
                }
            }
        }

        int u = get<0>(tempAnsEpoh);
        int v = get<1>(tempAnsEpoh);
        int w = get<2>(tempAnsEpoh);
        dynamicDijkstra(u, v, w);
        nBefore = nAfter;
        budget_w -= w;
        totalGain += gain;
#ifdef DEBUG
        cout << "round " << ++roundID << " choose (" << u << " , " << v << ")"
             << endl;
#endif
        //gain = -INF;
        gain = 0;
        commputeBound(targetID);
    }
#ifdef DEBUG
    cout << "The final gain = " << totalGain << endl;
    cout << "cost = " << cost << endl;
#endif

    //cout << "cost = " << cost << endl;

    if(totalGain <= 0) totalGain = 0;
    return totalGain;
}

struct EdgeWeightComparator {
    bool operator()(const tuple<int, int, int> &e1,
                    const tuple<int, int, int> &e2) const {
        return get<2>(e1) > get<2>(e2);
    }
};

int Graph::maxWeight(int targetID, pair<int,int> budget) {
    int nBefore = 0;
    int noww = 0;
    int budgetn = budget.first;
    int budgetw = budget.second;
    vector<tuple<int, int, int>> el;
    for (int i = 0; i < nVertices; ++i) {
        for (auto j : graph[i].adj) {
            if (i < j.first) {
                el.push_back(make_tuple(i, j.first, j.second));
            }
        }
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif

    priority_queue<tuple<int, int, int>, vector<tuple<int, int, int>>,
                   EdgeWeightComparator>
        min_heap;

    for (const auto &edge : el) {
        if (min_heap.size() < budgetn) {
            min_heap.push(edge);
        } else if (get<2>(edge) > get<2>(min_heap.top())) {
            min_heap.pop();
            min_heap.push(edge);
        }
    }

    vector<tuple<int, int, int>> largest_edges;
    while (!min_heap.empty()) {
        tuple<int, int, int> largest_edge = min_heap.top();
        min_heap.pop();

        int u = get<0>(largest_edge);
        int v = get<1>(largest_edge);
        int w = get<2>(largest_edge);
        if (noww + w > budgetw) continue; // skip if weight exceeds budget
        noww += w;
        dynamicDijkstra(u, v, w);
#ifdef DEBUG
        cout << "(" << u << ", " << v << "), ";
#endif
    }
#ifdef DEBUG
    cout << endl;
#endif

    int nafter = 0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nafter++;
    }
#ifdef DEBUG
    cout << "After = " << nafter << ", gain = " << nafter - nBefore << endl;
#endif
    int gain = nafter - nBefore;
#ifdef DEBUG
    cout << "maximum gain is " << gain << endl;
#endif

    return gain;
}

int Graph::chooseNbr(int targetID, pair<int,int> budget) {
    int nBefore = 0;
    int noww=0;
    int budgetn = budget.first;
    int budgetw = budget.second;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif

    vector<bool> visited(nVertices, false);
    queue<int> q;
    set<tuple<int, int, int>> el;
    visited[targetID] = true;
    q.push(targetID);

    while (!q.empty() && el.size() < budgetn) {
        int u = q.front();
        q.pop();

        for (const auto &edge : graph[u].adj) {
            int v = edge.first;
            int w = edge.second;
            if(noww+w > budgetw) continue; // skip if weight exceeds budget
            if (!visited[v]) {
                visited[v] = true;
                q.push(v);
            }

            tuple<int, int, int> e1 = make_tuple(u, v, w);
            tuple<int, int, int> e2 = make_tuple(v, u, w);
            if (el.find(e1) == el.end() && el.find(e2) == el.end()) {
                el.insert(e1);
                if (el.size() >= budgetn) {
                    break;
                }
            }
            noww += w;
        }
    }

    for (const auto &edge : el) {
        int u = get<0>(edge);
        int v = get<1>(edge);
        int w = get<2>(edge);
        dynamicDijkstra(u, v, w);
#ifdef DEBUG
        cout << "(" << u << ", " << v << "), ";
#endif
    }
#ifdef DEBUG
    cout << endl;
#endif

    int nafter = 0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nafter++;
    }
#ifdef DEBUG
    cout << "After = " << nafter << ", gain = " << nafter - nBefore << endl;
#endif
    int gain = nafter - nBefore;
#ifdef DEBUG
    cout << "maximum gain is " << gain << endl;
#endif
    return gain;
}
