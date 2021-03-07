//
//  Graph.hpp
//  pscan- index
//
//  Created by 孟令凯 on 2021/2/23.
//

#ifndef Graph_hpp
#define Graph_hpp

#include <iostream>
#include <stdio.h>
#include "string"
#include <ext/hash_map>
#include <vector>
#include <list>
#include <fstream>
#include <queue>
#include <set>

using namespace std;
using namespace __gnu_cxx;


class Graph {
private:
    string dir; //input graph directory
    
    int n, in_m, out_m,m,degree_m; //number of nodes and edges of the graph

    int eps_a1, eps_b1,eps_a2, eps_b2, miu; // eps_a2/eps_b2 = eps^2
    
    vector<int> pstart,in_pstart,out_pstart; //offset of neighbors of nodes 节点相邻点偏移量
    vector<int> edges,in_edges1,out_edges1; //adjacent ids of edges 相邻边
    vector<int> reverse; //the position of reverse edge in edges 反向边在边中的位置
    vector<int> min_cn1, min_cn2; //minimum common neighbor: -2 means not similar; -1 means similar; 0 means not sure; > 0 means the minimum common neighbor
    vector<int> degree,in_degree,out_degree;
    
    vector<int> similar_degree, effective_degree;
    
//    vector<int> cNum, fNum;
    
//    hash_map<int,hash_map<int,int>> cNum;
//    hash_map<int,hash_map<int,int>> fNum;
    
    vector<int> pa,rank;
    
    vector<int> cid;
    
    hash_map<int,vector<int>> out_edges;
    hash_map<int,vector<int>> in_edges;
    
    vector<pair<int,int>> noncore_cluster;

public:
    Graph(const string _dir) ;
    void read_graph() ;
    void pSCAN(const char *eps_s1,const char *eps_s2, int miu) ;
    void getNeiNum();
    void get_eps(const char *eps_s1,const char *eps_s2);
    void prune_and_cross_link(int eps_a1, int eps_b1,int eps_a2, int eps_b2, int miu, int *cores, int &cores_e) ;
//    int compute_common_neighbor_lowerbound(int u,int v,int eps_a2,int eps_b2, int k);
    int binary_search(vector<int> &array, int b, int e, int val);
    int find_root(int u);
    int check_common_neighbor(int u, int v, int c1,int c2);

    int similar_check_OP(int u, int idx, int eps_a1, int eps_b1, int eps_a2, int eps_b2);
    vector<int> get_f_and_c_num(int a,int b, int c);
    int BinarySearch(vector<int> &a, int value, int n);
    void my_union(int u,int v) ;
    
    void cluster_noncore_vertices(int eps_a1, int eps_b1,int eps_a2, int eps_b2, int mu) ;
    
    void output();

    

private:
    
    
};


#endif /* Graph_hpp */

