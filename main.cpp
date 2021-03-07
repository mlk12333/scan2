//
//  main.cpp
//  pscan- index
//
//  Created by 孟令凯 on 2021/2/23.
//

#include <iostream>
#include "string"
#include <fstream>
#include "Graph.hpp"

using namespace std;


int main(int argc, const char * argv[]) {
    // insert code here...
    if(argc<5) return 0;
    std::cout << "Hello, World!\n";
    
//    string str = "/Users/milk/test_data/wiki-Talk.txt";
    
    string f = "/Users/milk/test_data/index/";
    
//    string str2 = "0.4";
//    string str3 = "0.4";
//    string str4 = "4";
//    argv[2] = &str2[0];
//    argv[3] = &str3[0];
//    argv[4] = &str4[0];
    
    clock_t startTime,endTime;
 
    startTime = clock();//计时开始
    Graph graph(argv[1]);
    graph.read_graph();
    graph.pSCAN(argv[2], argv[3], atoi(argv[4]));
    graph.output();
    
    endTime = clock();//计时结束
    cout << "The time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    return 0;
}

