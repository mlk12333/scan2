//
//  main.cpp
//  scan2
//
//  Created by 孟令凯 on 2021/3/4.
//

#include <iostream>
#include "string"
#include <fstream>
#include "Graph.hpp"
using namespace std;

int main(int argc, const char * argv[]) {
    if(argc<5) return 0;
    // insert code here...
    clock_t startTime,endTime;
    startTime = clock();//计时开始
    
    std::cout << "Hello, World!\n";
//    string str = "/Users/milk/test_data/wiki-Talk.txt";
    string str = argv[1];
    string f = "/Users/milk/test_data/index/";
    Graph graph(str,f);
    std::cout << "Hello, World!\n";
    
//    graph.creatIndex();
    graph.creatIndex(atof(argv[2]),atof(argv[3]),atoi(argv[4]));
    std::cout << "Hello, World!\n";
    return 0;
}



