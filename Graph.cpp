//
//  Graph.cpp
//  scan2
//
//  Created by 孟令凯 on 2021/3/4.
//

#include "Graph.hpp"
#include <vector>

Graph::Graph(string str,string f){
    this->str = str;
    this->f = f;
}

void Graph::creatIndex(double parameter1,double parameter2, int parameterNum){
    ifstream infile;   //输入流
    int u,v;
    
    n=0;
    
    infile.open(str, ios::in);
    if (!infile.is_open())
        cout<<"Open file failure"<<endl;
    while (!infile.eof())            // 若未到文件结束一直循环
    {
        infile >> u >> v;
        m++;
        if(u>n) n = u;
        if(v>n) n = v;
        out_edges[u].push_back(v);
        in_edges[v].push_back(u);
        
//        out_edges[v].push_back(u);
//        in_edges[u].push_back(v);
    }
    infile.close();
    n++;
    out_degree.resize(n);
    in_degree.resize(n);
    for(int i = 0;i<n;i++){
        out_degree[i] = (int)out_edges[i].size();
        in_degree[i] = (int)in_edges[i].size();
    }
    
    clock_t startTime,endTime;
    
    startTime = clock();//计时开始
    
    std::cout << "read!\n";
    endTime = clock();//计时结束
    cout << "The time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    neiNum.resize(n);
    
    getNeiNum();
    
    std::cout << "getNeiNum\n";
    endTime = clock();//计时结束
    cout << "The time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    
//    getTriangles();
    getCycleTti();
    getFlowTti();
    
    std::cout << "getCycleTti!\n";
    endTime = clock();//计时结束
    cout << "The time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    
//    rankByDegree();
//    cout<<counting()<<endl;
    
    buildNeighborOrder();
    
    std::cout << "buildNeighborOrder!\n";
    endTime = clock();//计时结束
    cout << "The time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    
    
//    writeNOrder();
//    writeCOrder();
    writeCOrder2(parameterNum);
    
    std::cout << "riteCOrder!\n";
    endTime = clock();//计时结束
    cout << "The time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

    query(parameter1, parameter2, parameterNum);
    
    std::cout << "query!\n";
    
    endTime = clock();//计时结束
    cout << "The time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

}

void Graph::readIndex(){
    ifstream infile;   //输入流
    infile.open(f+"nei_num.txt", ios::in);
    if (!infile.is_open())
        cout<<"Open file failure"<<endl;
    infile>>n;
//    neiNum = new int(n);
    neiNum.resize(n);
    int i = 0;
    while (!infile.eof())            // 若未到文件结束一直循环
    {
        infile>>neiNum[i];
        i++;
    }
    infile.close();
    
    infile.open(f+"n-order1.txt", ios::in);
    int a,b;
    while (!infile.eof())            // 若未到文件结束一直循环
    {
        if(infile>>a>>b){
            vector<vector<double>> single;
            
            for(int j = 0;j<b;j++){
                vector<double> scare(2,0);
                infile>>scare[0]>>scare[1];
//                cout<<scare[0]<<" "<<scare[1]<<endl;
                single.push_back(scare);
            }
            
            neighbor_order1.push_back(single);
        }
    }
    infile.close();
    
    infile.open(f+"n-order2.txt", ios::in);
    
    while (!infile.eof())            // 若未到文件结束一直循环
    {
        if(infile>>a>>b){
            vector<vector<double>> single;
            
            for(int j = 0;j<b;j++){
                vector<double> scare(2,0);
                infile>>scare[0]>>scare[1];
//                cout<<scare[0]<<" "<<scare[1]<<endl;
                single.push_back(scare);
            }
            
            neighbor_order2.push_back(single);
        }
    }
    infile.close();
    
    int corei = 2;
    
    while (true) {
        fstream file;
        file.open(f+"c-order1-"+to_string(corei)+".txt",ios::in);
        if(!file){
            break;
        }
        
        vector<vector<double>> single;
        
        while (!file.eof())            // 若未到文件结束一直循环
        {
            vector<double> scare(2,0);
            if(file>>scare[0]>>scare[1]) single.push_back(scare);
        }
        core_order1.push_back(single);
        
        file.close();
        file.open(f+"c-order2-"+to_string(corei)+".txt",ios::in);
        if(!file){
            cout<<"文件不存在"<<endl;
            break;
        }

        vector<vector<double>> single2;
        
        while (!file.eof())            // 若未到文件结束一直循环
        {
            vector<double> scare2(2,0);
            if(file>>scare2[0]>>scare2[1]) single2.push_back(scare2);
        }
        core_order2.push_back(single2);
        file.close();
        corei++;
    }
}

void Graph::buildNeighborOrder(){
    vector<int> start(n,0);
    for(int i = 0;i<n;i++){
        for(int j = start[i];j<neiNum[i];j++){
//            double simi1 = getSimi1(i,neighbor_order1[i][j][1]);
//            double simi2 = getSimi2(i,neighbor_order2[i][j][1]);
            
            double simi1 = (cNum[i][neighbor_order1[i][j][1]] + 4)/(2*pow((neiNum[i]+1)*(neiNum[neighbor_order1[i][j][1]]+1), 0.5));
            double simi2 = (fNum[i][neighbor_order2[i][j][1]] + 12)/(6*pow((neiNum[i]+1)*(neiNum[neighbor_order2[i][j][1]]+1), 0.5));
//            double simi1 = 1;
//            double simi2 = 1;
            neighbor_order1[i][j][0] = simi1;
            neighbor_order2[i][j][0] = simi2;

            neighbor_order1[neighbor_order1[i][j][1]][start[neighbor_order2[i][j][1]]][0] = simi1;
            neighbor_order2[neighbor_order2[i][j][1]][start[neighbor_order2[i][j][1]]][0] = simi2;
            
            start[neighbor_order2[i][j][1]]++;
        }
        
        sort(neighbor_order1[i].rbegin(),neighbor_order1[i].rend());
        sort(neighbor_order2[i].rbegin(),neighbor_order2[i].rend());
    }
}

void Graph::getNeiNum(){
    m2 = 0;
    
    for(int i = 0;i<n;i++){//如果考虑有空下的点则不能这样写
//        cout<<i<<endl;
       
        int num = 0;
        
        int l1;
        int l2;
        
        vector<vector<double>> neighbor1;
        
        vector<vector<double>> neighbor2;
        
        if(in_edges.find(i)!=in_edges.end()){
//            sort(in_edges[i].begin(), in_edges[i].end());
            l2 = (int)in_edges[i].size();
        }else l2 = 0;
        
        
        if(out_edges.find(i)!=out_edges.end()){
            l1 = (int)out_edges[i].size();
        }else l1 = 0;
        
        int i1 = 0, i2 = 0;
        
        
        
        while(i1<l1 && i2<l2){
            vector<double> single1,single2;
            
            if(out_edges[i][i1] == in_edges[i][i2]){
                
                
                
                single1.push_back(-1);
                single2.push_back(-1);
                single1.push_back(out_edges[i][i1]);
                single2.push_back(out_edges[i][i1]);
                
                i1++;
                i2++;
                num++;
            }else if(out_edges[i][i1] < in_edges[i][i2]){
                
               
                
                single1.push_back(-1);
                single2.push_back(-1);
                single1.push_back(out_edges[i][i1]);
                single2.push_back(out_edges[i][i1]);
                i1++;
                num++;
            }else{
                single1.push_back(-1);
                single2.push_back(-1);
                single1.push_back(in_edges[i][i2]);
                single2.push_back(in_edges[i][i2]);
                i2++;
                num++;
            }
            
            neighbor1.push_back(single1);
            neighbor2.push_back(single2);
        }
        while(i1<l1){
            
            
            
            vector<double> single1,single2;
            single1.push_back(-1);
            single2.push_back(-1);
            single1.push_back(out_edges[i][i1]);
            single2.push_back(out_edges[i][i1]);
            neighbor1.push_back(single1);
            neighbor2.push_back(single2);
            i1++;
            num++;
        }
        while(i2<l2){
            
           
            
            vector<double> single1,single2;
            single1.push_back(-1);
            single2.push_back(-1);
            single1.push_back(in_edges[i][i2]);
            single2.push_back(in_edges[i][i2]);
            neighbor1.push_back(single1);
            neighbor2.push_back(single2);
            i2++;
            num++;
        }
        
        neighbor_order1.push_back(neighbor1);
        neighbor_order2.push_back(neighbor2);
        
        neiNum[i] = num;
    }
    
//    ofstream fout(f+"nei_num.txt",ios::out); //创建待写入数据文件
//    fout<<n<<"\n";
//    for(int j = 0;j<n;j++){
//        fout<<neiNum[j]<<"\n";
//    }
//    fout.close();
    
}
int Graph::BinarySearch(vector<int> &a, int value, int n)
{
    int low, high, mid;
    low = 0;
    high = n-1;
    while(low<=high)
    {
        mid = (low+high)/2;
        if(a[mid]==value)
            return 1;
        if(a[mid]>value)
            high = mid-1;
        if(a[mid]<value)
            low = mid+1;
    }
    return 0;
}

void Graph::writeNOrder(){
    ofstream fout(f+"n-order1.txt",ios::out); //创建待写入数据文件
    
    for(int i = 0; i < n; ++i) {
//        fout.width(2);  //设定宽度为2，默认右对齐
        fout<<i<<" "<<neiNum[i]<<"\n";  //依次写入数据，其他类型原理相同
        
        for(int j = 0;j<neiNum[i];j++){
            fout<<neighbor_order1[i][j][0]<<" "<<neighbor_order1[i][j][1]<<"\n";
        }
    }
    
    fout.close();  //关闭文件，写入成功
    
    ofstream fout2(f+"n-order2.txt",ios::out); //创建待写入数据文件
    
    for(int i = 0; i < n; ++i) {
//        fout.width(2);  //设定宽度为2，默认右对齐
        fout2<<i<<" "<<neiNum[i]<<"\n";  //依次写入数据，其他类型原理相同
        
        for(int j = 0;j<neiNum[i];j++){
            fout2<<neighbor_order2[i][j][0]<<" "<<neighbor_order2[i][j][1]<<"\n";
        }
    }
    
    fout2.close();  //关闭文件，写入成功
}

void Graph::writeCOrder(){
    int i = 2;
    while (true) {
        vector<vector<double>> c;
        for(int j=0;j<n;j++){
            if(neiNum[j]>=i-1){
                vector<double> single;
                single.push_back(neighbor_order1[j][i-2][0]);
                single.push_back(j);
                c.push_back(single);
            }
        }
        if(c.size() == 0) break;
        
        ofstream fout(f+"c-order1-"+to_string(i)+".txt",ios::out); //创建待写入数据文件
        sort(c.rbegin(), c.rend());
        core_order1.push_back(c);
        for(int j = 0;j<c.size();j++){
            fout<<c[j][0]<<" "<<c[j][1]<<"\n";
        }
        fout.close();
        i++;
    }
    i=2;
    
    while (true) {
        vector<vector<double>> c;
        for(int j=0;j<n;j++){
            if(neiNum[j]>=i-1){
                vector<double> single;
                single.push_back(neighbor_order2[j][i-2][0]);
                single.push_back(j);
                c.push_back(single);
            }
        }
        if(c.size() == 0) break;
        
        ofstream fout(f+"c-order2-"+to_string(i)+".txt",ios::out); //创建待写入数据文件
        sort(c.rbegin(), c.rend());
        
        core_order2.push_back(c);
        
        for(int j = 0;j<c.size();j++){
            fout<<c[j][0]<<" "<<c[j][1]<<"\n";
        }
        fout.close();
        i++;
    }

}

void Graph::writeCOrder2(int parameterNum){
    int i = parameterNum;
//    while (true) {
        vector<vector<double>> c;
        for(int j=0;j<n;j++){
            if(neiNum[j]>=i-1){
                vector<double> single;
                single.push_back(neighbor_order1[j][i-2][0]);
                single.push_back(j);
                c.push_back(single);
            }
        }
        if(c.size() == 0) return;
        sort(c.rbegin(), c.rend());
        core_order1.push_back(c);
//        i++;
//    }
//    i=2;
//
//    while (true) {
        vector<vector<double>> c1;
        for(int j=0;j<n;j++){
            if(neiNum[j]>=i-1){
                vector<double> single;
                single.push_back(neighbor_order2[j][i-2][0]);
                single.push_back(j);
                c1.push_back(single);
            }
        }
        if(c1.size() == 0) return;
        sort(c1.rbegin(), c1.rend());
        core_order2.push_back(c1);
//        i++;
//    }

}

void Graph::query(double parameter1,double parameter2, int parameterNum){
    vector<int> core = getAllCore(parameter1, parameter2, parameterNum-1);
    
    int ll = (int)core.size();
    
    vector<bool> allv(n,true);
    
    hash_map<int,vector<vector<int>>> cluster;
    
    int clusterNum = 0;
    
    for(int it = 0;it<ll;it++){
        
        vector<int> n_core;
        
        if(!allv[core[it]]) continue;
        clusterNum++;
        allv[core[it]] = false;
        
        queue<int> q;
        q.push(core[it]);
        
        while(!q.empty()){
            int m = q.front();
            vector<int> single;
            
            single.push_back(m);
            single.push_back(1);
            cluster[clusterNum].push_back(single);
            
            q.pop();
            vector<int> m_nei = getAllNei(parameter1, parameter2, m);
            int l = (int)m_nei.size();
            for(int i = 0;i<l;i++){
                if(!allv[m_nei[i]]) continue;
                if(BinarySearch(core,m_nei[i],(int)core.size())){//可以改进
                    q.push(m_nei[i]);
                    allv[m_nei[i]] = false;
                }else{
                    allv[m_nei[i]] = false;
                    
                    n_core.push_back(m_nei[i]);
                    
                    vector<int> single2;
                    
                    single2.push_back(m_nei[i]);
                    single2.push_back(0);
                    cluster[clusterNum].push_back(single2);
                }
            }
        }
        
        int n_corel = (int)n_core.size();
        for(int n_corei = 0;n_corei<n_corel;n_corei++){
            allv[n_core[n_corei]] = true;
        }
    }
    
//    ofstream fout(f+"/output/"+to_string(parameter1)+"-"+to_string(parameter2)+"-"+to_string(parameterNum)+".txt",ios::out); //创建待写入数据文件
//    fout<<"c/n"<<" "<<"v"<<" "<<"clu"<<"\n";
    
    cout<<"c/n"<<" "<<"v"<<" "<<"clu"<<"\n";

    for(hash_map<int,vector<vector<int>>>::iterator it = cluster.begin(); it != cluster.end();++it){
        
        int clul = (int)it->second.size();
        
//        fout<<it->first<<" num:........................."<<clul<<"\n";
        
        cout<<it->first<<" num:........................."<<clul<<"\n";

        
        for(int i = 0;i<clul;i++){
            string cc;
            if(it->second[i][1]) cc="c";
            else cc = "n";
//            fout<<cc<<" "<<to_string(it->second[i][0])<<" "<<it->first<<"\n";
            cout<<cc<<" "<<to_string(it->second[i][0])<<" "<<it->first<<"\n";
        }
    }
//    fout.close();
    
}

vector<int> Graph::getAllNei(double parameter1,double parameter2,int v){
    int l = neiNum[v];
    vector<int> single;
    vector<int> out;
    for(int i = 0;i<l;i++){
        if(neighbor_order1[v][i][0]<parameter1) break;
        single.push_back(neighbor_order1[v][i][1]);
    }
    
    for(int i = 0;i<l;i++){
        if(neighbor_order2[v][i][0]<parameter2) break;
        single.push_back(neighbor_order2[v][i][1]);
    }
    
    sort(single.begin(),single.end());
    
    int l2 = (int)single.size();
    
    for(int i = 0;i<l2-1;i++){
        if(single[i] == single[i+1]){
            out.push_back(single[i]);
            i++;
        }
    }
    return out;
}

vector<int> Graph::getAllCore(double parameter1,double parameter2, int corei){
    int l = (int)core_order1[0].size();
    vector<int> single;
    vector<int> out;
    for(int i = 0;i<l;i++){
        if(core_order1[0][i][0]<parameter1) break;
        single.push_back(core_order1[0][i][1]);
    }
    
    for(int i = 0;i<l;i++){
        if(core_order2[0][i][0]<parameter2) break;
        single.push_back(core_order2[0][i][1]);
    }
    
    sort(single.begin(),single.end());
    
    int l2 = (int)single.size();
    
    for(int i = 0;i<l2-1;i++){
        if(single[i] == single[i+1]){
            out.push_back(single[i]);
            i++;
        }
    }
    return out;
}



void Graph::getTriangles(){
    
    for(int i=0;i<n;i++){
        hash_map<int,int> neighbors;
        for(int j = 0;j<out_degree[i];j++){
            neighbors[out_edges[i][j]] = out_edges[i][j];
//            cout<<out_edges[i][j]<<"  ";
        }
//        for(int j = 0;j<in_degree[i];j++){
//            neighbors2[in_edges[i][j]] = in_edges[i][j];
//        }
        
        for(int j = 0;j<out_degree[i];j++){
//            cout<<out_degree[out_edges[i][j]]<<" "<<out_edges[i][j]<<" "<<out_degree[i]<<endl;
            if((out_degree[out_edges[i][j]]<out_degree[i])||(out_degree[out_edges[i][j]] == out_degree[i] && out_edges[i][j]<i)){
                for(int k = 0;k<out_degree[out_edges[i][j]];k++){
                    if(neighbors.find(out_edges[out_edges[i][j]][k]) != neighbors.end()){
                        
//                        cout<<"sanjiaoxing1:"<<i<<" "<<out_edges[i][j]<<" "<<out_edges[out_edges[i][j]][k]<<endl;
                        int a = i;
                        int b =out_edges[i][j];
                        int c =out_edges[out_edges[i][j]][k];
                        fNum[a][b]++;
                        fNum[a][c]++;
                        fNum[b][a]++;
                        fNum[b][c]++;
                        fNum[c][b]++;
                        fNum[c][a]++;
                    }
                }
            }
        }
        
        for(int j = 0;j<in_degree[i];j++){
            if((out_degree[in_edges[i][j]]<out_degree[i])||(out_degree[in_edges[i][j]] == out_degree[i] && in_edges[i][j]<i)){
                
                for(int k = 0;k<out_degree[i];k++){
                    
                }
                for(int k = 0;k<out_degree[in_edges[i][j]];k++){
                    if(neighbors.find(out_edges[in_edges[i][j]][k]) != neighbors.end()){
//                        cout<<"sanjiaoxing2:"<<i<<" "<<in_edges[i][j]<<" "<<out_edges[in_edges[i][j]][k]<<endl;
                        int a = i;
                        int b =in_edges[i][j];
                        int c =out_edges[in_edges[i][j]][k];
                        
                        fNum[a][b]++;
                        fNum[a][c]++;
                        fNum[b][a]++;
                        fNum[b][c]++;
                        fNum[c][b]++;
                        fNum[c][a]++;
                    }
                }
            }
        }
    }
    
}

bool Graph::judgeCycletrai(int a,int b,int c){
    return true;
}

void Graph::getCycleTti(){
    for(int i = 0;i<n;i++){
        for(int j =0;j<out_degree[i];j++){
            if(i<out_edges[i][j]){
                vector<int> commonP = intersection(in_edges[i], out_edges[out_edges[i][j]], in_degree[i], out_degree[out_edges[i][j]]);
                int l =(int)commonP.size();
                for(int k = 0;k<l;k++){
                    if(i<commonP[k]){
//                        cout<<i<<" "<<out_edges[i][j]<<" "<<commonP[k]<<endl;
                        int a = i;
                        int b =out_edges[i][j];
                        int c =commonP[k];

                        cNum[a][b]++;
                        cNum[a][c]++;
                        cNum[b][a]++;
                        cNum[b][c]++;
                        cNum[c][b]++;
                        cNum[c][a]++;
                                             
                    }
                }
            }
        }
    }
    
}

void Graph::getFlowTti(){
    for(int i = 0;i<n;i++){
        for(int j =0;j<out_degree[i];j++){
                vector<int> commonP = intersection(out_edges[i], out_edges[out_edges[i][j]], out_degree[i], out_degree[out_edges[i][j]]);
            int l =(int)commonP.size();
                for(int k = 0;k<l;k++){
//                        cout<<i<<" "<<out_edges[i][j]<<" "<<commonP[k]<<endl;
                        int a = i;
                        int b =out_edges[i][j];
                        int c =commonP[k];
                    
                    fNum[a][b]++;
                    fNum[a][c]++;
                    fNum[b][a]++;
                    fNum[b][c]++;
                    fNum[c][b]++;
                    fNum[c][a]++;
                }
            
        }
    }
}

vector<int> Graph::intersection(vector<int> &A, vector<int> &B, int a, int b){
    vector<int> out;
    int i = 0, j =0;
    while(i<a && j<b){
        if(A[i] == B[j]){
            out.push_back(A[i]);
            i++;
            j++;
        }else if(A[i]<B[j]){
            i++;
        }else{
            j++;
        }
    }
    return out;
}

