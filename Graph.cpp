//
//  Graph.cpp
//  pscan- index
//
//  Created by 孟令凯 on 2021/2/23.
//

#include "Graph.hpp"
Graph::Graph(const string _dir) {
    dir = _dir;

    n = in_m = out_m =m=degree_m= 0;

    eps_a1 = eps_b1 = eps_a2 = eps_b2 = miu = 0;
    
}

void Graph::read_graph(){
    ifstream infile;   //输入流
    int u,v;
    
    infile.open(dir, ios::in);
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
        
    }
    infile.close();
    n++;
    out_degree.resize(n);
    in_degree.resize(n);
    degree.resize(n);
    pstart.resize(n+1);
    in_pstart.resize(n+1);
    out_pstart.resize(n+1);
    in_pstart[0] = 0;
    out_pstart[0] = 0;


    for(int i = 0;i<n;i++){
        out_degree[i] = (int)out_edges[i].size();
        in_degree[i] = (int)in_edges[i].size();
        
        in_pstart[i+1] = in_pstart[i] + in_degree[i];
        out_pstart[i+1] = out_pstart[i] + out_degree[i];
        
        for(int j = 0;j<in_degree[i];j++) in_edges1.push_back(in_edges[i][j]);
        for(int j = 0;j<out_degree[i];j++) out_edges1.push_back(out_edges[i][j]);
    }
    
//    for(int i = 0;i<m;i++) cout<<in_edges1[i]<<" ";
//    cout<<endl;
//
//    for(int i = 0;i<m;i++) cout<<out_edges1[i]<<" ";
//    cout<<endl;

    getNeiNum();
}

void Graph::getNeiNum(){
    
    pstart[0] = 0;
    
    for(int i = 0;i<n;i++){//如果考虑有空下的点则不能这样写
//        cout<<i<<endl;
       
        int num = 0;
        
        int l1;
        int l2;
        if(in_edges.find(i)!=in_edges.end()){
//            sort(in_edges[i].begin(), in_edges[i].end());
            l2 = (int)in_edges[i].size();
        }else l2 = 0;
        
        
        if(out_edges.find(i)!=out_edges.end()){
            l1 = (int)out_edges[i].size();
        }else l1 = 0;
        
        int i1 = 0, i2 = 0;
        
        while(i1<l1 && i2<l2){
            if(out_edges[i][i1] == in_edges[i][i2]){
                edges.push_back(out_edges[i][i1]);
                i1++;
                i2++;
                num++;
            }else if(out_edges[i][i1] < in_edges[i][i2]){
                edges.push_back(out_edges[i][i1]);
                i1++;
                num++;
            }else{
                edges.push_back(in_edges[i][i2]);
                i2++;
                num++;
            }        }
        while(i1<l1){
            edges.push_back(out_edges[i][i1]);
            i1++;
            num++;
        }
        while(i2<l2){
            edges.push_back(in_edges[i][i2]);

            i2++;
            num++;
        }
        degree[i] = num;
        pstart[i+1] = pstart[i] + degree[i];
        ++ degree[i];
    }
    
    degree_m = (int)edges.size();
    reverse.resize(degree_m);
    min_cn1.resize(degree_m);
    min_cn2.resize(degree_m);
    
//    cNum.resize(degree_m,4);
//    fNum.resize(degree_m,12);
    
//    for(int i = 0;i<degree_m;i++) cout<<edges[i]<<" ";
//    cout<<endl;
}

void  Graph::pSCAN(const char *eps_s1,const char *eps_s2, int _miu){
    get_eps(eps_s1,eps_s2);
    miu = _miu;
    
    similar_degree.resize(n,0);
    effective_degree.resize(n);
    for(int i = 0;i < n;i ++) effective_degree[i] = degree[i]-1;
    
    pa.resize(n);//pa和rank用来找core
    rank.resize(n);
    for(int i = 0;i < n;i ++) {
        pa[i] = i;
        rank[i] = 0;
    }
    
    int *edge_buf = new int[n];
    int *cores = new int[n];
    int cores_n = 0;
    
    
    clock_t startTime,endTime;
    
    startTime = clock();//计时开始
   
    prune_and_cross_link(eps_a1, eps_b1, eps_a2, eps_b2, miu, cores, cores_n);
    printf("\t*** Finished prune and cross link!\n");
    
    endTime = clock();//计时结束
    cout << "The prune_and_cross_link time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
    
    
//    for(int corei=0;corei<n;corei++) cout<<cores[corei]<<" ";
//    cout<<endl;

    int *bin_head = new int[n];//这两个方便排序
    int *bin_next = new int[n];
    for(int i = 0;i < n;i ++) bin_head[i] = -1;
    
//    for(int i = 0;i < n;i ++) cout<<effective_degree[i]<<endl;
    
    int max_ed = 0;
    for(int i = 0;i < n;i++) if(effective_degree[i] >= miu) {
        int ed = effective_degree[i];
        if(ed > max_ed) max_ed = ed;
        bin_next[i] = bin_head[ed];
        bin_head[ed] = i;
    }
    
//    for(int i = 0;i < n;i ++) cout<<bin_head[i]<<"  ";
//    cout<<endl;
//    for(int i = 0;i < n;i ++) cout<<bin_next[i]<<"  ";
//    cout<<endl;
    
    
    while(true) {
//        for(int corei=0;corei<n;corei++) cout<<cores[corei]<<" ";
//        cout<<endl;
        int u = -1;
        if(cores_n) u = cores[-- cores_n];
        else {
            while(max_ed >= miu&&u == -1) {
                for(int x = bin_head[max_ed];x != -1;) {
                    int tmp = bin_next[x];
                    int ed = effective_degree[x];
                    if(ed == max_ed) {
                        u = x;
                        bin_head[max_ed] = bin_next[x];
                        break;
                    }
                    else if(ed >= miu) {
                        bin_next[x] = bin_head[ed];//超了，强制变回max_ed
                        bin_head[ed] = x;
                    }
                    x = tmp;
                }
                if(u == -1) {//度最大的找完了
                    bin_head[max_ed] = -1;
                    -- max_ed;
                }
            }
        }
        if(u == -1) break;//最大的度已经不满足要求了
        
//        vector<vector<int>> edge_buf2;
        
        int edge_buf_n = 0;
        for(int j = pstart[u+1]-1;j >= pstart[u];j --) {//查看当前点的所有邻居
            if(min_cn1[j] == -2) continue;
//            if(find_root(u) != find_root(edges[j])) edge_buf[edge_buf_n ++] = j;
            if(similar_degree[u] < miu||find_root(u) != find_root(edges[j])){
                
                
                edge_buf[edge_buf_n++] = j;
//                vector<int> edge_buf3;
//                edge_buf3.push_back(degree[j]);
//                edge_buf3.push_back(j);
//
//                edge_buf2.push_back(edge_buf3);
            }             //edge_buf存储还没有确定关系的点以及他的邻居，edge_buf_n是个数
        }
        
//        int edge_buf_n = edge_buf2.size();
//        sort(edge_buf2.rbegin(), edge_buf2.rend());


        int i = 0;
        while(similar_degree[u] < miu&&effective_degree[u] >= miu&&i < edge_buf_n) {//还有机会相似的点
            int idx = edge_buf[i];//u的邻居的一部分
            if(min_cn1[idx] != -1) {//-1已经合格 不许再谈
#ifdef _DEBUG_
                if(min_cn1[idx] == 0) printf("WA min_cn!\n");
#endif
                int v = edges[idx];//具体的点
                
                
                
                
                min_cn1[idx] = min_cn1[reverse[idx]] = min_cn2[idx] = min_cn2[reverse[idx]] = similar_check_OP(u, idx, eps_a1, eps_b1, eps_a2, eps_b2);
//                cout<<min_cn1[idx]<<endl;
                
                if(min_cn1[idx] == -1) ++ similar_degree[u];
                else -- effective_degree[u];

                if(effective_degree[v] >= 0) {
                    if(min_cn1[idx] == -1) {
                        ++ similar_degree[v];

                        if(similar_degree[v] == miu) cores[cores_n ++] = v;
                    }
                    else -- effective_degree[v];
                }
            }

            ++ i;
        }

        effective_degree[u] = -1;//标记为已检查过的点

        if(similar_degree[u] < miu) continue;//没希望成为核心了

        for(int j = 0;j < edge_buf_n;j++) {
            int idx = edge_buf[j];
            if(min_cn1[idx] == -1&&similar_degree[edges[idx]] >= miu)
                my_union(u, edges[idx]);//核心抱团
        }

        while(i < edge_buf_n) {//检查剩余的点
            int idx = edge_buf[i];
            int v = edges[idx];
            if(min_cn1[idx] < 0||similar_degree[v] < miu||find_root(u) == find_root(v)) {//没有确定相似，但已经是核心并且没抱在一起才能跳过continue，是核心，但是还不确定想不相思
                ++ i;
                continue;
            }

            min_cn1[idx] = min_cn1[reverse[idx]] = min_cn2[idx] = min_cn2[reverse[idx]] = similar_check_OP(u, idx, eps_a1, eps_b1, eps_a2, eps_b2);


            if(effective_degree[v] >= 0) {//说明还没有检查过
                if(min_cn1[idx] == -1) {
                    ++ similar_degree[v];

                    if(similar_degree[v] == miu) cores[cores_n ++] = v;
                }
                else -- effective_degree[v];
            }

            if(min_cn1[idx] == -1) my_union(u,v);//已经是核心并且相似 那就抱一起吧

            ++ i;
        }
        //printf(")\n");
    }

    delete[] edge_buf; edge_buf = NULL;
    delete[] cores; cores = NULL;
    delete[] bin_head; bin_head = NULL;
    delete[] bin_next; bin_next = NULL;
    
    cluster_noncore_vertices(eps_a1,eps_b1,eps_a2,eps_b2,miu);
}


void Graph::get_eps(const char *eps_s1,const char *eps_s2) {
    int i = 0, eps_a = 0, eps_b = 1;
    while(eps_s1[i] != '\0'&&eps_s1[i] != '.') {
        eps_a = eps_a*10 + (eps_s1[i]-'0');
        ++ i;
    }

    if(eps_s1[i] == '.') {
        ++ i;
        while(eps_s1[i] != '\0') {
            eps_a = eps_a*10 + (eps_s1[i]-'0');
            eps_b *= 10;
            ++ i;
        }
    }

    if(eps_a > eps_b||eps_b > 100||eps_a <= 0) {
        printf("??? Wrong eps format: %d/%d\n", eps_a, eps_b);
        exit(1);
    }

    eps_a1 = eps_a * eps_a;
    eps_b1 = eps_b * eps_b;
//    cout<<eps_a1<<" "<<eps_b1<<endl;
    
    i = 0;
    eps_a = 0;
    eps_b = 1;
   while(eps_s2[i] != '\0'&&eps_s2[i] != '.') {
       eps_a = eps_a*10 + (eps_s2[i]-'0');
       ++ i;
   }

   if(eps_s2[i] == '.') {
       ++ i;
       while(eps_s2[i] != '\0') {
           eps_a = eps_a*10 + (eps_s2[i]-'0');
           eps_b *= 10;
           ++ i;
       }
   }

   if(eps_a > eps_b||eps_b > 100||eps_a <= 0) {
       printf("??? Wrong eps format: %d/%d\n", eps_a, eps_b);
       exit(1);
   }

   eps_a2 = eps_a * eps_a;
   eps_b2 = eps_b * eps_b;
//    cout<<eps_a2<<" "<<eps_b2<<endl;

}

void Graph::prune_and_cross_link(int eps_a1, int eps_b1,int eps_a2, int eps_b2, int miu, int *cores, int &cores_e){
    
    for(int i = 0;i < n;i ++) { //must be iterating from 0 to n-1
        for(int j = pstart[i];j < pstart[i+1];j ++) {//i点所有邻居
            if(edges[j] < i) {
                if(min_cn1[j] == 0) min_cn1[j] = -2;
                if(min_cn2[j] == 0) min_cn2[j] = -2;

                continue; //this edge has already been checked
            }
            
            int v = edges[j];//v是i的邻居
            int a = degree[i], b = degree[v];//此时度是加一的
            if(a > b) swap(a, b);
            
            if((((long long)a)*eps_b2 < ((long long)b)*eps_a2) || (((long long)a)*eps_b1 < ((long long)b)*eps_a1)) {//近似计算
                min_cn1[j] = -2;
                min_cn2[j] = -2;

                -- effective_degree[i];
                -- effective_degree[v];
            }
            else {
                
                int c1 = ceil(sqrtl(((long double)a*(long double)b*(long double)eps_a1)/(long double)eps_b1) * (long double)2);

                int c2 = ceil(sqrtl(((long double)a*(long double)b*(long double)eps_a2)/(long double)eps_b2) * (long double)6);

                if(c1 <= 4 && c2<=12) {//直接认定相似
                    min_cn1[j] = -1;
                    min_cn2[j] = -1;
                    
                    ++ similar_degree[i];
                    ++ similar_degree[v];
                    
                    if(similar_degree[i] == miu)
                        cores[cores_e ++] = i;
                    if(similar_degree[v] == miu)
                        cores[cores_e ++] = v;
                }
                else{
                    min_cn1[j] = c1;
                    min_cn2[j] = c2;
                }
                
            }
            
            if(min_cn1[j] != -2) {
                int r_id = binary_search(edges, pstart[v], pstart[v+1], i);
                reverse[j] = r_id;//记录反向边的位置
                reverse[r_id] = j;

                min_cn1[r_id] = min_cn1[j];//最低要求是一样的
                min_cn2[r_id] = min_cn2[j];//最低要求是一样的
            }
        }
    }
}

int Graph::binary_search(vector<int> &array, int b, int e, int val) {
#ifdef _DEBUG_
    if(e < b) printf("??? WA1 in binary_search\n");
#endif
    -- e;
    if(array[e] < val) return e+1;
    while(b < e) {
        int mid = b + (e-b)/2;
        if(array[mid] >= val) e = mid;
        else b = mid+1;
    }
#ifdef _DEBUG_
    if(array[e] < val) printf("??? WA2 in binary_search\n");
#endif
    return e;
}

int Graph::find_root(int u) {
    int x = u;
    while(pa[x] != x) x = pa[x];
//    cout<<endl;
//    for(int iu = 0;iu<n;iu++) cout<<pa[iu]<<" ";
//    cout<<endl;
//    for(int iu = 0;iu<n;iu++) cout<<rank[iu]<<" ";
//    cout<<endl;

    while(pa[u] != x) {
        int tmp = pa[u];
        pa[u] = x;
        u = tmp;
    }

    return x;
}


int Graph::check_common_neighbor(int u, int v, int c1,int c2) {
    
    int cn1 = 4, cn2 = 12;
    

    if(degree[u] > degree[v]) swap(u,v);

    int du = degree[u]+1, dv = degree[v]+1;
    int i = pstart[u], j = pstart[v];
    while(i < pstart[u+1]&&j < pstart[v+1]&&(cn1 < c1 || cn2 < c2)&&(2*du>=c1&&2*dv>=c1)&&(6*du>=c2&&6*dv>=c2)) {
        if(edges[i] < edges[j]) {
            -- du;
            ++ i;
        }
        else if(edges[i] > edges[j]) {
            -- dv;
            ++ j;
        }
        else {
            
            vector<int> cn11 = get_f_and_c_num(u,v,edges[j]);

            cn1 = cn1 + cn11[0];
            cn2 = cn2 + cn11[1];
            ++ i;
            ++ j;
        }
    }

    if(cn1 >= c1&&cn2 >= c2) return -1;
    return -2;
}



int Graph::similar_check_OP(int u, int idx, int eps_a1, int eps_b1, int eps_a2, int eps_b2) {
    int v = edges[idx];

#ifdef _DEBUG_
    if(min_cn1[idx] == -1||min_cn1[idx] == -2) printf("??? WA in similar_check\n");
#endif

    if(min_cn1[idx] == 0) {//肯定还没计算最低标准
        int du = degree[u], dv = degree[v];
        min_cn1[idx] = min_cn1[reverse[idx]] = ceil(sqrtl(((long double)du*(long double)dv*(long double)eps_a1)/(long double)eps_b1) * (long double)2);
        min_cn2[idx] = min_cn2[reverse[idx]] = ceil(sqrtl(((long double)du*(long double)dv*(long double)eps_a2)/(long double)eps_b2) * (long double)6);
    

#ifdef _DEBUG_
        if(du < c||dv < c) return -2;
#endif

        if(min_cn1[idx] <= 4 && min_cn2[idx]<=12) return -1;
        
    }
    
    return check_common_neighbor(u, v, min_cn1[idx], min_cn2[idx]);
}

vector<int> Graph::get_f_and_c_num(int a,int b, int c){
//    int a_to_b = 1;
//    int b_to_a = 1;
//    int a_to_c = 1;
//    int c_to_a = 1;
//    int c_to_b = 1;
//    int b_to_c = 1;
    
    int a_to_b = BinarySearch(out_edges[a],b,(int)out_edges[a].size());
    int b_to_a = BinarySearch(out_edges[b],a,(int)out_edges[b].size());
    int a_to_c = BinarySearch(out_edges[a],c,(int)out_edges[a].size());
    int c_to_a = BinarySearch(out_edges[c],a,(int)out_edges[c].size());
    int c_to_b = BinarySearch(out_edges[c],b,(int)out_edges[c].size());
    int b_to_c = BinarySearch(out_edges[b],c,(int)out_edges[b].size());
    
    int all =a_to_b+b_to_a+a_to_c+c_to_a+c_to_b+b_to_c;
    
    int cNum = 0;
    int fNum = 0;
    
    if(all == 6){
        cNum = cNum + 2;
        fNum = fNum + 6;    }
    else if(all == 5){
        cNum = cNum + 1;
        fNum = fNum + 3;
    }else{
        if(a_to_b && b_to_c && c_to_a) cNum++;
        if(a_to_c && c_to_b && b_to_a) cNum++;
        if(a_to_b && c_to_b && c_to_a) fNum++;
        if(a_to_b && a_to_c && c_to_b) fNum++;
        if(a_to_b && b_to_c && a_to_c) fNum++;
        if(a_to_c && b_to_c && b_to_a) fNum++;
        if(b_to_a && b_to_c && c_to_a) fNum++;
        if(b_to_a && c_to_b && c_to_a) fNum++;
    }
    vector<int> cn;
    cn.push_back(cNum);
    cn.push_back(fNum);
    return cn;
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

void Graph::my_union(int u, int v) {
    int ru = find_root(u);
    int rv = find_root(v);

    if(ru == rv) return ;

    if(rank[ru] < rank[rv]) pa[ru] = rv;
    else if(rank[ru] > rank[rv]) pa[rv] = ru;
    else {
        pa[rv] = ru;
        ++ rank[ru];
    }
}

void Graph::cluster_noncore_vertices(int eps_a1, int eps_b1,int eps_a2, int eps_b2, int mu){
    cid.resize(n,n);//用cid分离出有几个群体

    for(int i = 0;i < n;i ++) {
        
//        cout<<similar_degree[i]<<endl;
        if(similar_degree[i] >= mu) {
            int x = find_root(i);
//            cout<<cid[i]<<endl;
            if(i < cid[x]) cid[x] = i;
        }
        
    }
        
    
//    for(int cidi=0;cidi<n;cidi++) cout<<cid[cidi]<<" ";
//    cout<<endl;

    noncore_cluster.clear();
    noncore_cluster.reserve(n);
    for(int i = 0;i < n;i ++)
    if(similar_degree[i] >= mu) {
        for(int j = pstart[i];j < pstart[i+1];j ++)
        if(similar_degree[edges[j]] < mu) {//只需要找非核心的
            if(min_cn1[j] >= 0) {//还没有相似性检查
                min_cn1[j] = min_cn2[j] = similar_check_OP(i, j, eps_a1, eps_b1, eps_a2, eps_b2);
                if(reverse[reverse[j]] != j) printf("WA cluster_noncore\n");
                min_cn1[reverse[j]] = min_cn1[j];
                min_cn2[reverse[j]] = min_cn2[j];
                if(min_cn1[j] == -1) {
                    ++ similar_degree[i];
                    ++ similar_degree[edges[j]];
                }
            }

            if(min_cn1[j] == -1) noncore_cluster.push_back(make_pair(cid[pa[i]], edges[j]));
        }
    }
    
    

}

void Graph::output(){
    printf("\t*** Start write result into disk!\n");
    printf("c/n vertex_id cluster_id\n");

    int mu = miu;
    for(int i = 0;i < n;i ++) if(similar_degree[i] >= mu) {
        printf("c %d %d\n", i, cid[pa[i]]);//先输出核心
    }

    sort(noncore_cluster.begin(), noncore_cluster.end());
    noncore_cluster.erase(unique(noncore_cluster.begin(), noncore_cluster.end()), noncore_cluster.end());
    for(int i = 0;i < noncore_cluster.size();i ++) {
        printf("n %d %d\n", noncore_cluster[i].second, noncore_cluster[i].first);
    }

}

