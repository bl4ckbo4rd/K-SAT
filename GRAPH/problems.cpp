#include "problems.h"


void f_Instance(Graph& G, int N, int M, int p){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
}


void f_plantedInstance(Graph& G, int N, int M, int p){
    
    //here we explicitely define a planted configuration
    //vector<int> ps = make_vector<int>() << 0 << 0 << 1 << 1 << 1 << 0 << 0 << 1 << 1;
    //or we can create it at random
    
    vector<bool> ps;
    ps.reserve(N);
    for (int i = 0; i < N; ++i){
        if (2*(double)rand()/RAND_MAX-1 > 0)
            ps.push_back(1);
        else
            ps.push_back(0);
    }
    
    //we build an ErdosRenyi factor graph for which ps is a solution
    G.plantedErdosRenyi(M,ps);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
    
    cout << "print 1 if the planted solution satisfy the formula: " << G.check(ps) << endl;
    cout << G.check(ps) << endl;;
}


void f_BPsweep(Graph& G, int N, int M, int p){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
    
    bool verbose = 1;
    BP It(G,verbose);
    
    //we make a single BP sweep
    It.BPsweep();
    
}


void f_BPiteration(Graph& G, int N, int M, int p){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
    
    bool verbose = 1;
    BP It(G,verbose);
    
    //we iterate the BP equation until convergence
    //specifying eps and T
    double eps = 0.0001;
    int T = 100;
    
    It.BPiteration(eps, T);

    
}


void f_BPguidedDecimation(Graph& G, int N, int M, int p, int TT){
    
    //we build an ErdosRenyi factor graph
    G.ErdosRenyi(M);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
    bool verbose = 1;
    //int       TT = N;
    int       T  = 100;
    double   eps = 0.001;

    BP It(G,verbose);

    vector <int> v_bias;
    vector <bool> v_q;
    
    It.initDecimation(v_bias,v_q);
    It.BPguidedDecimation(eps,T,TT);
    
}
