#include "problems.h"

int main (int argc, char *argv[]){

    //N : # of Nodes
    //M : # of Factors
    //p : # of variables that enter in a factor
    
    int N, M, p, seed, TT;

    
    if(argc==6){
        int i = 1;
        N     = atoi(argv[i++]);
        M     = atoi(argv[i++]);
        p     = atoi(argv[i++]);
        seed  = atoi(argv[i++]);
        TT    = atoi(argv[i++]);
    }
    else{
        cout << "argument: N, M, p, seed, TT" << endl;
        return 0;
    }
    
    
    srand(seed);
    
    Graph G(N,p);
    //f_Instance(G, N, M, p);
    //f_plantedInstance(G, N, M, p);
    //f_BPsweep(G, N, M, p);
    //f_BPiteration(G, N, M, p);
    f_BPguidedDecimation(G, N, M, p, TT);
    
    return 1;
}
