#include "problems.h"

int main (int argc, char *argv[]){

    //N : # of Nodes
    //M : # of Factors
    //p : # of variables that enter in a factor
    
    int N, M, seed;
    int p=3;

    
    if(argc==4){
        int i = 1;
        N     = atoi(argv[i++]);
        M     = atoi(argv[i++]);
        seed  = atoi(argv[i++]);
    }
    else{
        cout << "argument: N, M, seed" << endl;
    }
    
    
    srand(seed);
    
    Graph G(N,p);
    f_BPsweep(G,N,M,p);

    
    return 1;
}
