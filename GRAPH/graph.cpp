#include "graph.h"



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Node
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



int Node::numberOfFactors(){
    
  return v_fac.size();
    
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Factor
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



//input variables:
//p_f : index of the factor
//p_p : # of the variables entering in a factor

Factor::Factor(int p_f, vector<bool> p_v_J, int p_p){
    
    f   = p_f;
    v_J = p_v_J;
    p   = p_p;

    v_node.reserve(p);
    
};



int Factor::numberOfNodes(){
    
    return v_node.size();
    
};



//input variables:
//i,j,k : indeces of the nodes entering in the factor

double Factor::clause(int i, int j, int k){
    
    bool J1 = v_J[0];
    bool J2 = v_J[1];
    bool J3 = v_J[2];
    
    //this is the T=0 clause
    return (2*J1 - 1) * (J1 - i) | (2*J2 - 1) * (J2 - j) | (2*J3 - 1) * (J3 - k);
    
};



//input variables:
//i,j,k : indeces of the nodes entering in the factor

//given a triplet of variables entering in a clause, we can construct several planted "clauses".
//a clause of the SAT problem is specified by three values J1, J2, J3, being 0 or 1 if the corresponding variables appear unnegated or negated in the clause,
//drawn at random in the non-planted setting.
//for instance, consider the clause x1 | !x2 | !x3 where ! indicates the negation. in this case we set J1=1, J2=0, J3=0.
//it is easy to see that the only choice of J's for which the triplet do not satisfy the clause is 0, 1, 1, and more generally, for a particular realization of (x1,x2,x3) we need
//to exclude J1=x1, J2=x2, J3=x3.

void Factor::plantedClause(int i, int j, int k){
    
    vector<bool> temp_v_J, p_J;
    bool J1, J2, J3;
    int flag=1;
    
    while (flag){

        if (2*(double)rand()/RAND_MAX-1 > 0)
            J1=1;
        else
            J1=0;
    
        if (2*(double)rand()/RAND_MAX-1 > 0)
            J2=1;
        else
            J2=0;
    
        if (2*(double)rand()/RAND_MAX-1 > 0)
            J3=1;
        else
            J3=0;

        p_J = make_vector<bool>() << i << j << k;
        temp_v_J = make_vector<bool>() << J1  << J2  << J3;
    
        if (p_J != temp_v_J)
            flag=0;
    }

    //this is the T=0 clause
    v_J = temp_v_J;
    
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Graph
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



//input variables:
//p_N : # of nodes
//p_p : # of variables entering in a factor

Graph::Graph(int p_N, int p_p) : N ( p_N ), p ( p_p ) {
  v.reserve(N);
  for (int i = 0; i < N; ++i){
    v.push_back (Node(i));
  }
};



int Graph::numberOfTotalFactors(){
  return F.size();
};



int Graph::numberOfTotalNodes(){
  return N;
};



void Graph::factorsOfNode(int i){
  cout << "node " << i << " has " << v[i].numberOfFactors() << " factors: " << endl;
  for (vector<int>::iterator it = v[i].v_fac.begin() ; it != v[i].v_fac.end(); ++it)
    cout << *it << endl;
};



void Graph::nodesOfFactor(int a){
  cout << "factor " << a << " has " << F[a]->numberOfNodes() << " nodes: " << endl;
  for (vector<int>::iterator it = F[a]->v_node.begin() ; it != F[a]->v_node.end(); ++it)
    cout << *it << endl;
};



//input variables:
//p_a    : factor index;
//p_v_J  : p=3-component vector whose components are 0 or 1 depending if the corresponding variable is negated or not
//v_da   : nodes attached to the factor

int Graph::addFactor(int p_a, vector<bool> p_v_J, vector<int> v_da){
    
    vector <int> v1 = v_da;
    std::sort(v1.begin(), v1.end());
    
    int flag=1;
    
    //before adding a factor we check that the factor does not already exist.
    //if the factor already exists, flag is set to 0.
    //to this aim it is sufficient to check that the first node of v_da does not appear in a factor with the two other nodes of v_da.
    
    for(vector<int>::iterator it_a = v[v_da[0]].v_fac.begin(); it_a != v[v_da[0]].v_fac.end(); it_a++){
        
        vector <int> v2 = F[*it_a]->v_node;
        std::sort(v2.begin(), v2.end());
        
        if (v1 == v2){
            flag=0;
        }
    }
    
    if(flag){
        
        //a is a pointer to the derived class FactorSat
        Factor* a = new Factor(p_a,p_v_J,p);
        
        for (int i=0; i<v_da.size(); i++)
            a->v_node.push_back(v_da[i]);
        
        F.push_back(a);
        
        for(int i=0; i<v_da.size(); i++)
            v[v_da[i]].v_fac.push_back(p_a);
    }
    
    return flag;
    
};



void Graph::graphStructure(){
    
    cout << endl;
    cout << endl;
    
    cout << "structure of the graph: " << endl;
    cout << endl;
    
    for (int i=0; i<N; i++)
        factorsOfNode(i);
    
    for (int a=0; a<M; a++)
        nodesOfFactor(a);
    
    cout << endl;
    cout << endl;
    
    cout << "factors types: " << endl;
    for (int a=0; a<M; a++)
        cout << a << " J's: " << F[a]->v_J[0] << " " << F[a]->v_J[1] << " " << F[a]->v_J[2] << endl;
    
    cout << endl;
    cout << endl;
    
};



//input variables:
//p_M: # of factors

void Graph::ErdosRenyi(int p_M){
    
    M = p_M;
    
    vector<int> v;
    vector <bool> v_J;
    int i,j,k;
    bool J1, J2, J3;
    int a = 0;
    int flag;
    
    while (a<M){
        i = rand() % N ;
        j = rand() % N ;
        k = rand() % N ;
        
        if (i!=j && i!=k && j!=k){
            v = make_vector<int>() << i << j << k;
            
            if (2*(double)rand()/RAND_MAX-1 > 0)
                J1=1;
            else
                J1=0;
            
            if (2*(double)rand()/RAND_MAX-1 > 0)
                J2=1;
            else
                J2=0;
            
            if (2*(double)rand()/RAND_MAX-1 > 0)
                J3=1;
            else
                J3=0;
            
            v_J = make_vector<bool>() << J1 << J2 << J3;
            
            flag=addFactor(a, v_J, v);
            if (flag) a++;
        }
    }
    
};



//input variables:
//p_M : # of factors
//ps  : planted configuration

void Graph::plantedErdosRenyi(int p_M, vector<int> &ps){
    
    M = p_M;
    
    //ps is the planted solution;
    
    vector<int> v;
    vector <bool> v_J;
    int i,j,k;
    int a = 0;
    bool J1, J2, J3;
    int flag;
    
    while (a<M){
        i=rand() % N ;
        j=rand() % N ;
        k=rand() % N ;
        
        if (i!=j && i!=k && j!=k){
            v   = make_vector<int>() << i  << j  << k;
            v_J = make_vector<bool>() << 0 << 0 << 0;
            flag=addFactor(a, v_J, v);
            
            if (flag){
                a++;
                //ptr is the pointer to the last added factor
                //please remind that F is a vector of pointers to Factor's.
                //in order to use this pointer as a pointer to XorsatFactor we need to use the cast operator
                Factor* ptr = F.back();
                F[ptr->f]->plantedClause(ps[i],ps[j],ps[k]);
                
                //J1 = ((FactorSat*)F[ptr->f])->v_J[0];
                //J2 = ((FactorSat*)F[ptr->f])->v_J[1];
                //J3 = ((FactorSat*)F[ptr->f])->v_J[2];
                //cout << ps[i] << " " << ps[j] << " " << ps[k]  << " " << ((2*J1 - 1) * (J1 - ps[i]) | (2*J2 - 1) * (J2 - ps[j]) | (2*J3 - 1) * (J3 - ps[k])) << endl;
               
            }
        }
    }
};



//input variables:
//ps  : planted configuration

bool Graph::check(vector<int> &ps){
    
    bool prod=1;
    int i, j, k;
    bool J1, J2, J3;
    for (int a=0; a<M; a++){
        i = F[a]->v_node[0];
        j = F[a]->v_node[1];
        k = F[a]->v_node[2];
        
        J1 = F[a]->v_J[0];
        J2 = F[a]->v_J[1];
        J3 = F[a]->v_J[2];
        
        prod*=((2*J1 - 1) * (ps[i] - (1-J1)) | (2*J2 - 1) * (ps[j] - (1-J2)) | (2*J3 - 1) * (ps[k] - (1-J3)));
        
    }
    
    return prod;
}




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Messages
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



Messages::Messages(Graph& p_G) : G (p_G) {
   
    N=G.numberOfTotalNodes();
    M=G.numberOfTotalFactors();
    
    xi_NodeToFac.resize(N);
    Hxi_FacToNode.resize(M);
    p_bias.resize(N);
    marginal.resize(N);
    
    initMessages();
    
};



void Messages::initMessages(){
    
    //Init messages xi_NodeToFac to 0.5
    //which basically means to init messages eta_NodeToFac to uniform distributions
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        int size_di = G.v[it_i->n].numberOfFactors();  //di is the set of factors attached to node i. size_di is the number of such factors
        xi_NodeToFac[it_i->n].resize(size_di);
        
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            int index_a = distance (G.v[it_i->n].v_fac.begin(), it_a);
            
            xi_NodeToFac[it_i->n][index_a] = 0.5;
            
        }
    }
    
    //resize Hxi_FacToNode
    for(vector<Factor*>::iterator it_a=G.F.begin(); it_a !=G.F.end(); ++it_a){
        int size_da = G.F[(*it_a)->f]->numberOfNodes(); //da is the set of nodes attached to factor a. size_da is the number of such nodes and should be p
        Hxi_FacToNode[(*it_a)->f].resize(size_da);
    }
    
    //set no bias, i.e set bias to the uniform distribution
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        p_bias[it_i->n] = 0.5;
    }
    
    //resize marginals
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        marginal[it_i->n].resize(2);
    }
    
};



void Messages::HxiUpdate(){
    
    double prod_xi;
    
    for (vector<Factor*>::iterator it_b = G.F.begin() ; it_b != G.F.end(); ++it_b){
        for (vector<int>::iterator it_i = G.F[(*it_b)->f]->v_node.begin() ; it_i != G.F[(*it_b)->f]->v_node.end(); ++it_i){
            //cout << "Hxi from factor " << (*it_b)->f << " to node " << *it_i << endl;
            //compute message nu_b_to_i
            int index_i = distance (G.F[(*it_b)->f]->v_node.begin(), it_i);
            vector <vector <double> > in_eta;
            
            prod_xi = 1.;
            
            for(vector<int>::iterator it_j = G.F[(*it_b)->f]->v_node.begin(); it_j != G.F[(*it_b)->f]->v_node.end(); ++it_j){
                if(*it_j!=*it_i){
                    //cout << "j: " << *it_j << endl;
                    vector<int>::iterator it = find(G.v[*it_j].v_fac.begin(), G.v[*it_j].v_fac.end(), (*it_b)->f);
                    int index_b = distance (G.v[*it_j].v_fac.begin(), it);
                    
                    prod_xi *= xi_NodeToFac[*it_j][index_b];
                    
                }
            }
            
            //Hxi_b_to_i
            Hxi_FacToNode[(*it_b)->f][index_i] = ( 1 - prod_xi ) / ( 2 - prod_xi );
            
        }
    }
    
};



void Messages::xiUpdate (){
    
    double prod_HxiS_1, prod_HxiS_2;
    double prod_HxiU_1, prod_HxiU_2;
    
    int index_i;
    bool Jai, Jbi;
    
    for (vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            //cout << "xi from node " << it_i->n << " to factor " << *it_a << endl;
            //compute message eta_i_to_a
            int index_a = distance (G.v[it_i->n].v_fac.begin(), it_a);
            
            prod_HxiS_1 = 1;
            prod_HxiU_1 = 1;
            
            prod_HxiS_2 = 1;
            prod_HxiU_2 = 1;

            
            vector<int>::iterator it = find(G.F[*it_a]->v_node.begin(), G.F[*it_a]->v_node.end(), it_i->n);
            index_i = distance (G.F[*it_a]->v_node.begin(), it);
            
            Jai = G.F[*it_a]->v_J[index_i];
            
            for(vector<int>::iterator it_b = G.v[it_i->n].v_fac.begin(); it_b != G.v[it_i->n].v_fac.end(); ++it_b){
                if(*it_b!=*it_a){
                    //cout << "b: " << *it_b << endl;
                    vector<int>::iterator it = find(G.F[*it_b]->v_node.begin(), G.F[*it_b]->v_node.end(), it_i->n);
                    index_i = distance (G.F[*it_b]->v_node.begin(), it);
                    
                    Jbi = G.F[*it_b]->v_J[index_i];
                    
                    if (Jbi == Jai){
                        prod_HxiS_1 *= Hxi_FacToNode[*it_b][index_i];
                        prod_HxiS_2 *= (1 - Hxi_FacToNode[*it_b][index_i]);
                    }
                    else{
                        prod_HxiU_1 *= (1 - Hxi_FacToNode[*it_b][index_i]);
                        prod_HxiU_2 *= Hxi_FacToNode[*it_b][index_i];
                    }
                    
                }
                
            }
            
            //xi_i_to_a
            if (Jai == 1)
                xi_NodeToFac[it_i->n][index_a] = p_bias[it_i->n]*prod_HxiS_1*prod_HxiU_1 / ( p_bias[it_i->n]*prod_HxiS_1*prod_HxiU_1 + (1-p_bias[it_i->n])*prod_HxiS_2*prod_HxiU_2 );
            else
                xi_NodeToFac[it_i->n][index_a] = (1-p_bias[it_i->n])*prod_HxiS_1*prod_HxiU_1 / ( (1-p_bias[it_i->n])*prod_HxiS_1*prod_HxiU_1 + p_bias[it_i->n]*prod_HxiS_2*prod_HxiU_2 );
            
        }
    }
    
};



void Messages::nodeMarginals(int verbose=0){
    
    double prod_HxiS_1, prod_HxiS_2;
    double prod_HxiU_1, prod_HxiU_2;
    
    bool Jai;
    
    //for each node i
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        vector <double> marg(2,0.);
        
        prod_HxiS_1 = 1;
        prod_HxiU_1 = 1;

        prod_HxiS_2 = 1;
        prod_HxiU_2 = 1;

        //we compute the product of the messages nu_FacToNode[a][i] coming from the factor attached to i
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            vector<int>::iterator it = find(G.F[*it_a]->v_node.begin(), G.F[*it_a]->v_node.end(), it_i->n);
            int index_i = distance (G.F[*it_a]->v_node.begin(), it);
            
            Jai = G.F[*it_a]->v_J[index_i];
            
            if (Jai == 1){
                prod_HxiS_1 *= Hxi_FacToNode[*it_a][index_i];
                prod_HxiS_2 *= (1 - Hxi_FacToNode[*it_a][index_i]);
            }
            else{
                prod_HxiU_1 *= (1 - Hxi_FacToNode[*it_a][index_i]);
                prod_HxiU_2 *= Hxi_FacToNode[*it_a][index_i];
            }
            
        }
        
        marg[1] = prod_HxiS_1 * prod_HxiU_1;
        marg[0] = prod_HxiS_2 * prod_HxiU_2;
        
        //and we take into account the bias on the node i
        marg[1] *= p_bias[it_i->n];
        marg[0] *= (1-p_bias[it_i->n]);
        
        //finally we normalize the marginal
        double sum=0.;
        for (int k=0; k<2; k++){
            sum+=marg[k];
        }
        
        if (!sum){
            if(verbose)
                cout << "node " << it_i->n << " receives conflicting messages, see below: " << endl;
        }
        else{
            for (int k=0; k<2; k++){
                marg[k]/=sum;
            }
            marginal[it_i->n]=marg;
        }
        
    }
    

};
 
 

//the three following methods of the class Messages print the state of the BP algorithm
void Messages::nuState(){
    
    bool Jbi;
    
    cout << endl;
    cout << "---*---*---*---printing messages from factors to nodes---*---*---*---" << endl;
    
    for (vector<Factor*>::iterator it_b = G.F.begin() ; it_b != G.F.end(); ++it_b){
        for (vector<int>::iterator it_i = G.F[(*it_b)->f]->v_node.begin() ; it_i != G.F[(*it_b)->f]->v_node.end(); ++it_i){
            int index_i = distance (G.F[(*it_b)->f]->v_node.begin(), it_i);
            
            Jbi = G.F[(*it_b)->f]->v_J[index_i];

            cout << "message from factor " << (*it_b)->f << " to node " << *it_i << " Jbi= " << Jbi << endl;

            if (Jbi == 1)
                cout << 1-Hxi_FacToNode[(*it_b)->f][index_i] << " " <<  Hxi_FacToNode[(*it_b)->f][index_i]  << endl;
            else
                cout << Hxi_FacToNode[(*it_b)->f][index_i]   << " " << 1-Hxi_FacToNode[(*it_b)->f][index_i] << endl;
            
        }
    }
    
};



void Messages::etaState(){

    bool Jai;
    
    cout << endl;
    cout << "---*---*---*---printing messages from nodes to factors---*---*---*---" << endl;

    for (vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            cout << "message from node " << it_i->n << " to factor " << *it_a << endl;
            int index_a = distance (G.v[it_i->n].v_fac.begin(), it_a);
            
            vector<int>::iterator it = find(G.F[*it_a]->v_node.begin(), G.F[*it_a]->v_node.end(), it_i->n);
            int index_i = distance (G.F[*it_a]->v_node.begin(), it);
            
            Jai = G.F[*it_a]->v_J[index_i];
        
            if (Jai == 1)
                cout << 1-xi_NodeToFac[it_i->n][index_a] << " " <<  xi_NodeToFac[it_i->n][index_a]  << endl;
            else
                cout << xi_NodeToFac[it_i->n][index_a]   << " " << 1-xi_NodeToFac[it_i->n][index_a] << endl;
                
    }
  }

};



void Messages::marginalState(){

  cout << endl;
  cout << "---*---*---*---printing marginals---*---*---*---" << endl;

  for (vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
    cout << "marginal of node " << it_i->n << endl;
    cout << marginal[it_i->n][0] << " " << marginal[it_i->n][1] << endl;
  }

};




/*



 
 void Messages::setHardBias(vector<int>& v_bias, vector<int>& v_q){
 
 //v_bias contains the indices of nodes that have to be biased
 //v_q contains the color towards which they have to be biased
 
 int size=v_bias.size();
 if (size != v_q.size()){
 cout << "error: v_q and v_bias have to have the same size!" << endl;
 return;
 }
 else{
 for (int i=0; i<size; i++){
 int n=v_bias[i];
 vector <double> b(q,0);
 //we set the bias towards the index color specified in v_q[i]
 b[v_q[i]]=1;
 bias[n]=b;
 }
 }
 
 };
 */


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class BP
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



BP::BP(Graph& p_G, bool p_verbose=0) : G (p_G), mess (G), verbose(p_verbose) { N=G.N; };



void BP::BPsweep(){
    
    mess.HxiUpdate();
    mess.xiUpdate();
    
    mess.nodeMarginals(verbose);
    
    if(verbose)
        BPprint();
    
};



void BP::BPprint(){
    
    mess.nuState();
    mess.etaState();
    mess.marginalState();
    
};



bool BP::BPiteration(double eps, int T){
    
    int    t = 0;
    int    f = 1;
    double tmp_eps = 1;
    
    initPreviousMarginals();
    
    while (tmp_eps > eps && t < T){
        BPsweep();
        
        if (t == 2) {
            f = findFrustratedSpins();
            if (f == 0){
                cout << "****************** ERROR: ************************************* the variable set by the last calling of warningPropagation causes contraddictions" << endl;
                break;
            }
        }
        
        tmp_eps = compareMarginals();
        storePreviousMarginals();
        if(verbose){
            cout << endl;
            cout << "BP iteration: at time t=" << t << " the maximum error between current and previous marginals is " << tmp_eps << endl;
        }
        
        t++;
    }
    
    return f;
    
}



void BP::initPreviousMarginals(){
    
    prev_marginal.resize(N);
    
    for(int i = 0; i < N; i++){
        prev_marginal[i] = 0.;
    }
    
};



void BP::storePreviousMarginals(){
    
    for(int i = 0; i < N; i++){
        prev_marginal[i] = mess.marginal[i][0];
    }
    
};



double BP::compareMarginals(){
    
    double tmp, max = 0.;
    
    for(int i = 0; i < N; i++){
        tmp = abs(prev_marginal[i] - mess.marginal[i][0]);
        if (tmp > max)
            max = tmp;
    }
    
    return max;
    
};



bool BP::findFrustratedSpins(){
    
    vector<int>().swap(frustratedSpins);
    
    bool f=1;
    
    
    double prod_HxiS_1, prod_HxiS_2;
    double prod_HxiU_1, prod_HxiU_2;
    
    bool Jai;
    
    //for each node i
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        vector <double> marg(2,0.);
        
        prod_HxiS_1 = 1;
        prod_HxiU_1 = 1;
        
        prod_HxiS_2 = 1;
        prod_HxiU_2 = 1;
        
        //we compute the product of the messages nu_FacToNode[a][i] coming from the factor attached to i
        for (vector<int>::iterator it_a = G.v[it_i->n].v_fac.begin(); it_a != G.v[it_i->n].v_fac.end(); ++it_a){
            vector<int>::iterator it = find(G.F[*it_a]->v_node.begin(), G.F[*it_a]->v_node.end(), it_i->n);
            int index_i = distance (G.F[*it_a]->v_node.begin(), it);
            
            Jai = G.F[*it_a]->v_J[index_i];
            
            if (Jai == 1){
                prod_HxiS_1 *= mess.Hxi_FacToNode[*it_a][index_i];
                prod_HxiS_2 *= (1 - mess.Hxi_FacToNode[*it_a][index_i]);
            }
            else{
                prod_HxiU_1 *= (1 - mess.Hxi_FacToNode[*it_a][index_i]);
                prod_HxiU_2 *= mess.Hxi_FacToNode[*it_a][index_i];
            }
            
        }
        
        marg[1] = prod_HxiS_1;
        marg[0] = prod_HxiS_2;
        
        
        //and we take into account the bias on the node i
        marg[0] *= (1-mess.p_bias[it_i->n]);
        marg[1] *= mess.p_bias[it_i->n];
        
        
        //we check if marginal is normalizable. if it is not, the messages is receiving conflicting messages
        double sum=0.;
        for (int k=0; k<2; k++){
            sum+=marg[k];
        }
        if (!sum){
            frustratedSpins.push_back(it_i->n);
            f=0;
        }
        
    }
    
    return f;
    
};



/*
void BP::initDecimation(vector<int>& v_bias, vector<int>& v_q){
    
    for(int i=0; i<N; i++)
        notFixedSpins.push_back(i);
    
    if(v_bias.size()){

        int size=v_bias.size();
        if (size != v_q.size()){
            cout << "error: v_q and v_bias have to have the same size!" << endl;
            return;
        }
        else{
            for (int i=0; i<size; i++){
                int n=v_bias[i];
                fixedSpins.push_back(n);
                fixedValues.push_back(v_q[i]);
                notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), n), notFixedSpins.end());
                
                //we bias eta's towards the index color specified in v_q[i]
 
                //for (vector<int>::iterator it_a = G.v[n].v_fac.begin(); it_a != G.v[n].v_fac.end(); ++it_a){
                //    int index_a = distance (G.v[n].v_fac.begin(), it_a);
                //    vector <double> eta(q,0);
                //    eta[v_q[i]]=1;
                //    mess.eta_NodeToFac[n][index_a]=eta;
                //}
                
            }
        }

        mess.setHardBias(v_bias, v_q);
        
    }

};

 */









/*
bool BP::warningDecimation(int verbose=0){
    
    //Pictorially this is how this method works:
    
    //time t-1:
    //set bias and fix a group G of variables
    //time t:
    //update nu  (using eta at time t-1)
    //update eta (using nu  at time t and bias at time t-1)
    //at this time step, no variables get fixed because of the set G, but some other variables
    //may get fixed because of the variables that we fixed at time t-2.
    //time t+1:
    //update nu  (using eta at time t)
    //use these nu's to set bias and fix other variables
    
    //it is clear then that the variables fixed at time t-1 have some effect on the others after 2 iterations, not 1.
    
    //g is the number of spins that became frozen (decimated) at each iteration.
    //g_past is the number of spins that became frozen (decimated) at the previous iteration.
    //we do keep track of both for the reason explained above.
    //f is a boolean variable that is set to 0 when frustrated variables are found
    
    int  g=1, g_past=1;
    int  t=0;
    bool f1=1;
    bool f2=1;
    
    while (f1==1 && (g+g_past>0)){
    
        f1 = BPsweep(verbose);
        
        f2 = findFrustratedSpins();
        
        //f1 and f2 have to be equal
        assert(f1==f2);

        if (f1 == 0){
            cout << "****************** ERROR: ************************************* frustrated variables are found in warningPropagation" << endl;
            return f1;
            break;
        }

        
        g_past = g;
        g=fixSpins(verbose);

        t++;
    
    }
    
    return f1;
    
}


void BP::BPguidedDecimation(int TT, int verbose=0){
    
    initPreviousMarginals();
    
    double eps = 0.001;
    double f;
    int    T = 100;
    int    t=0;
    
    while (t < TT && notFixedSpins.size()>0){

    
        if(verbose){
            cout << "--------------------------------------time: " << t << "------------------------------------ " << endl;
            cout << endl;
            cout << "frozen variables:" << " (size=" << fixedSpins.size()    << ")" << endl;
            vec_print(fixedSpins);
            cout << "free variables:"   << " (size=" << notFixedSpins.size() << ")" << endl;
            vec_print(notFixedSpins);
            cout << endl;
            cout << "run BP until convergence: " << endl;
        }
        
        
        f = BPiteration(eps, T, verbose);
        
        if (f == 0)
            break;
        
        vector<int> v_bias;
        vector<int> v_q;
    
        findMostBiased(v_bias, v_q);
        
        if(v_bias.size() != 0)
            mess.setHardBias(v_bias,v_q);
    
        if (verbose){
            cout << endl;
            cout << "*************** printing the nodes that gets frozen at this time step: " << endl;
            cout << "*************** we freeze the most biased node: " << v_bias[0] << " towards the color " << v_q[0] <<  endl;
        }
        
        fixedSpins.push_back(v_bias[0]);
        fixedValues.push_back(v_q[0]);
        notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), v_bias[0]), notFixedSpins.end());
        
        if(verbose)
            cout << "--------------------------------------------------------------------------------------------------------------------------- call warningDecimation" << endl;
    
        f = warningDecimation(verbose);
        
        if(verbose){
            cout << "--------------------------------------------------------------------------------------------------------------------------- warning propagated" << endl;
            if(f)
                cout << "--------------------------------------------------------------------------------------------------------------------------- no contraddictions found" << endl;
            else
                cout << "--------------------------------------------------------------------------------------------------------------------------- contraddictions found" << endl;
            cout << endl;
        }
        
        if (f == 0){
            break;
        }
        t++;
        
    }
    
    if (verbose){
        cout << endl;
        cout << "----------------------------------END OF ITERATION----------------------------------------" << endl;
        cout << "(maybe partial) solution: " << endl;
        cout << "size of the solution: "     << fixedSpins.size() << endl;
        vec_print(fixedSpins);
        vec_print(fixedValues);
    }


};
 
*/




/*
void BP::findMostBiased(vector<int> & v_bias, vector<int>& v_q){
    
    int    i_max = notFixedSpins[0];
    int    col = 0;
    double tmp, max = 0.;
    
    for (vector<int>::iterator it_i = notFixedSpins.begin(); it_i != notFixedSpins.end(); ++it_i){
        tmp = mess.marginal[*it_i][0]-mess.marginal[*it_i][1];
        if (abs(tmp) > max){
            i_max = *it_i;
            
            if (tmp > 0) col = 0;
            else col = 1;
            
            max = abs(tmp);
        }
        
    }
    
    v_bias.push_back(i_max);
    v_q.push_back(col);
    
};

void BP::fixRandomSpin(int verbose){
    
    vector<int> v_bias;
    vector<int> v_q;
    bool col;
    
    int N_free = notFixedSpins.size();
    int i = rand() % N_free ;
    if (2*(double)rand()/RAND_MAX-1 > 0)
        col = 1;
    else
        col = 0;
    
    int n = notFixedSpins[i];
    
    v_bias.push_back(n);
    v_q.push_back(col);
    fixedSpins.push_back(n);
    fixedValues.push_back(col);
    notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), n), notFixedSpins.end());
    
    if(v_bias.size() != 0)
        mess.setHardBias(v_bias,v_q);
    
    if (verbose){
        cout << "*************** printing the nodes that gets randomly frozen: " << endl;
        for(int i = 0; i < v_bias.size(); i++)
            cout << v_bias[i] << " " << v_q[i] << endl;
    }
    
};

*/



    
    
/*

int BP::fixSpins(int verbose){

    //this is function needs to be generalized to deal with a general q
    
    //v_bias contains the indices of nodes that have to be biased
    //v_q contains the color towards which they have to be biased

    vector<int> v_bias;
    vector<int> v_q;


    //for each factor a
    for (vector<Factor*>::iterator it_a = G.F.begin() ; it_a != G.F.end(); ++it_a){
        //we look at all the nodes to which a is sending messages
        for (vector<int>::iterator it_i = G.F[(*it_a)->f]->v_node.begin() ; it_i != G.F[(*it_a)->f]->v_node.end(); ++it_i){
            //if this not is not fustrated (i.e. it is not receiving conflicting messages from factors), we fix it
            bool flag = find(frustratedSpins.begin(), frustratedSpins.end(), *it_i) != frustratedSpins.end();
            if(!flag){
                
                int index_i = distance (G.F[(*it_a)->f]->v_node.begin(), it_i);
                //if this message contains a clear indication of the color towards which
                //the node should be biased, we add the node to v_bias and the color to v_q,
                //after having checked that the node has not been frozen yet.
                if (mess.nu_FacToNode[(*it_a)->f][index_i][0]==0.){
                    bool contains = find(fixedSpins.begin(), fixedSpins.end(), *it_i) != fixedSpins.end();
                    if(!contains){
                        v_bias.push_back(*it_i);
                        v_q.push_back(1);
                        fixedSpins.push_back(*it_i);
                        fixedValues.push_back(1);
                        notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), *it_i), notFixedSpins.end());

                    }
                    continue;
                }
                if (mess.nu_FacToNode[(*it_a)->f][index_i][1]==0.){
                    bool contains = find(fixedSpins.begin(), fixedSpins.end(), *it_i) != fixedSpins.end();
                    if(!contains){
                        v_bias.push_back(*it_i);
                        v_q.push_back(0);
                        fixedSpins.push_back(*it_i);
                        fixedValues.push_back(0);
                        notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), *it_i), notFixedSpins.end());
                    }
                    continue;
                }
            }
        }
    }
    
    if(v_bias.size() != 0)
        mess.setHardBias(v_bias,v_q);
    

    if (verbose){
        cout << "*************** printing the nodes that gets frozen at this time step: " << endl;
        for(int i=0; i<v_bias.size(); i++)
            cout << v_bias[i] << " " << v_q[i] << endl;

        cout << endl;
        cout << endl;
    }
    
    return v_bias.size();

};

*/

void vec_print(vector<int>& vec){
    for (int i=0; i<vec.size(); i++)
        cout << vec[i] << ' ';
    cout << endl;
}

