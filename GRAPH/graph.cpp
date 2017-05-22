#include "graph.h"



//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Node
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



int Node::numberOfFactors(){
    
  return v_factors.size();
    
};



//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Factor
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



//input variables:
//p_f : index of the factor
//p_p : # of the variables entering in a factor

Factor::Factor(int p_f, vector<bool> p_v_J){
    
    f   = p_f;
    v_J = p_v_J;
    p   = v_J.size();

    v_nodes.reserve(p);
    
};



int Factor::numberOfNodes(){
    
    return v_nodes.size();
    
};



bool Factor::clause(){
    
    if(v_values == v_J)
        return 0;
    else
        return 1;
    
};



//given a triplet of variables entering in a clause, we can construct several planted "clauses".
//a clause of the SAT problem is specified by three values J1, J2, J3, being 0 or 1 if the corresponding variables appear unnegated or negated in the clause,
//drawn at random in the non-planted setting.
//for instance, consider the clause x1 | !x2 | !x3 where ! indicates the negation. in this case we set J1=1, J2=0, J3=0.
//it is easy to see that the only choice of J's for which the triplet do not satisfy the clause is 0, 1, 1, and more generally, for a particular realization of (x1,x2,x3) we need
//to exclude J1=x1, J2=x2, J3=x3.

void Factor::plantedClause(){
    
    int  flag=1;
    
    while (flag){

        vector<bool> temp_v_J;
        
        for (int k = 0; k < p ; ++k){
        
            if (2*(double)rand()/RAND_MAX-1 > 0)
                temp_v_J.push_back(1);
            else
                temp_v_J.push_back(0);
    
        }
        
        if (v_values != temp_v_J){
            v_J = temp_v_J;
            flag=0;
        }
        
    }
    
};



//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Graph
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



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
  for (vector<int>::iterator it = v[i].v_factors.begin() ; it != v[i].v_factors.end(); ++it)
    cout << *it << endl;
};



void Graph::nodesOfFactor(int a){
  cout << "factor " << a << " has " << F[a].numberOfNodes() << " nodes: " << endl;
  for (vector<int>::iterator it = F[a].v_nodes.begin() ; it != F[a].v_nodes.end(); ++it)
    cout << *it << endl;
};



//input variables:
//p_a    : factor index;
//p_v_J  : p-component vector whose components are 1 or 0 depending if the corresponding variable is negated or not
//v_da   : nodes attached to the factor

int Graph::addFactor(int p_a, vector<bool> p_v_J, vector<int> v_da){
    
    vector <int> v1 = v_da;
    std::sort(v1.begin(), v1.end());
    
    int flag = 1;
    
    //before adding a factor we check that the factor does not already exist.
    //if the factor already exists, flag is set to 0.
    //to this aim it is sufficient to check that the first node of v_da does not appear in a factor with the p-1 other nodes of v_da.
    
    for(vector<int>::iterator it_a = v[v_da[0]].v_factors.begin(); it_a != v[v_da[0]].v_factors.end(); it_a++){
        
        vector <int> v2 = F[*it_a].v_nodes;
        std::sort(v2.begin(), v2.end());
        
        if (v1 == v2)
            flag=0;
        
    }
    
    if(flag){
        
        Factor a(p_a,p_v_J);
        
        for (int i = 0; i < v_da.size(); ++i)
            a.v_nodes.push_back(v_da[i]);
        
        F.push_back(a);
        
        for(int i = 0; i < v_da.size(); ++i)
            v[v_da[i]].v_factors.push_back(p_a);
    
    }
    
    return flag;
    
};



//input variables:
//p_M: # of factors

void Graph::ErdosRenyi(int p_M){
    
    M = p_M;
    
    
    int  a = 0;
    int  flag;
    
    while (a < M){
        
        vector<int>  v_f;               //vector of nodes entering in a factor
        vector<bool> v_J;               //vector of J values attached to nodes entering in a factor
        
        for (int i = 0; i < p; ++i)
            v_f.push_back(rand() % N );
        
        //check for no repetitions
        
        sort (v_f.begin(), v_f.end());
        vector<int>::iterator it = unique (v_f.begin(), v_f.end());
        v_f.resize(distance(v_f.begin(),it) );
        
        if (v_f.size() == p){
            
            for (int i = 0; i < p; ++i){
                if (2*(double)rand()/RAND_MAX-1 > 0)
                    v_J.push_back(1);
                else
                    v_J.push_back(0);
            }
        
            flag = addFactor(a, v_J, v_f);
            if (flag) a++;
            
        }
    
    }
    
};



//input variables:
//p_M : # of factors
//ps  : planted configuration

void Graph::plantedErdosRenyi(int p_M, vector<bool> &ps){
    
    M = p_M;
    
    for (int i = 0; i < N; ++i)
        v[i].value = ps[i];
    
    
    int  a = 0;
    int  flag;
    
    while (a < M){
        
        vector<int>  v_f;           //vector of nodes entering in a factor
        vector<bool> v_J;           //vector of J values attached to nodes entering in a factor
        vector<bool> val;           //values of the binary variables entering in a factor

        for (int i = 0; i < p; ++i){
            int n  = rand() % N;
            val.push_back(ps[n]);
            v_f.push_back(n);
        }
        
        //check for no repetitions
        
        sort (v_f.begin(), v_f.end());
        vector<int>::iterator it = unique (v_f.begin(), v_f.end());
        v_f.resize(distance(v_f.begin(),it) );
        
        if (v_f.size() == p){
        
            for (int i = 0; i < p; ++i)
                v_J.push_back(0);
            
            flag = addFactor(a, v_J, v_f);
            
            if (flag){
                
                a++;
                //ptr is a reference to the last added factor
                Factor ptr = F.back();
                F[ptr.f].v_values = val;
                F[ptr.f].plantedClause();
                
            }
        }
    }
};



void Graph::graphStructure(){
    
    cout << endl;
    cout << endl;
    
    cout << "structure of the graph: " << endl;
    cout << endl;
    
    for (int i = 0; i < N; i++)
        factorsOfNode(i);
    
    for (int a = 0; a < M; a++)
        nodesOfFactor(a);
    
    cout << endl;
    cout << endl;
    
    cout << "factors types: " << endl;
    for (int a = 0; a < M; a++){
        cout << a << " J's: ";
        for (int k = 0; k < p; ++k)
            cout << F[a].v_J[k] << " ";
        cout << endl;
    }
    
    
    cout << endl;
    cout << endl;
    
};



//input variables:
//ps  : configuration

bool Graph::check(vector<bool> &ps){
    
    bool prod = 1;
    
    for (int a = 0; a < M; ++a){
        
        if (F[a].v_values == F[a].v_J)
            prod = 0;
        
    }
    
    return prod;
}




//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Messages
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



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
        
        for (vector<int>::iterator it_a = G.v[it_i->n].v_factors.begin(); it_a != G.v[it_i->n].v_factors.end(); ++it_a){
            int index_a = distance (G.v[it_i->n].v_factors.begin(), it_a);
            
            xi_NodeToFac[it_i->n][index_a] = 0.5;
            
        }
    }
    
    //resize Hxi_FacToNode
    for(vector<Factor>::iterator it_a=G.F.begin(); it_a !=G.F.end(); ++it_a){
        int size_da = G.F[it_a->f].numberOfNodes(); //da is the set of nodes attached to factor a. size_da is the number of such nodes and should be p
        Hxi_FacToNode[it_a->f].resize(size_da);
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
    
    for (vector<Factor>::iterator it_b = G.F.begin() ; it_b != G.F.end(); ++it_b){
        for (vector<int>::iterator it_i = G.F[it_b->f].v_nodes.begin() ; it_i != G.F[it_b->f].v_nodes.end(); ++it_i){
            //cout << "Hxi from factor " << (*it_b).f << " to node " << *it_i << endl;
            //compute message nu_b_to_i
            int index_i = distance (G.F[it_b->f].v_nodes.begin(), it_i);
            vector <vector <double> > in_eta;
            
            prod_xi = 1.;
            
            for(vector<int>::iterator it_j = G.F[it_b->f].v_nodes.begin(); it_j != G.F[it_b->f].v_nodes.end(); ++it_j){
                if(*it_j != *it_i){
                    //cout << "j: " << *it_j << endl;
                    vector<int>::iterator it = find(G.v[*it_j].v_factors.begin(), G.v[*it_j].v_factors.end(), it_b->f);
                    int index_b = distance (G.v[*it_j].v_factors.begin(), it);
                    
                    prod_xi *= xi_NodeToFac[*it_j][index_b];
                    
                }
            }
            
            //Hxi_b_to_i
            Hxi_FacToNode[it_b->f][index_i] = ( 1 - prod_xi ) / ( 2 - prod_xi );
            
        }
    }
    
};



void Messages::xiUpdate (){
    
    double prod_HxiS_1, prod_HxiS_2;
    double prod_HxiU_1, prod_HxiU_2;
    
    int index_i;
    bool Jai, Jbi;
    
    for (vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
        for (vector<int>::iterator it_a = G.v[it_i->n].v_factors.begin(); it_a != G.v[it_i->n].v_factors.end(); ++it_a){
            //cout << "xi from node " << it_i->n << " to factor " << *it_a << endl;
            //compute message eta_i_to_a
            int index_a = distance (G.v[it_i->n].v_factors.begin(), it_a);
            
            prod_HxiS_1 = 1;
            prod_HxiU_1 = 1;
            
            prod_HxiS_2 = 1;
            prod_HxiU_2 = 1;

            
            vector<int>::iterator it = find(G.F[*it_a].v_nodes.begin(), G.F[*it_a].v_nodes.end(), it_i->n);
            index_i = distance (G.F[*it_a].v_nodes.begin(), it);
            
            Jai = G.F[*it_a].v_J[index_i];
            
            for(vector<int>::iterator it_b = G.v[it_i->n].v_factors.begin(); it_b != G.v[it_i->n].v_factors.end(); ++it_b){
                if(*it_b!=*it_a){
                    //cout << "b: " << *it_b << endl;
                    vector<int>::iterator it = find(G.F[*it_b].v_nodes.begin(), G.F[*it_b].v_nodes.end(), it_i->n);
                    index_i = distance (G.F[*it_b].v_nodes.begin(), it);
                    
                    Jbi = G.F[*it_b].v_J[index_i];
                    
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



bool Messages::nodeMarginals(int verbose=0){
    
    double prod_HxiS_1, prod_HxiS_2;
    double prod_HxiU_1, prod_HxiU_2;
    
    bool Jai;
    
    bool f = 1;
    
    //for each node i
    for(vector<Node>::iterator it_i=G.v.begin(); it_i !=G.v.end(); ++it_i){
        vector <double> marg(2,0.);
        
        prod_HxiS_1 = 1;
        prod_HxiU_1 = 1;

        prod_HxiS_2 = 1;
        prod_HxiU_2 = 1;

        //we compute the product of the messages nu_FacToNode[a][i] coming from the factor attached to i
        for (vector<int>::iterator it_a = G.v[it_i->n].v_factors.begin(); it_a != G.v[it_i->n].v_factors.end(); ++it_a){
            vector<int>::iterator it = find(G.F[*it_a].v_nodes.begin(), G.F[*it_a].v_nodes.end(), it_i->n);
            int index_i = distance (G.F[*it_a].v_nodes.begin(), it);
            
            Jai = G.F[*it_a].v_J[index_i];
            
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
        for (int k = 0; k < 2; ++k){
            sum += marg[k];
        }
        
        if (!sum){
            f = 0;
            if(verbose)
                cout << "node " << it_i->n << " receives conflicting messages, see below: " << endl;
        }
        else{
            for (int k = 0; k < 2; ++k){
                marg[k] /= sum;
            }
            marginal[it_i->n] = marg;
        }
        
    }
    
    return f;

};
 
 

//the three following methods of the class Messages print the state of the BP algorithm
void Messages::nuState(){
    
    bool Jbi;
    
    cout << endl;
    cout << "---*---*---*---printing messages from factors to nodes---*---*---*---" << endl;
    
    for (vector<Factor>::iterator it_b = G.F.begin() ; it_b != G.F.end(); ++it_b){
        for (vector<int>::iterator it_i = G.F[it_b->f].v_nodes.begin() ; it_i != G.F[it_b->f].v_nodes.end(); ++it_i){
            int index_i = distance (G.F[it_b->f].v_nodes.begin(), it_i);
            
            Jbi = G.F[it_b->f].v_J[index_i];

            cout << "message from factor " << it_b->f << " to node " << *it_i << " Jbi= " << Jbi << endl;

            if (Jbi == 1)
                cout << 1-Hxi_FacToNode[it_b->f][index_i] << " " <<  Hxi_FacToNode[it_b->f][index_i]  << endl;
            else
                cout <<  Hxi_FacToNode[it_b->f][index_i]  << " " << 1-Hxi_FacToNode[it_b->f][index_i] << endl;
            
        }
    }
    
};



void Messages::etaState(){

    bool Jai;
    
    cout << endl;
    cout << "---*---*---*---printing messages from nodes to factors---*---*---*---" << endl;

    for (vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
        for (vector<int>::iterator it_a = G.v[it_i->n].v_factors.begin(); it_a != G.v[it_i->n].v_factors.end(); ++it_a){
            
            cout << "message from node " << it_i->n << " to factor " << *it_a << endl;
            int index_a = distance (G.v[it_i->n].v_factors.begin(), it_a);
            
            vector<int>::iterator it = find(G.F[*it_a].v_nodes.begin(), G.F[*it_a].v_nodes.end(), it_i->n);
            int index_i = distance (G.F[*it_a].v_nodes.begin(), it);
            
            Jai = G.F[*it_a].v_J[index_i];
        
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



//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class BP
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



BP::BP(Graph& p_G, bool p_verbose=0) : G (p_G), mess (G), verbose(p_verbose) { N=G.N; };



bool BP::BPsweep(){
    
    mess.HxiUpdate();
    mess.xiUpdate();
    
    bool f = mess.nodeMarginals(verbose);
    
    if(verbose)
        BPprint();
    
    return f;
    
};



void BP::BPprint(){
    
    mess.nuState();
    mess.etaState();
    mess.marginalState();
    
};



bool BP::BPiteration(double eps, int T){
    
    int    t = 0;
    double tmp_eps = 1;
    
    bool f = 1;
    
    initPreviousMarginals();
    
    while (tmp_eps > eps && t < T && f == 1){
        f = BPsweep();
        
        tmp_eps = compareMarginals();
        storePreviousMarginals();
        if(verbose){
            cout << endl;
            cout << "BP iteration: at time t=" << t << " the maximum error between current and previous marginals is " << tmp_eps << endl;
        }
        
        t++;
        
        if (f == 0)
            return 0;
    }
    
    if (tmp_eps > eps)
        f = 0;
    
    return f;
    
};



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



void BP::initDecimation(vector<int>& v_bias, vector<bool>& v_q){
    
    for(int i = 0; i < N; ++i)
        notFixedSpins.push_back(i);
    
    if(v_bias.size()){

        int size = v_bias.size();
        if (size != v_q.size()){
            cout << "error: v_q and v_bias have to have the same size!" << endl;
            return;
        }
        else{
            for (int i = 0; i < size; ++i){
                int n = v_bias[i];
                fixedSpins.push_back(n);
                fixedValues.push_back(v_q[i]);
                notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), n), notFixedSpins.end());
                
            }
        }

        setHardBias(v_bias, v_q);
        
    }

};



void BP::setHardBias(vector<int>& v_bias, vector<bool>& v_q){
    
    //v_bias contains the indices of nodes that have to be biased
    //v_q contains the color towards which they have to be biased
    
    int size = v_bias.size();
    if (size != v_q.size()){
        cout << "error: v_q and v_bias have to have the same size!" << endl;
        return;
    }
    else{
        for (int i = 0; i < size; ++i){
            int n = v_bias[i];
            //we set the bias towards the value specified in v_q[i]
            if (v_q[i] == 0)
                mess.p_bias[n] = 0.;
            else
                mess.p_bias[n] = 1.;
        }
    }
    
};



void BP::BPguidedDecimation(double eps, int T, int TT){
    
    
    bool f;
    int  t = 0;
    
    while (t < TT && notFixedSpins.size() > 0){

    
        if(verbose){
            cout << "-------------------------------external time: " << t << "------------------------------------ " << endl;
            cout << endl;
            cout << "frozen variables:" << " (size=" << fixedSpins.size()    << ")" << endl;
            vec_print(fixedSpins);
            cout << "free variables:"   << " (size=" << notFixedSpins.size() << ")" << endl;
            vec_print(notFixedSpins);
            cout << endl;
            cout << "run BP until convergence: " << endl;
        }
        
        f = BPiteration(eps, T);
        
        if (f == 0){
            cout << endl;
            cout << "************************************----- contraddictions found -----************************************" << endl;
            break;
        }
        
        vector<int> v_bias;
        vector<bool> v_q;
    
        findMostBiased(v_bias, v_q);
        
        if(v_bias.size() != 0)
            setHardBias(v_bias,v_q);
        
        if (verbose){
            cout << endl;
            cout << "*************** printing the nodes that gets frozen at time step: " << t << endl;
            cout << "*************** we freeze the most biased node: " << endl;
            for (int i = 0; i < v_bias.size(); ++i){
                cout << v_bias[i] << " towards the color " << v_q[i] <<  endl;
            }
        }
        
        for (int i = 0; i < v_bias.size(); ++i){
            fixedSpins.push_back(v_bias[i]);
            fixedValues.push_back(v_q[i]);
            notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), v_bias[i]), notFixedSpins.end());
        }
        
        
        cout << endl;
        cout << "fixed spins" << endl;
        vec_print(fixedSpins);
        cout << "fixed values" << endl;
        vec_print(fixedValues);
        cout << "bias" << endl;
        vec_print(mess.p_bias);
        cout << "number of free variables" << endl;
        cout << notFixedSpins.size() << endl;
        

        if(verbose)
            cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- call warningDecimation" << endl;
        
        int tw = 0;
        f = warningDecimation(tw);
        
         
        if(verbose){
            cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------- warning propagated for " << tw << " time steps" << endl;
            if(f)
                cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- no contraddictions found" << endl;
            else
                cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- contraddictions found" << endl;
            cout << endl;
        }
        
         
        if (f == 0){
            cout << endl;
            cout << "************************************----- contraddictions found -----************************************" << endl;
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


 
void BP::findMostBiased(vector<int> & v_bias, vector<bool>& v_q){
    
    int    i_max = notFixedSpins[0];
    bool   col;
    double tmp, max = 0.;
    
    
    for (vector<int>::iterator it_i = notFixedSpins.begin(); it_i != notFixedSpins.end(); ++it_i){
        tmp = mess.marginal[*it_i][0] - mess.marginal[*it_i][1];
        if (abs(tmp) > 0.999){
            
            if (tmp > 0)
                col = 0;
            else
                col = 1;
            
            v_bias.push_back(*it_i);
            v_q.push_back(col);
            
        }
        
    }
    
    if (v_bias.size() == 0){
        
        for (vector<int>::iterator it_i = notFixedSpins.begin(); it_i != notFixedSpins.end(); ++it_i){
            tmp = mess.marginal[*it_i][0] - mess.marginal[*it_i][1];
            if (abs(tmp) > max){
                
                if (tmp > 0)
                    col = 0;
                else
                    col = 1;
                
                max = abs(tmp);
                
                i_max = *it_i;
                
            }
            
        }
        
        v_bias.push_back(i_max);
        v_q.push_back(col);

    }
    
    
};



bool BP::warningDecimation(int &t){
    
    //Pictorially this is how this method works:
    
    //iteration time t-1:
    //set bias and fix a group G(t-1) of variables (generally one)
    //iteration time t:
    //update nu  (using eta at iteration time t-1)
    //update eta (using nu  at iteration time t and bias at iteration time t-1)
    //if we set a bias on variable j at iteration time t-1, at time t the eta's from j
    //to all the factors to which it is connected will feel this bias.
    //by the way, no extra variables get fixed because of the set G(t-1).
    //G(t) is made by variables that get fixed because of biases set at iteration time t-2.
    //iteration time t+1:
    //update nu  (using eta at iteration time t)
    //update eta (using nu  at iteration time t+1 and bias at iteration time t)
    
    //variables fixed at time t-1 have some effect on the others after 2 iterations, not 1.
    
    //storedMarginals1 and storedMarginals2 contain the marginals at time t and t+2.
    //at each time t of the iteration of the warnings, we fill these vectors with only the marginals of the variables that have not been fixed yet.
    
    vector <double> storedMarginals1;
    storedMarginals1.resize(N);
    vector <double> storedMarginals2;
    storedMarginals2.resize(N);
    
    bool flag = 1;
    bool f1, f2;
    
    for (vector<int>::iterator it_i = notFixedSpins.begin(); it_i != notFixedSpins.end(); ++it_i)
        storedMarginals1[*it_i] = mess.marginal[*it_i][0];

    f1 = BPsweep();
    
    if (f1 == 0)
        return 0;
    
    t = t + 1;
    
    while (flag == 1 && notFixedSpins.size() > 0){

        flag = 0;
        
        f2 = BPsweep();
    
        if (f2 == 0)
            return 0;
        
        for (vector<int>::iterator it_i = notFixedSpins.begin(); it_i != notFixedSpins.end(); ++it_i)
            storedMarginals2[*it_i] = mess.marginal[*it_i][0];
        
        vector <int> v_bias;
        vector <bool> v_q;

     
        for (vector<int>::iterator it_i = notFixedSpins.begin(); it_i != notFixedSpins.end(); ++it_i){
            
            if (storedMarginals2[*it_i] == 0 && storedMarginals1[*it_i] != 0){
                int i = notFixedSpins[*it_i];
                fixedSpins.push_back(i);
                fixedValues.push_back(1);
                notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), i), notFixedSpins.end());
                flag = 1;
                v_bias.push_back(i);
                v_q.push_back(1);
            }
            if (storedMarginals2[*it_i] == 1 && storedMarginals1[*it_i] != 1){
                int i = notFixedSpins[*it_i];
                fixedSpins.push_back(i);
                fixedValues.push_back(0);
                notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), i), notFixedSpins.end());
                flag = 1;
                v_bias.push_back(i);
                v_q.push_back(0);
            }
            
        }
        
        if (flag == 1)
            setHardBias(v_bias,v_q);
        
        cout << "------------------------------------------- updating during warning decimation -------------------------------------------- " << endl;
        cout << endl;
        cout << "fixed spins" << endl;
        vec_print(fixedSpins);
        cout << "fixed values" << endl;
        vec_print(fixedValues);
        cout << "bias" << endl;
        vec_print(mess.p_bias);
        cout << "number of free variables" << endl;
        cout << notFixedSpins.size() << endl;
        cout << endl;
        vec_print(storedMarginals1);
        vec_print(storedMarginals2);
        cout << endl;
        
        storedMarginals1 = storedMarginals2;
        
        t++;
    }

    return 1;
    
}