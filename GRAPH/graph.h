#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <cassert>

using namespace std;


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------// class representing a node
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Node{
public:

    Node(int p_n) : n ( p_n ) {};
  
    int n;                                                  //n is the index of the node
    
    bool value;                                             //value of the binary variable
    
    vector <int> v_factors;                                 //v_factors contains the indices of factors attached to n
    
    int numberOfFactors();                                  //this method returns the number of factors to which the node n is attached
    
};



//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------------------------------------------------------------------// class representing a factor
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Factor{
public:
    
    Factor(int p_f, vector<bool> p_v_J);
    
    int f;                                                  //f is the index of the factor
    int p;                                                  //p is the number of variables entering in a factor. p=3 in the 3-SAT problem

    vector <int> v_nodes;                                   //v_nodes contains the indices of nodes attached to f
                                                            //it is specified by the method addFactor of the class Graph
    
    vector <bool> v_values;                                 //v_values contains the values of the binary variables attached to the factor
    
    vector <bool> v_J;                                      //v_J is a set of p binary variables specifying the SAT clause.
                                                            //consider the clause x1 | !x2 | !x3 where ! indicates the negation.
                                                            //be Jk the k-component of the vector v_J.
                                                            //in this case we set J1=0, J2=1, J3=1.
                                                            //it is specified by the constructor

    bool clause();                                          //this method implements the SAT clause represented by a factor.
                                                            //it returns 1 if the clause is satified, and 0 if it is not.
    
    void plantedClause();                                   //this method sets the vector v_J according to v_nodes.
                                                            //given a set of p variables specified in v_nodes,
                                                            //there are several options to specify a v_J such that the clause is satisfied.
                                                            //the method picks one of them at random.
    
    int numberOfNodes();                                    //this method returns the numer of variables to which a factor is attached. this number is p by definition

};



//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------------------------------------------------------------------// class representing a factor graph
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Graph{
public:
    
    Graph(int, int);
    
    int N;
    int p;
    int M;
    
    vector <Node> v;                                        //v contains all the nodes of the graph
                                                            //this vector is filled by the constructor.
    
    vector <Factor> F;                                      //F contains pointers to all the factors of the graph.
                                                            //this vector is filled by the addFactor method.
    
    
    int numberOfTotalFactors();                             //this method returns the total number of factors in the graph
    int numberOfTotalNodes();                               //this method returns the total number or nodes in the graph

    void factorsOfNode(int);                                //this method returns the factors attached to the input node
    void nodesOfFactor(int);                                //this method returns the nodes attached to the input factor
    
    void ErdosRenyi(int);                                   //this method build an ER graph
                                                            //input variables:
                                                            //M : # of factors
    
    void plantedErdosRenyi(int, vector<bool> &);            //this method build an ER graph for which the second argument is a solution
                                                            //inputs:
                                                            //M  : # of factors
                                                            //ps : planted solution

    
    int addFactor(int, vector<bool>, vector<int>);          //the addFactor method implement the operations needed when adding a factor in the graph
                                                            //namely one needs to create a FactorSat object a, store its neighour variables and for each of them
                                                            //add the index of a as a neghbouring factor.
                                                            //input variables: factor index, p-component vector whose components are 1 or 0 if the corresponding variable
                                                            //appear negated or unnegated in the clause, vector of nodes attached to the factor.
                                                            //output: 1 if the operation could be done, 0 if not (the factor already exists).
    
    
    bool check(vector<bool> &);                             //this method checks that a given configuration verifies all the clauses.
                                                            //output: 1 if it does, 0 if it does not.

    void graphStructure();                                  //this method prints the structure of the graph

    
};



//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------------------------------------// class representing the BP messages
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Messages{
public:
    
    Messages(Graph &);

    int N;
    int M;
    Graph G;

    
    vector < vector <double> > xi_NodeToFac;                    //xi_NodeToFac contains the messages from node i to factor a expressing the cavity extimation
                                                                //of the probability that x_i is equal to Ja_i in the graph where a is removed.
    vector < vector <double> > Hxi_FacToNode;                   //Hxi_FacToNode contains the messages from factor b to node i expressing the cavity extimation
                                                                //of the probability that x_i is equal to Jb_i in the graph where i only appears in factor b.
    //Loosely speaking
    //xi_NodeToFac[i][a] is the message from node i to factor a
    //while
    //Hxi_FacToNode[b][i] is the message from node b to factor i
    
    vector <double> p_bias;                                     //p_bias is the prior probability that a variable is 1.

    vector < vector <double> > marginal;                        //marginal probability for each spin
  
    
    
    void initMessages();                                        //this method initializes xi messages to umbiased values. p_biases are set
                                                                //to the 0.5 and the marginals are set to the uniform distribution as well.
                                                                //it is invoked inside the constructor

    //xiUpdate and HxiUpdate are the two functions that implement a BP sweep on the whole graph:
    void HxiUpdate();
    void xiUpdate();

    bool nodeMarginals(int);                                    //this method computes (do not print) node marginals
                                                                //input variables:
                                                                //verbose: set it to 1 to print error messages (if any) when a node receives
                                                                //conflicting messages from neighbouring factors
                                                                //output: 1 if no conflicting messages are found, 0 otherwise.

    //these are the printing functions
    void etaState();                                            //this method returns the state of the eta's, i.e. messages NodeToFac
    void nuState();                                             //this method returns the state of the nu's, i.e. messages FacToNode
    void marginalState();                                       //this method returns the state of the marginals

    
    
};



//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------------------------------------// class reprensenting BP iteration
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class BP{
public:
 
    BP(Graph &, bool);
 
    int N;
    bool verbose;                                               //set this variable to 1 to have a verbose version of the methods of the class
    
    Graph G;
    Messages mess;
 
    vector <double> prev_marginal;                              //this vector is used to store the values of the k=0 component of the marginals at time t-1
                                                                //we need to keep memory of them when iterating the BP equation to check
                                                                //their convergence.
    
    bool BPsweep();                                             //this method updates all the nu's and all the eta's in the graph.
                                                                //output: 0 if at least one node receives conflicting messages, 1 otherwise
    
    void BPprint();                                             //this method prints the BP messages and the marginals

    bool BPiteration(double, int);                              //this method iterates BP equations by calling BP_sweep until convergence
                                                                //input variables:
                                                                //eps     : this value sets the convergence quality. set it to 10^-3.
                                                                //T       : maximum iteration time. set it to ~ N.
                                                                //output: 0 if at least one node receives conflicting messages, or the algorithm did not converge;
                                                                //1 otherwise
    
    void initPreviousMarginals();                               //this method initializes prev_marginals to zeros.
    
    void storePreviousMarginals();                              //this method stores the k=0 component of the marginals at time t-1 in the vector prev_marginal
    
    double compareMarginals();                                  //this method compares marginals at time t and marginal at time t-1
                                                                //the output is the maximum value of the absolute value of the differences between
                                                                //the k-components of the marginals at time t and t-1.
    
    void initDecimation(vector<int> &, vector<bool> &);         //this method set an hard bias on some specified variables
                                                                //the first vector contains the biased variables
                                                                //the second vector contais the colours toward which the biased nodes are biased.
                                                                //it also fills the vectors fixedSpins and fixedValues and notFixedSpins.
                                                                //if no specified variable are given, the last vector is made by all the variables of the system
    
    vector <int> fixedSpins;                                    //this vector contains the spin that get frozen (decimated) along the decimation process.
                                                                //it is filled by the method initDecimation and by the method fixSpins during the decimation process.
    
    vector <bool> fixedValues;                                  //this vector contains the values of the spin that get frozen (decimated) along the decimation process.
                                                                //it is filled by the method initDecimation and by the method fixSpins during the decimation process.
    
    vector <int> notFixedSpins;                                 //this vector contains the indices of the spin that are not frozen.
                                                                //at t=0, it is formed by all the spins, unless we decide to run the decimation process
                                                                //from a specific set of fixed nodes by feeding initDecimation with it.
                                                                //As t increases and the decimation process continues, its size decreases.
                                                                //it is filled by the method initDecimation and by the method fixSpins during the decimation process.

    void setHardBias(vector<int>&, vector<bool>&);              //set an hard bias on the nodes specified by the first vector according to the value specified in the second.
                                                                //default input vectors are empty vectors.
                                                                //input variables: vector of nodes to be biased, vector of colors to towards which node biases have to be biased
    
    void findMostBiased(vector<int>&, vector<bool>&);           //after having ran the BP equation till convergence, we find the most biased variables.
                                                                //by this we mean variables for which | p[0]-p[1] | > 0.999 or,if none, the variable with the largest
                                                                //absolute value of the difference p[0]-p[1].
                                                                //the search is carefully made only on the variables that have not been fixed yet.
                                                                //this method is called inside BPguidedDecimation. 

    bool warningDecimation(int &);                              //this method runs BPsweep twice in sequence to see if the spin that we decimate are able to
                                                                //fix other spins in the graph. it returns 0 if at least one node receives conflicting messages, 1 otherwise
                                                                //input variable: dummy variable set to 0 that is then used in BPguidedDecimation to count the number of time
                                                                //steps spent iterating the warnings.
                                                                //output: 1 if no contraddictions are found, 0 otherwise.
    
    void BPguidedDecimation(double, int, int);                  //this method uses the strategy of a decimation based on the hints given by BP to find solutions to a SAT instance
                                                                //at each time step (external time step) it runs BP until convergence and looks for the most biased variables.
                                                                //it fixes this(these) variable(s), runs warning propagation to see if the fixing of this variable implies the
                                                                //fixing of other variables, and then it continues.
                                                                //it stops or when contraddictions are found, or when a solution is found.
                                                                //input variables: eps, #of internal time steps (the number of time steps for which BP is iterated at most to find)
                                                                //and #of external time steps, which, in order to find complete solutions, should be O(N).
    
};




//------------------------------------------------------------------------------------------------------------------------------------------------------------------// useful functions

//this function defines allows to fill a vector in one single line
//it has been downloaded from https://gist.github.com/pablomtz/5577626
//examples:
//vector<int> v = make_vector<int>() << 0 << 0 << 1 << 1 << 1 << 0 << 0 << 1 << 1;

template <typename T>
class make_vector {
 public:
  typedef make_vector<T> my_type;
  my_type& operator<< (const T& val) {
    data_.push_back(val);
    return *this;
  }
  operator std::vector<T>() const {
    return data_;
  }
 private:
  std::vector<T> data_;
};

//this function print all the elements of a vector.
template <class T>
void vec_print(vector<T>& vec){
    for (int i=0; i<vec.size(); i++)
        cout << vec[i] << ' ';
    cout << endl;
}

