#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <cassert>

using namespace std;



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------// class representing a node
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Node{
public:

    Node(int p_n) : n ( p_n ) {};
  
    int n;                                                  //n is the index of the node
    vector <int> v_fac;                                     //v_fac contains the indices of factors attached to n
    
    int d;                                                  //this value is determined by the Leaf Removal Algorithm.
                                                            //it is 1 if the variable belongs to the 2-core structure, 0 otherwise
    
    int numberOfFactors();                                  //this method returns the number of factors to which the node n is attached
    
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------------// class representing a factor
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Factor{
public:
    
    Factor(int p_f, vector<bool> p_v_J, int p_p);
    
    int f;                                                  //f is the index of the factor
    int p;                                                  //p is the number of variables entering in a factor. p=3 in the 3-XORSAT problem

    vector <int> v_node;                                    //v_node contains the indices of nodes attached to f
    
    int numberOfNodes();                                    //this method returns the numer of variables to which a factor is attached. this number is p by definition

    vector <bool> v_J;                                      //v_J is a set of p binary variables specifying the SAT clause.
                                                            //consider the clause x1 | !x2 | !x3 where ! indicates the negation.
                                                            //be Jk the k-component of the vector v_J.
                                                            //in this case we set J1=0, J2=1, J3=1 and we write the clause as
                                                            //(2 J1 - 1) [ J1 - x1 ] | (2 J2 - 1) [ J2 - x2 ] | (2 J3 - 1) [ J3 - x3 ]
                                                            //in other words, we just send xk to wk = (2 Jk - 1) [ Jk - xk ]
                                                            //and we see that if Jk is 0, then wk = xk while if Jk=1, then wk = 1 - xk = ! xk

    double clause(int, int, int);                           //this method implements the SAT clause represented by a factor. input: the p=3 variables involved in each clause
    
    void plantedClause(int, int, int);                      //this method looks at the input variables and sets the vector v_J accordingly.
                                                            //input: the p=3 variables involved in each clause
    
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------// class representing a factor graph
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Graph{
public:
    
    Graph(int, int);
    
    int N;
    int p;
    int M;
    
    vector <Node> v;                                        //v contains all the nodes of the graph
                                                            //this vector is filled by the constructor.
    
    vector <Factor*> F;                                     //F contains pointers to all the factors of the graph.
                                                            //this vector is filled in the classes derived from Graph, that specify the model under consideration.
    
    
    int numberOfTotalFactors();                             //this method returns the total number of factors in the graph
    int numberOfTotalNodes();                               //this method returns the total number or nodes in the graph

    void factorsOfNode(int);                                //this method returns the factors attached to the input node
    void nodesOfFactor(int);                                //this method returns the nodes attached to the input factor
    
    void ErdosRenyi(int);                                   //this method build an ER graph
                                                            //input variables:
                                                            //M : # of factors
    
    void plantedErdosRenyi(int, vector<int> &);             //this method build an ER graph for which the second argument is a solution
                                                            //inputs:
                                                            //M  : # of factors
                                                            //ps : planted solution

    void graphStructure();                                  //this method prints the structure of the graph
    
    int addFactor(int, vector<bool>, vector<int>);          //the addFactor method implement the operations needed when adding a factor in the graph
                                                            //namely one needs to create a FactorSat object a, store its neighour variables and for each of them
                                                            //add the index of a as a neghbouring factor.
                                                            //input variables: factor index, p=3-component vector whose components are 0 or 1 if the corresponding variable
                                                            //appear negated or unnegated in the clause, vector of nodes attached to the factor.
                                                            //output: is 1 if the operation could be done, 0 if not (the factor already exists).
    
    
    bool check(vector<int> &);                              //this method checks that the planted solution verifies all the clauses.
                                                            //output: 1 if it does, 0 if it does not.

    
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------------------------------------------------------// class reprensenting the BP messages
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

class Messages{
public:
    
    Messages(Graph &);

    int N;
    int M;
    Graph G;

    
    vector < vector <double> > xi_NodeToFac;                    //xi_NodeToFac contains the messages from node i to factor a expressing the cavity extimation
                                                                //of the probability that x_i is equal to Ja_i in the graph where a is removed.
    vector < vector <double> > Hxi_FacToNode;                   //Hxi_FacToNode contains the messages from factors b to node i expressing the cavity extimation
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

    void nodeMarginals(int);                                    //this method computes (do not print) node marginals
                                                                //input variables:
                                                                //verbose: set it to 1 to print error messages (if any) when a node receives
                                                                //conflicting messages from neighbouring factors

    //these are the printing functions
    void etaState();                                            //this method returns the state of the eta's, ie. messages NodeToFac
    void nuState();                                             //this method returns the state of the nu's, ie. messages FacToNode
    void marginalState();                                       //this method returns the state of the marginals

    
    
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------// class reprensenting BP algorithms
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

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
    
    void BPsweep();                                             //this method updates all the nu's and all the eta's in the graph.
                                                                //it returns 0 if at least one node receives conflicting messages, 1 otherwise
                                                                //input variable:
                                                                //verbose: set it to 1 to print the messages, 0 otherwise
    
    void BPprint();                                             //this method prints the BP messages and the marginals

    bool BPiteration(double, int);                              //this method iterates BP equations by calling BP_sweep until convergence
                                                                //input variables:
                                                                //eps     : this value sets the convergence quality. set it to 10^-3.
                                                                //T       : maximum iteration time. set it to ~ N.
                                                                //output  : it returns 0 if conflicting messages are found, 1 otherwise
    
    void initPreviousMarginals();                               //this method initializes prev_marginals to zeros.
    void storePreviousMarginals();                              //this method stores the k=0 component of the marginals at time t - 1 in the vector prev_marginal
    double compareMarginals();                                  //this method compares marginals at time t and marginal at time t-1
                                                                //the output is the maximum value of the absolute value of the differences between
                                                                //the k-components of the marginals at time t and t-1.
    
    bool findFrustratedSpins();                                 //this method evaluates, for each node, all the messages from its factors (like nodeMarginals) and seeks if they contain
                                                                //contraddictory information. In this case, this spin is frustrated.
                                                                //when frustrated spins are found the method returns 0.
    
    void initDecimation(vector<int> &, vector<int> &);          //this method set an hard bias on some specified variables
                                                                //the first vector contains the biased variables
                                                                //the second vector contais the colours toward which the biased nodes are biased.
    
    vector <int> fixedSpins;                                    //this vector contains the spin that get frozen (decimated) along the decimation process.
                                                                //it is filled by the method initDecimation and by the method fixSpins during the decimation process.
    
    vector <int> fixedValues;                                   //this vector contains the values of the spin that get frozen (decimated) along the decimation process.
                                                                //it is filled by the constructor and by the method fixSpins during the decimation process.
    
    vector <int> notFixedSpins;                                 //this vector contains the indices of the spin that are not frozen.
                                                                //at t=0, it is formed by all the spins. As t increases and the decimation process continues, its size
                                                                //decreases.
                                                                //it is filled by the constructor and by the method fixSpins during the decimation process.
    
    vector <int> frustratedSpins;                               //this vector contains the spins that receive conflicting messages from the factors.
                                                                //it is filled by the method findFrustratedSpins().
    
    
    
    
    
    
    int fixSpins(int);                                          //this method fix spins values, by calling setHardBias.
                                                                //it is invoked in the method warningDecimation,
                                                                //when a node receives a clear message from at least one if its factors.
                                                                //it also fills the vector fixedSpins and, when doing this, erase nodes from the vector notFixedSpins.
                                                                //input variable:
                                                                //verbose: set it to 1 to print the nodes that get frozen
    

    
    
    
    bool warningDecimation(int);                                //this method runs BPsweep twice to see if the spin that we decimate are able to fix other spins
                                                                //in the graph.
                                                                //it returns 0 if at least one node receives conflicting messages, 1 otherwise
    
    


    void BPguidedDecimation(int, int);

    

    void findMostBiased(vector<int>&, vector<int>&);
    
    

    
    

};

 /*
 
 void setHardBias(vector<int> & , vector<int> &);            //this method defines hard biases towards color q for each node contained
 //in the first input vector and according to the color index contained in the second input vector.
 //default input vectors are empty vectors.
 //input variables: vector of nodes to be biased, vector of colors to towards which node biases have to be biased
 
 
 */

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------// useful functions

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
void vec_print(vector<int>& vec);


