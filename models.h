#include<vector>
#include "gurobi_c++.h"
using namespace std;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//		TEST INSTANCE GENERATOR
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* Generates a random k-median instance on n points. 
	Put points on n-by-n grid and use rectilinear distances */
vector< vector<int> > generateRandomKMedianInstance(int n);

/* Draws points given x- and y-coordinates */
void draw_instance(const vector<int> &x, const vector<int> &y);



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//		SOLVE k-MEDIAN IP AND LP RELAXATION
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void solve_LP_relaxation(const vector<vector<int> > &c, int k);
void solve_kMedian_with_IP(const vector<vector<int> > &c, int k);
void solve_Gurobi_model(const vector<vector<int> > &c, int k, bool solveIP);



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//		LAGRANGIAN FUNCTIONS
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/* struct for (Lagrangian-based) branch-and-bound nodes */
struct bbnode
{
	vector<int> F;	// F[j]=0 means x_jj=0 fixed; F[j]=1 means x_jj=1 fixed; F[j]=-1 means x_jj is free/unfixed.
	double *alpha;	// lagrangian multipliers, warm started from multipliers of parent node
	int sizeF0;		// |F_0| = number of vars fixed to zero
	int sizeF1;		// |F_1| = number of vars fixed to one
	int branchvar;	// variable to branch on next
	double parentLB;// L(parent) = Lagrangian bound coming from the parent node
};

/* Compute the objective function value when S is the set of medians (S in boolean representation) */
double obj(const vector<vector<int>> &c, const vector<bool> &S, int k);

/* Compute the Lagrangian relaxation bound L(F_0,F_1) at a branch-and-bound node b=(F_0,F_1) */
double solve_Lagrangian_relaxation(const vector<vector<int> > &c, int k, bbnode &b, vector<bool> &bestS, int &branchvar);
double solve_Lagrangian_relaxation(const vector<vector<int> > &c, int k); // for root relaxation: b = (empty,empty)

/* Compute the Lagrangian relaxation bound L(F_0,F_1,alpha) when multipliers alpha are fixed */
void solveInnerProblem(const double* alpha, int k, const vector<vector<int>>& c, vector<vector<double>>& c_bar, vector<double>& C, 
	double* grad, double& f_val, vector<bool>& S, const vector<int> &F, int &branchvar);

/* Overall function for solving k-median problem with lagrangian techniques */
void solve_kMedian_with_Lagrangian(const vector<vector<int> > &c, int k);



