#include "models.h"
#include "ralg.h"
#include <algorithm>  // used for sorting
#include <queue>	// priority queue used for B&B 

double solve_Lagrangian_relaxation(const vector<vector<int> > &c, int k)
{
	int n = int(c[0].size());
	bbnode b;
	b.alpha = new double[n];
	for (int i = 0; i < n; ++i)	b.alpha[i] = 1.; // whatever
	b.F = vector<int>(n, -1);
	vector<bool> bestS;
	int branchvar;

	return solve_Lagrangian_relaxation(c, k, b, bestS, branchvar);
}

double solve_Lagrangian_relaxation(const vector<vector<int> > &c, int k, bbnode &b, vector<bool> &bestS, int &branchvar)
{
	double LB = -INFINITY;
	int n = int(c[0].size());

	vector<double> C(n, 0);
	vector<vector<double>> c_bar(n, vector<double>(n));
	vector<bool> S(n); // boolean representation of the current solution (just the medians)
	double * bestAlpha = new double[n];

	auto cb_grad_func = [&c, k, &C, &c_bar, &S, &LB, n, &bestS, &b, &branchvar](const double* alpha, double& f_val, double* grad)
	{
		// compute L(F_0,F_1,alpha) for fixed alpha
		solveInnerProblem(alpha, k, c, c_bar, C, grad, f_val, S, b.F, branchvar);
		if (f_val > LB)
		{
			LB = f_val;
			bestS = S;
		}
		return true;
	};

	ralg_options opt = defaultOptions; opt.is_monotone = false;
	if (count(b.F.begin(), b.F.end(), -1) == n) // root node
	{
		opt.itermax = 1000;
		opt.output_iter = 100;
	}
	else // not the root node
	{
		opt.itermax = 200;
		opt.output = false;
	}

	// lower bound from lagrangian
	LB = ralg(&opt, cb_grad_func, n, b.alpha, bestAlpha, RALG_MAX);

	delete[] bestAlpha;
	return LB;
}

void solveInnerProblem(const double* alpha, int k, const vector<vector<int>>& c, vector<vector<double>>& c_bar,
	vector<double>& C, double* grad, double& f_val, vector<bool>& S, const vector<int> &F, int &branchvar)
{
	int n = int(c[0].size());

	// reset medians
	for (int i = 0; i < n; ++i)
		S[i] = false;

	// update c_bar
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			c_bar[i][j] = c[i][j] - alpha[i];

	// recompute C_j, the contribution to objective of opening facility j
	for (int j = 0; j < n; ++j)
	{
		C[j] = c_bar[j][j];
		for (int i = 0; i < n; ++i)
			if (i != j && c_bar[i][j] < 0)
				C[j] += c_bar[i][j];
	}

	// select k smallest
	vector<int> C_indices(C.size());
	for (size_t i = 0; i < C.size(); ++i)
		C_indices[i] = i;
	std::sort(C_indices.begin(), C_indices.end(), [&C, &F](int i1, int i2) {
		if (F[i1] == F[i2])
			return C[i1] < C[i2];
		if (F[i1] == 1) return true;
		if (F[i2] == 0) return true;
		return false;
	});

	// compute inner problem function value (f_val)
	f_val = 0.;
	for (int i = 0; i < n; ++i)
		f_val += alpha[i];

	for (int j = 0; j < k; ++j)
	{
		int v = C_indices[j];
		f_val += C[v];
		S[v] = true;
	}

	// find variable to branch on next
	branchvar = -1;
	for (int j = k - 1; j >= 0; --j)
	{
		int v = C_indices[j];
		if (F[v] != 1)
			branchvar = v;
	}

	// compute subgradient
	for (int i = 0; i < n; ++i)
	{
		grad[i] = 1.;
		for (int j = 0; j < k; ++j)
			if (i == C_indices[j] || c_bar[i][C_indices[j]] < 0)
				grad[i] -= 1;
	}
}

double obj(const vector<vector<int>> &c, const vector<bool> &S, int k)
{
	if (count(S.begin(), S.end(), true) != k) return DBL_MAX;

	double objval = 0;
	int n = int(c[0].size());
	vector<int> intS;

	// cost to assign S -> S
	for (int i = 0; i < n; ++i) 
	{
		if (S[i])
		{
			intS.push_back(i);
			objval += c[i][i];
		}
	}

	// cost to assign V\S -> S
	for (int i = 0; i < n; ++i) 
	{
		if (S[i]) continue;

		double minval = DBL_MAX;
		for (int p = 0; p < k; ++p)
		{
			int j = intS[p];
			if (c[i][j] < minval)
				minval = c[i][j];
		}
		objval += minval;
	}

	return objval;
}

void solve_kMedian_with_Lagrangian(const vector<vector<int> > &c, int k)
{
	int n = int(c[0].size());

	auto cmp = [](pair<double, bbnode> left, pair<double, bbnode> right) { return left.first < right.first; };
	priority_queue< pair<double, bbnode>, vector <pair<double, bbnode>>, decltype(cmp) > B(cmp);

	// construct root node and add to BB tree
	bbnode rootnode;
	rootnode.sizeF0 = 0;
	rootnode.sizeF1 = 0;
	rootnode.F.resize(n, -1);
	rootnode.alpha = new double[n];
	for (int i = 0; i < n; ++i)
		rootnode.alpha[i] = 1.; // whatever

	B.push(make_pair(0., rootnode)); // 0 is just a placeholder. doesn't matter for root node.

	// perform BB
	vector<bool> S;			// feasible solution obtained from current BB node
	vector<bool> bestS;		// best feasible solution found so far (just the medians)
	double UB = DBL_MAX;	// objective value of bestS
	int nodesExplored = 0;
	int nodesTotal = 1;
	int nodesSaved = 0;

	while (!B.empty())
	{
		bbnode b = B.top().second;
		B.pop();

		if (b.parentLB > UB) {
			nodesSaved++;
			continue; // this can be pruned without computing new lagrangian bound
		}
		nodesExplored++;
		
		// compute a lower bound, from Lagrangian
		double LB = solve_Lagrangian_relaxation(c, k, b, S, b.branchvar);

		double tempUB = obj(c, S, k);
		if (tempUB < UB)
		{
			bestS = S;
			UB = tempUB;
		}
		if (nodesExplored == 1)
		{
			cout << "\t" << "  #BB nodes" << "\t" << endl;
			cout << "\t" << "expl" << "\t" << "total" << "\t" << "L(parent)" << "\t" << "L(F0,F1)" << "\t" << "UB" << "\t" << "|F0|" << "\t" << "|F1|" << "\t" << "obj(S)" << endl;
			cout << "\t" << nodesExplored << "\t" << nodesTotal << "\t" << "NA" << "\t\t" << LB << "\t\t" << UB << "\t" << b.sizeF0 << "\t" << b.sizeF1 << "\t" << tempUB << endl;
		}
		else
			cout << "\t" << nodesExplored << "\t" << nodesTotal << "\t" << b.parentLB << "\t\t" << LB << "\t\t" << UB << "\t" << b.sizeF0 << "\t" << b.sizeF1 << "\t" << tempUB << endl;


		if (LB < UB && b.branchvar>=0)
		{
			int v = b.branchvar;
			// do variable fixing? 
			// TODO

			// create left branch
			bbnode leftbranch;
			leftbranch.parentLB = LB;
			leftbranch.sizeF0 = b.sizeF0 + 1;
			leftbranch.sizeF1 = b.sizeF1;
			leftbranch.F = b.F;
			leftbranch.F[v] = 0;	// fix x_vv = 0
			leftbranch.alpha = new double[n];
			for (int i = 0; i < n; ++i)
				leftbranch.alpha[i] = b.alpha[i];

			// create right branch
			bbnode rightbranch;
			rightbranch.parentLB = LB;
			rightbranch.sizeF0 = b.sizeF0;
			rightbranch.sizeF1 = b.sizeF1 + 1;
			rightbranch.F = b.F;
			rightbranch.F[v] = 1;	// fix x_vv = 1
			rightbranch.alpha = new double[n];
			for (int i = 0; i < n; ++i)
				rightbranch.alpha[i] = b.alpha[i];

			// add branches to B
			B.push(make_pair(-LB, leftbranch));	// best bound (estimated, based on parent)
			B.push(make_pair(-LB, rightbranch));

			//int depth = b.sizeF0 + b.sizeF1 + 1;
			//B.push(make_pair(-depth, leftbranch));    // breadth first search
			//B.push(make_pair(-depth, rightbranch));

			//B.push(make_pair(depth, leftbranch));	// depth first search
			//B.push(make_pair(depth, rightbranch));

			//B.push(make_pair(LB, leftbranch));	// worst bound (estimated, based on parent)
			//B.push(make_pair(LB, rightbranch));

			//B.push(make_pair(rand(), leftbranch));	// random search
			//B.push(make_pair(rand(), rightbranch));
			nodesTotal += 2;
		}
	}
	cout << "Final BB solution : ";
	for (int i = 0; i < n; ++i)
		if (bestS[i])
			cout << i << " ";
	cout << " has objective value " << UB << " ?= " << obj(c, bestS, k) << endl;
	cout << "Visited this many BB nodes = " << nodesExplored << endl;
	cout << "Created this many BB nodes = " << nodesTotal << endl;
	cout << "Avoided Lagrangian at this many BB nodes = " << nodesSaved << endl;
	return;
}
