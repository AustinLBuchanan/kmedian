#include "models.h"

void solve_LP_relaxation(const vector<vector<int> > &c, int k)
{
	solve_Gurobi_model(c, k, false);
}

void solve_kMedian_with_IP(const vector<vector<int> > &c, int k)
{
	solve_Gurobi_model(c, k, true);
}

void solve_Gurobi_model(const vector<vector<int> > &c, int k, bool solveIP)
{
	int n = int(c[0].size());

	try
	{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		model.set(GRB_DoubleParam_TimeLimit, 3600.); // 1 hour time limit
		model.set(GRB_DoubleParam_NodefileStart, 4); // use at most 4 GB RAM
		model.set(GRB_IntParam_Method, 3);			 // use concurrent method (3) to solve root LP
		model.set(GRB_DoubleParam_MIPGap, 0);		 // force Gurobi to prove optimality

		// create n^2 variables, each var x[i][j] has objective coefficient c[i][j]
		GRBVar** x = new GRBVar*[n];
		for (int i = 0; i < n; ++i)
		{
			x[i] = new GRBVar[n];
			for (int j = 0; j < n; ++j)
				x[i][j] = model.addVar(0., 1., c[i][j], GRB_CONTINUOUS); // Gurobi assumes minimization objective by default.
		}
		model.update();

		if (solveIP)
			for (int j = 0; j < n; ++j)
				x[j][j].set(GRB_CharAttr_VType, GRB_BINARY);

		// pick k different medians
		GRBLinExpr expr = 0;
		for (int j = 0; j < n; ++j)
			expr += x[j][j];
		model.addConstr(expr == k);

		// each i is assigned to exactly one j
		for (int i = 0; i < n; ++i)
		{
			GRBLinExpr constr = 0;
			for (int j = 0; j < n; ++j)
				constr += x[i][j];

			model.addConstr(constr == 1);
		}

		// coupling constraints. if j is closed, don't assign i to it
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				if (i != j)
					model.addConstr(x[i][j] <= x[j][j]);

		model.optimize();
		double LB = model.get(GRB_DoubleAttr_ObjBound);
		cerr << LB << " ";

		if (solveIP)
		{
			cout << "Gurobi says to pick these vertices as medians: ";
			for (int i = 0; i < n; ++i)
				if (x[i][i].get(GRB_DoubleAttr_X) > 0.5)
					cout << i << " ";
			cout << endl;
		}

	}
	catch (GRBException e) {
		printf("Error code = %d\n", e.getErrorCode());
		printf("%s\n", e.getMessage().c_str());
	}
	catch (const char* msg) {
		printf("Exception with message : %s\n", msg);
	}
	catch (...) {
		printf("Exception during optimization\n");
	}
}


