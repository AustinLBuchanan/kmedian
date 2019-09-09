#include "gurobi_c++.h"
#include "models.h"
#include <chrono>

int main(int argc, char *argv[])
{
	if (argc < 3) 
	{
		cout << "\t Codes for solving random k-median instance on n points." << endl;
    cout << "\t Solves LP relaxation and then finds IP optimal solution (with Gurobi)." << endl;
    cout << "\t Solves Lagrangian bound (with Shor's r-algorithm) and performs Lagrangian-based B&B." << endl;
		cout << "\t Usage: n k " << endl;
		return 0;
	}
	
	int n = atoi(argv[1]);
	int k = atoi(argv[2]);
	cout << "Read n = " << n << ", k = " << k << endl;
	cerr << n << " " << k << " ";
	if ( !(1<=k && k<=n) ) 
	{
		cerr << "Error: must have 1 <= k <= n." << endl;
		return 0;
	}

	//	generate random k-median instance 
	//	for n-by-n grid with rectilinear distances
	vector< vector<int> > c = generateRandomKMedianInstance(n);
	
	// solve LP relaxation using Gurobi
	auto lp_start = chrono::steady_clock::now();
	solve_LP_relaxation(c, k);
	chrono::duration<double> lp_duration = chrono::steady_clock::now() - lp_start;
	cerr << lp_duration.count() << " ";

	// solve k-median using IP/Gurobi
	auto ip_start = chrono::steady_clock::now();
	solve_kMedian_with_IP(c, k);
	chrono::duration<double> ip_duration = chrono::steady_clock::now() - ip_start;
	cerr << ip_duration.count() << " ";

	// solve Lagrangian relaxation model using Lykhovyd's implementation of Shor's r-algorithm
	auto lagrangian_start = chrono::steady_clock::now();
	solve_Lagrangian_relaxation(c, k);
	chrono::duration<double> lagrangian_duration = chrono::steady_clock::now() - lagrangian_start;
	cerr << lagrangian_duration.count() << endl;

	// solve k-median using Lagrangian-based B&B
	auto bb_start = chrono::steady_clock::now();
	solve_kMedian_with_Lagrangian(c, k);
	chrono::duration<double> bb_duration = chrono::steady_clock::now() - bb_start;
	cerr << bb_duration.count() << endl;

	return 0; 
}
