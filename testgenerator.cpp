#include "models.h"

vector< vector<int> > generateRandomKMedianInstance(int n)
{
	vector<int> x_coordinate(n);
	vector<int> y_coordinate(n);
	for (int i = 0; i < n; ++i)
	{
		x_coordinate[i] = rand() % n;
		y_coordinate[i] = rand() % n;
	}

	//draw_instance(x_coordinate, y_coordinate);

	vector< vector<int> > c(n, vector<int>(n));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			// compute c[i][j] as rectlinear distance
			if (x_coordinate[i] > x_coordinate[j])
				c[i][j] = x_coordinate[i] - x_coordinate[j];
			else
				c[i][j] = x_coordinate[j] - x_coordinate[i];

			if (y_coordinate[i] > y_coordinate[j])
				c[i][j] += y_coordinate[i] - y_coordinate[j];
			else
				c[i][j] += y_coordinate[j] - y_coordinate[i];

		}
	}

	return c;
}

void draw_instance(const vector<int> &x, const vector<int> &y)
{
	int n = x.size();
	cout << endl << "~~ drawing points on n-by-n grid ~~ \n";
	
	vector< vector<int> > numGridPoints(n, vector<int>(n, 0));
	for (int i = 0; i < n; ++i)
		numGridPoints[x[i]][y[i]]++;

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
			cout << numGridPoints[i][j] << " ";
		cout << endl;
	}
}

