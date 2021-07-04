#include<iostream>

using namespace std;

void GaussEli(double* A, double* b, int n)
{
	for (int k = 1; k <= n - 1; k++)
	{
		double a = A[(k - 1) * n + (k - 1)];
		for (int i = k + 1; i <= n; i++)
		{
			double aa = A[(i - 1) * n + (k - 1)] / a;
			for (int j = k; j <= n; j++)
			{
				A[(i - 1) * n + (j - 1)] -= A[(k - 1) * n + (j - 1)] * aa;
			}
			b[i - 1] -= b[k - 1] * aa;
		}
	}
}

void GaussEliPivot(double* A, double* b, int* piv, int n)
{
	for (int p = 1; p <= n; p++)
		piv[p - 1] = p - 1;

	for (int k = 1; k <= n - 1; k++)
	{

		int MaxID = k - 1;
		for (int p = k + 1; p <= n; p++)
		{
			if (abs(A[piv[MaxID] * n + (k - 1)]) < abs(A[piv[p - 1] * n + (k - 1)]))
				MaxID = p - 1;
		}
		if (MaxID != (k - 1))
		{
			int temp;
			temp = piv[k - 1];
			piv[k - 1] = piv[MaxID];
			piv[MaxID] = temp;
		}
		double a = A[piv[k - 1] * n + (k - 1)];
		for (int i = k + 1; i <= n; i++)
		{
			double aa = A[piv[i - 1] * n + (k - 1)] / a;
			for (int j = k; j <= n; j++)
			{
				A[piv[i - 1] * n + (j - 1)] -= A[piv[k - 1] * n + (j - 1)] * aa;
			}
			b[piv[i - 1]] -= b[piv[k - 1]] * aa;
		}
	}
}

void SolveUy(double* U, double* y, double* x, int n)
{
	for (int i = n; i >= 1; i--)
	{
		x[i - 1] = y[i - 1] / U[(i - 1) * n + (i - 1)];
		for (int j = 1; j <= i - 1; j++)
		{
			y[j - 1] -= U[(j - 1) * n + (i - 1)] * x[i - 1];
		}
	}
}

void SolveUyPivot(double* U, double* y, double* x, int* piv, int n)
{
	for (int i = n; i >= 1; i--)
	{
		x[piv[i - 1]] = y[piv[i - 1]] / U[piv[i - 1] * n + (i - 1)];
		for (int j = 1; j <= i - 1; j++)
		{
			y[piv[j - 1]] -= U[piv[j - 1] * n + (i - 1)] * x[piv[i - 1]];
		}
	}
}

int main()
{
	int n;
	cout << "Input n=";
	cin >> n;
	double* A = new double[n * n];
	double* b = new double[n];
	double* x = new double[n];
	double* xx = new double[n];
	int* piv = new int[n];
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
		{
			cout << "Input A[" << i << "][" << j << "]=" << endl;
			cin >> A[(i - 1) * n + (j - 1)];
		}

	for (int i = 1; i <= n; i++)
	{
		cout << "Input b[" << i << "]=" << endl;
		cin >> b[i - 1];
	}

	GaussEliPivot(A, b, piv, n);
	cout << "U=" << endl;
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
			cout << A[piv[i - 1] * n + (j - 1)] << "  ";
		cout << endl;
	}
	cout << "L^(-1)b=" << endl;
	for (int i = 1; i <= n; i++)
		cout << b[piv[i - 1]] << endl;

	SolveUyPivot(A, b, x, piv, n);
	cout << "x=" << endl;
	for (int i = 1; i <= n; i++)
		cout << x[piv[i - 1]] << endl;


	/*GaussEli(A, b, n);
	cout << "U=" << endl;
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
			cout << A[(i - 1) * n + (j - 1)] << "  ";
		cout << endl;
	}
	cout << "L^(-1)b=" << endl;
	for (int i = 1; i <= n; i++)
		cout << b[i-1]<< endl;

	SolveUy(A, b, x, n);
	cout << "x=" << endl;
		for (int i = 1; i <= n; i++)
			cout << x[i - 1] << endl;

	*/
	return 0;
}