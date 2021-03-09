#define _USE_MATH_DEFINES 

#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <vector>
#include <fstream>
#include <string> 
#define M_PI       3.14159265358979323846

using std::cout;
using std::vector;
using namespace std::chrono;
using namespace std;

double A[1000000];
double B[1000000];
typedef double (*function)(double);

double task1(function func, double x1, double x2, vector<int> _split) {
	vector<int> vec = _split;
	cout.precision(20);

	double f1 = func(x1);
	double f2 = func(x2);

	//метод трапеций вариант c включённым векторизатором №17
	cout << "with vectorizator\n";
	double sum = 0;

	for (int count_iter : vec) {
		const double dx = (x2 - x1) / count_iter;
		high_resolution_clock::time_point t1 = high_resolution_clock::now();

		for (int i = 1; i < count_iter; ++i)
			B[i] = x1 + dx * i;
		for (int i = 1; i < count_iter; ++i) {
			A[i] = double(6) / sqrt(B[i] * (2 - B[i]));
			//A[i] = func(B[i]);
		}
		for (int i = 1; i < count_iter; ++i) {
			sum += A[i];
		}
		sum = (sum + (f1 + f2) / 2) * dx;
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> duration = (t2 - t1);
		cout << "iterations :" << count_iter << "\nsum: " << sum << "\nprecision: " << abs(M_PI - sum) << "\ntime: " << duration.count() << "\n\n";
	}
	return sum;
}

double task2(function func, double x1, double x2, vector<int> _split) {
	vector<int> vec = _split;
	cout.precision(20);

	double f1 = func(x1);
	double f2 = func(x2);

	//метод трапеций вариант с отключением векторизации вариант №17
	cout << "with vectorizator\n";
	double sum = 0;

	for (int count_iter : vec) {
		const double dx = (x2 - x1) / count_iter;
		high_resolution_clock::time_point t1 = high_resolution_clock::now();

#pragma loop(no_vector)
		for (int i = 1; i < count_iter; ++i)
			B[i] = x1 + dx * i;
#pragma loop(no_vector)
		for (int i = 1; i < count_iter; ++i) {
			A[i] = double(6) / sqrt(B[i] * (2 - B[i]));
			//A[i] = func(B[i]);
		}
		for (int i = 1; i < count_iter; ++i) {
			sum += A[i];
		}
		sum = (sum + (f1 + f2) / 2) * dx;
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> duration = (t2 - t1);
		cout << "iterations :" << count_iter << "\nsum: " << sum << "\nprecision: " << abs(M_PI - sum) << "\ntime: " << duration.count() << "\n\n";
	}
	return sum;
}

double task3(function func, double x1, double x2, vector<int> _split) {
	vector<int> vec = _split;
	cout.precision(20);

	double f1 = func(x1);
	double f2 = func(x2);

	//метод трапеций вариант с параллелизатором №17
	cout << "with parallelizator\n";
	double sum = 0;

	for (int count_iter : vec) {
		const double dx = (x2 - x1) / count_iter;
		high_resolution_clock::time_point t1 = high_resolution_clock::now();

#pragma loop(hint_parallel(4))
		for (int i = 1; i < count_iter; ++i)
			B[i] = x1 + dx * i;
#pragma loop(hint_parallel(4))
		for (int i = 1; i < count_iter; ++i) {
			A[i] = double(6) / sqrt(B[i] * (2 - B[i]));
			//A[i] = func(B[i]);
		}
		for (int i = 1; i < count_iter; ++i) {
			sum += A[i];
		}
		sum = (sum + (f1 + f2) / 2) * dx;
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> duration = (t2 - t1);
		cout << "iterations :" << count_iter << "\nsum: " << sum << "\nprecision: " << abs(M_PI - sum) << "\ntime: " << duration.count() << "\n\n";
	}
	return sum;
}

double _SUM;
const int COUNT_THREADS = 4;

int run_integral(function func, double x1, double x2, int count_iter, int start_inx) {
	std::cout.precision(20);

	double dx = (x2 - x1) / count_iter;
	double sum = 0;
	if (start_inx == 1) {
		sum += (func(x1) + func(x2)) / 2;
	}
#pragma loop(no_vector)
	for (double value = x1 + start_inx * dx; value < x2; value += COUNT_THREADS * dx) {
		sum += double(6) / sqrt(value * (2 - value));
		//sum += func(value);
	}
	sum *= dx;

	_SUM += sum;
	return 0;
}

double task4(function func, double x1, double x2, vector<int> _split) {
	vector<int> splitting = _split;
	const int num_iterations = 5;
	duration<double> duration[num_iterations];

	vector<double> sum_split;
	for (int i = 0; ; ++i) {
		if (i == num_iterations) {
			sum_split.push_back(_SUM);
			_SUM = 0;
			break;
		}
		else if (i != 0) {
			//_SUM += func(x1) + func(x2);
			sum_split.push_back(_SUM);
			_SUM = 0;
		}
		high_resolution_clock::time_point startpotok = high_resolution_clock::now();
		thread th1(run_integral, func, x1, x2, splitting[i], 1);
		thread th2(run_integral, func, x1, x2, splitting[i], 2);
		thread th3(run_integral, func, x1, x2, splitting[i], 3);
		thread th4(run_integral, func, x1, x2, splitting[i], 4);
		th1.join();
		th2.join();
		th3.join();
		th4.join();
		high_resolution_clock::time_point endpotok = high_resolution_clock::now();
		duration[i] = endpotok - startpotok;
	}
	cout << endl;
	cout << scientific;
	cout << setprecision(20);
	for (int i = 0; i < num_iterations; ++i) {
		std::cout << "iterations :" << splitting[i] << "\nsum: " << sum_split[i] << "\nprecision: " << abs(M_PI - sum_split[i]) << "\ntime: " << duration[i].count() << "\n\n";
	}
	return 0;
}

double doptask(function func, double x1, double x2, vector<int> _split) {
	vector<int> vec = _split;
	cout.precision(20);

	double f1 = func(x1);
	double f2 = func(x2);

	//метод трапеций вариант без векторизатора №17
	cout << "with vectorizator\n";
	double sum = 0;

	for (int count_iter : vec) {
		double dx = (x2 - x1) / count_iter;
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		//#pragma loop(no_vector)
#pragma loop(hint_parallel(4))
		for (double value = x1 + dx; value < x2; value += dx) {
			//sum += func(value);
			sum += double(6) / sqrt(value * (2 - value));
		}
		sum = (sum + (f1 + f2) / 2) * dx;
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> duration = (t2 - t1);
		cout << "iterations :" << count_iter << "\nsum: " << sum << "\nprecision: " << abs(M_PI - sum) << "\ntime: " << duration.count() << "\n\n";
	}
	return sum;
}

int main()
{
	auto func = [](double x) {return double(6) / sqrt(x * (2 - x)); };
	double x1 = 0.5;
	double x2 = 1;
	vector<int> vec = { 100, 1000, 10000, 100000, 1000000 };
	cout << "TASK1 RESULTS:\n\n";
	task1(func, x1, x2, vec);
	cout << "TASK2 RESULTS:\n\n";
	task2(func, x1, x2, vec);
	cout << "TASK3 RESULTS:\n\n";
	task3(func, x1, x2, vec);
	cout << "TASK4 RESULTS:\n\n";
	task4(func, x1, x2, vec);
	cout << "DOP_TASK RESULTS:\n\n";
	doptask(func, x1, x2, vec);
}
