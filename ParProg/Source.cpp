#include <iostream> 
#include <ctime>  
#include <mpi.h>

using namespace::std;

int  curr_rank_proc; 
int  num_of_procs;
double  time_seq_work_alg = 0;
double  time_pp_work_alg = 0;



double** Create_and_init_matr(int size_matr_row, int size_matr_column)
{
	if (size_matr_row < 1 || size_matr_column < 1)
		return NULL;

	double** matr;



	matr = new double* [size_matr_row]; 
	for (int i = 0; i < size_matr_row; i++)
		matr[i] = new double[size_matr_column];


	srand((unsigned)time(NULL));
	if (size_matr_row > 10 || size_matr_column > 10)
	{
		for (int i = 0; i < size_matr_row; i++)
			for (int j = 0; j < size_matr_column; j++)
				matr[i][j] = (double)(rand() % 100) - 50.0; 
	}
	else
	{
		for (int i = 0; i < size_matr_row; i++)
		{
			cout << "Input double elements of " << i << " string " << endl;
			for (int j = 0; j < size_matr_column; j++)
				cin >> matr[i][j];
		}
	}
	return matr;
}


double* Convert_matrix_to_vector(double** matr, int size_row, int size_column)
{
	if (matr != NULL)
	{
		double* vec = new double[size_row * size_column];
		int index_vec = 0;
		for (int i = 0; i < size_row; i++)
			for (int j = 0; j < size_column; j++)
			{
				vec[index_vec] = matr[i][j];
				index_vec++;
			}
		return vec;
	}
	return NULL;
}

void Show_matr(double** matr, int size_matr_row, int size_matr_column)
{
	if (matr != NULL)
		for (int i = 0; i < size_matr_row; i++)
		{
			for (int j = 0; j < size_matr_column; j++)
				cout << matr[i][j] << " ";
			cout << endl;
		}
}
int main(int argc, char* argv[])
{
	int size_row = 5, size_column = 5;
	double** test_matr = NULL; 
	double* matrix_as_vector = NULL;

	double sum_el_seq = 0;
	double sum_el_pp = 0;

	double end_time_of_seq_alg = 0;
	double start_time_of_seq_alg = 0;
	double end_time_of_pp_alg = 0;
	double start_time_of_pp_alg = 0;


	double partial_summ = 0;
	int size_el = 1;
	int size_work_of_proc = 0; 

	MPI_Status stat;


	/* Начало MPI кода */

	MPI_Init(&argc, &argv); 
							

	MPI_Comm_size(MPI_COMM_WORLD, &num_of_procs); 
												  
	MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank_proc);  

	if (num_of_procs < 1)
	{
		cout << "Incorrect number of processes (at least 1 must be)" << endl;
		return 0;
	}

	if (curr_rank_proc == 0)
	{
		cout << "Input size of row: " << endl;
		cin >> size_row;
		cout << "Input size of column: " << endl;
		cin >> size_column;

		size_el = size_row * size_column;
		test_matr = Create_and_init_matr(size_row, size_column);
		matrix_as_vector = Convert_matrix_to_vector(test_matr, size_row, size_column); 

		if (test_matr == NULL || matrix_as_vector == NULL)
		{
			cout << "Incorrect input data, try again" << endl;
			return 0;
		}

		if (size_row * size_column < 1000)
		{
			cout << "Current matrix:" << endl;
			Show_matr(test_matr, size_row, size_column);
		}
		cout << endl;

		/* Подсчет суммы всех элементов матрицы (последовательная версия алгоритма): */

		start_time_of_seq_alg = clock();

		for (int i = 0; i < size_row; i++)
			for (int j = 0; j < size_column; j++)
				sum_el_seq += test_matr[i][j];

		end_time_of_seq_alg = clock();
		time_seq_work_alg = end_time_of_seq_alg - start_time_of_seq_alg;

		cout << "Sum of all elements in matrix is  " << sum_el_seq << " " << endl;
		cout << "Spent time on the implementation of this algorithm (Sequence version)" << time_seq_work_alg << " ms " << endl;
		cout << endl;
		cout << "Num of procs: " << num_of_procs << endl;



		start_time_of_pp_alg = MPI_Wtime();

		size_work_of_proc = size_el / num_of_procs; 

		for (int i = 1; i < num_of_procs; i++)
			MPI_Send(&size_el, 1, MPI_INT, i, 1, MPI_COMM_WORLD);

		for (int i = 1; i < num_of_procs; i++)
			MPI_Send(matrix_as_vector + size_work_of_proc * (i - 1), size_work_of_proc, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);

		cout << "Process with rang " << curr_rank_proc << " start own job" << endl;
		
		for (int i = size_work_of_proc * (num_of_procs - 1); i < size_el; i++)
			partial_summ += matrix_as_vector[i];

		sum_el_pp = partial_summ;

		for (int i = 1; i < num_of_procs; i++)
		{
			MPI_Recv(&partial_summ, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &stat);
			sum_el_pp += partial_summ;
		}

		end_time_of_pp_alg = MPI_Wtime();

		time_pp_work_alg = (end_time_of_pp_alg - start_time_of_pp_alg) * 1000;
		cout << "Sum of all elements in matrix is (Parallel version): " << sum_el_pp << endl;
		cout << "Spent time on the implementation of this algorithm (Parallel version)" << time_pp_work_alg << " ms" << endl;

		if (sum_el_pp == sum_el_seq)
			cout << "Results of parallel and sequence versions are identical! " << endl;
		else
			cout << "Results of parallel and sequence versions are not identical! " << endl;

		if (time_pp_work_alg <= time_seq_work_alg)
			cout << "Parallel version faster, then sequence" << endl;
		else
			cout << "Sequence version faster, then parallel" << endl;

		delete matrix_as_vector;
		delete[] test_matr;
	}
	else
	{
		double* recv_matrix_as_vector;

		MPI_Recv(&size_el, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);

		size_work_of_proc = size_el / num_of_procs;
		recv_matrix_as_vector = new double[size_work_of_proc];
		MPI_Recv(recv_matrix_as_vector, size_work_of_proc, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &stat);

		cout << "Process with rang " << curr_rank_proc << " start own job" << endl;

		for (int i = 0; i < size_work_of_proc; i++)
			partial_summ += recv_matrix_as_vector[i];

		MPI_Send(&partial_summ, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		delete[] recv_matrix_as_vector;
	}

	MPI_Finalize();

	/* Конец MPI кода */
	return 0;
}