
#include "mpi.h"
#include <iostream>
#include "math.h"
#include <ctime>

using namespace std;

//------ Последовательный алгоритм ------//
double* MVmul(double* pMatrix, double* pVector, const int& Size) {
	double* res = new double[Size];
	for (int i = 0; i < Size; ++i) {
		res[i] = 0;
		for (int j = 0; j < Size; j++)
			res[i] += pMatrix[Size * i + j] * pVector[j];
	}
	return res;
}
double* VVsub(double* pVectorLeft, double* pVectorRight, const int& Size) {
	double* res = new double[Size];
	for (int i = 0; i < Size; ++i) {
		res[i] = pVectorLeft[i] - pVectorRight[i];
	}
	return res;
}
double* VVsum(double* pVectorLeft, double* pVectorRight, const int& Size) {
	double* res = new double[Size];
	for (int i = 0; i < Size; ++i) {
		res[i] = pVectorLeft[i] + pVectorRight[i];
	}
	return res;
}
double ScalMul(double* pVectorLeft, double* pVectorRight, const int& Size) {
	double res = 0;
	for (int i = 0; i < Size; ++i) {
		res += pVectorLeft[i] * pVectorRight[i];
	}
	return res;
}
double* VDmul(double* pVector, double pDouble, const int& Size) {
	double* res = new double[Size];
	for (int i = 0; i < Size; ++i) {
		res[i] = pVector[i] * pDouble;
	}
	return res;
}
void ProcessInitializationGradient(double*& pVector, double*& pPrevResult, double*& pResult, double*& pPrevD, double*& pD,
	double*& pPrevG, double*& pG, int Size) {
	pPrevResult = new double[Size];
	pResult = new double[Size];
	pPrevD = new double[Size];
	pD = new double[Size];
	pPrevG = new double[Size];
	pG = new double[Size];
	for (int i = 0; i < Size; i++) {
		pPrevResult[i] = 0;
		pPrevD[i] = 0;
		pPrevG[i] = -pVector[i];
	}
}
void SoprGradSeq(double*& pMatrixGrad, double*& pVectorGrad, double*& pPrevResult, double*& pResult, double*& pPrevD, double*& pD,
	double*& pPrevG, double*& pG, double& s, double eps, int Size) {
	ProcessInitializationGradient(pVectorGrad, pPrevResult, pResult, pPrevD, pD, pPrevG, pG, Size);
	int flag = 1;
	while (flag) {
		//first
		double* pLeft = MVmul(pMatrixGrad, pPrevResult, Size);
		pG = VVsub(pLeft, pVectorGrad, Size);
		delete pLeft;
		//second
		double above = ScalMul(pG, pG, Size);
		double beyond = ScalMul(pPrevG, pPrevG, Size);
		double* pRight = VDmul(pPrevD, (double)(above / beyond), Size);
		pD = VVsub(pRight, pG, Size);
		delete pRight;
		//third
		above = ScalMul(pD, pG, Size);
		double* beyondr = MVmul(pMatrixGrad, pD, Size);
		beyond = ScalMul(pD, beyondr, Size);
		delete beyondr;
		s = -(double)(above / beyond);
		//fourth
		pRight = VDmul(pD, s, Size);
		pResult = VVsum(pPrevResult, pRight, Size);
		double sum = 0;
		for (int i = 0; i < Size; i++)
			sum += fabs(pResult[i] - pPrevResult[i]);
		delete pPrevG; delete pPrevD; delete pPrevResult;
		pPrevG = pG; pPrevD = pD; pPrevResult = pResult;
		if (sum < eps) flag = 0; 
	}
}


//------ Метод Гаусса ------//
int* pPivotPos;
int* pPivotIter;
void Show(double* pMatrix, double* pVector, int Size)
{
	cout << "\n";
	for (int i = 0; i < Size; i++)
	{
		cout << '|';
		for (int j = 0; j < Size; j++)
		{
			cout << pMatrix[Size * i + j] << " ";
		}
		cout << '|';
		cout << '=';
		cout << '|' << pVector[i] << '|' << "\n";
	}
}
int FindPivotRow(double* pMatrix, int Size, int Iter)
{
	int PivotRow = -1;
	double MaxValue = 0;
	for (int i = 0; i < Size; i++)
	{
		if ((pPivotIter[i] == -1) && (fabs(pMatrix[i * Size + Iter]) > MaxValue))
		{
			PivotRow = i;
			MaxValue = fabs(pMatrix[i * Size + Iter]);
		}
	}
	return PivotRow;
}
void SerialColumnElimination(double* pMatrix, double* pVector, int Size, int Iter, int PivotRow)
{
	double PivotValue, PivotFactor;
	PivotValue = pMatrix[PivotRow * Size + Iter];
	for (int i = 0; i < Size; i++)
	{
		if (pPivotIter[i] == -1)
		{
			PivotFactor = pMatrix[i * Size + Iter] / PivotValue;
			for (int j = Iter; j < Size; j++)
			{
				pMatrix[i * Size + j] -= PivotFactor * pMatrix[PivotRow * Size + j];
			}
			pVector[i] -= PivotFactor * pVector[PivotRow];
		}
	}

}
void SerialGaussianElimination(double* pMatrix, double* pVector, int Size)
{
	int Iter;
	int PivotRow;
	for (Iter = 0; Iter < Size; Iter++)
	{
		PivotRow = FindPivotRow(pMatrix, Size, Iter);
		pPivotPos[Iter] = PivotRow;
		pPivotIter[PivotRow] = Iter;
		SerialColumnElimination(pMatrix, pVector, Size, Iter, PivotRow);
	}
}
void SerialBackSubstitution(double* pMatrix, double* pVector, double* pResult, int Size)
{
	int RowIndex, Row;
	for (int i = Size - 1; i >= 0; --i)
	{
		RowIndex = pPivotPos[i];
		pResult[i] = pVector[RowIndex] / pMatrix[Size * RowIndex + i];
		pMatrix[RowIndex * Size + i] = 1;
		for (int j = 0; j < i; ++j)
		{
			Row = pPivotPos[j];
			pVector[Row] -= pMatrix[Row * Size + i] * pResult[i];
			pMatrix[Row * Size + i] = 0;
		}
	}
}
void SerialResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size)
{
	SerialGaussianElimination(pMatrix, pVector, Size);
	SerialBackSubstitution(pMatrix, pVector, pResult, Size);
}

//------ Параллельный алгоритм ------//
void RandomDataInitialization(double*& pMatrixGA, double*& pVectorGA, double*& pMatrixGR, double*& pVectorGR, int& Size)
{
	pMatrixGA = new double[Size * Size];
	pVectorGA = new double[Size];
	cout << "\nInput matrix:\n";
	srand(time(0));
	for (int i = 0; i < Size; i++)
	{
		for (int j = 0; j < Size; j++)
		{
			pMatrixGA[Size * i + j] = rand() % 100 + 50;
			pMatrixGA[Size * j + i] = pMatrixGA[Size * i + j];
			pMatrixGR[Size * i + j] = pMatrixGA[Size * i + j];
			pMatrixGR[Size * j + i] = pMatrixGA[Size * i + j];
		}
		pVectorGA[i] = rand() % 100 + 50;
		pVectorGR[i] = pVectorGA[i];
	}
	pPivotPos = new int[Size];
	pPivotIter = new int[Size];
	for (int i = 0; i < Size; ++i)
	{
		pPivotPos[i] = 0;
		pPivotIter[i] = -1;
	}

}
void ExactDataInitialization(double*& pVector, double*& pPrevProcResult, double*& pPrevProcD, double*& pPrevProcG,
	int* pSendInd, int Size, int RowNum, int ProcRank) {
	pPrevProcResult = new double[RowNum];
	pPrevProcD = new double[RowNum];
	pPrevProcG = new double[RowNum];
	for (int i = 0; i < RowNum; i++) {
		pPrevProcResult[i] = 0;
		pPrevProcD[i] = 0;
		pPrevProcG[i] = -pVector[pSendInd[ProcRank] / Size + i];
	}
}
void ProcessInitialization(double*& pMatrixGA, double*& pVectorGA, double*& pResultGA, double*& pMatrixGR, double*& pVectorGR, double*& pResultGR, double*& pProcRows, double*& pProcResult,
	int& Size, int& RowNum, int ProcRank, int ProcNum, double*& pPrevProcResult, double*& pPrevProcD, double*& pProcD, double*& pPrevProcG,
	double*& pProcG, double*& tmpVec, double& eps) {
	int RestRows; 

	if (ProcRank == 0) {
		do {
			cout << "\n\tEnter equations count: ";
			cin >> Size;
			if (Size < ProcNum) {
				cout << "\tCount of equations must be bigger than number of processes! \n ";
			}
		} while (Size < ProcNum);
		cout << "\tEnter a valid error value:"; 
		cin >> eps;
	}
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	RestRows = Size;
	for (int i = 0; i < ProcRank; i++)
		RestRows = RestRows - RestRows / (ProcNum - i);
	RowNum = RestRows / (ProcNum - ProcRank);

	pVectorGA = new double[Size];
	pResultGA = new double[Size];
	pVectorGR = new double[Size];
	pResultGR = new double[Size];
	pProcRows = new double[RowNum * Size];
	pProcResult = new double[RowNum];
	pProcD = new double[RowNum];
	pProcG = new double[RowNum];
	tmpVec = new double[RowNum];
	if (ProcRank == 0) {
		pMatrixGA = new double[Size * Size];
		pMatrixGR = new double[Size * Size];
		RandomDataInitialization(pMatrixGA, pVectorGA, pMatrixGR, pVectorGR, Size);
	}
}
void DataDistribution(double*& pMatrix, double*& pProcRows, double*& pVector, double*& pPrevProcResult, double*& pPrevProcD,
	double*& pPrevProcG, int*& pSendNum, int*& pSendInd, int Size, int RowNum, int ProcRank, int ProcNum) {
	int RestRows = Size;
	MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	pSendInd = new int[ProcNum];
	pSendNum = new int[ProcNum];
	RowNum = (Size / ProcNum);
	pSendNum[0] = RowNum * Size;
	pSendInd[0] = 0;
	for (int i = 1; i < ProcNum; i++) {
		RestRows -= RowNum;
		RowNum = RestRows / (ProcNum - i);
		pSendNum[i] = RowNum * Size;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
	}

	MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows,
		pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

	ExactDataInitialization(pVector, pPrevProcResult, pPrevProcD, pPrevProcG, pSendInd, Size, RowNum, ProcRank);
}
void ReceiveInfoCalculation(int*& pReceiveNum, int*& pReceiveInd, int Size, int ProcNum) {
	int i;
	int RestRows = Size; 
	pReceiveNum = new int[ProcNum];
	pReceiveInd = new int[ProcNum];

	pReceiveInd[0] = 0;
	pReceiveNum[0] = Size / ProcNum;
	for (i = 1; i < ProcNum; i++) {
		RestRows -= pReceiveNum[i - 1];
		pReceiveNum[i] = RestRows / (ProcNum - i);
		pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
	}
}
void ParallelMVmulCalculation(double* pProcRows, double* pVector, double* pProcResult, int Size, int RowNum) {
	int i, j;
	for (i = 0; i < RowNum; i++) {
		pProcResult[i] = 0;
		for (j = 0; j < Size; j++)
			pProcResult[i] += pProcRows[i * Size + j] * pVector[j];
	}
}
void ParallelVVsubCalculation(double* pProcLeft, double* pProcRight, double* pProcResult, int Size, int RowNum) {
	int i;
	for (i = 0; i < RowNum; i++) {
		pProcResult[i] = pProcLeft[i] - pProcRight[i];
	}
}
void ParallelVVsumCalculation(double* pProcLeft, double* pProcRight, double* pProcResult, int RowNum) {
	int i;
	for (i = 0; i < RowNum; i++) {
		pProcResult[i] = pProcLeft[i] + pProcRight[i];
	}
}
void ParallelVVmulCalculation(double* pProcLeft, double* pProcRight, double& Result, int RowNum) {
	int i;
	double ProcResult = 0;
	for (i = 0; i < RowNum; i++) {
		ProcResult += pProcLeft[i] * pProcRight[i];
	}
	MPI_Allreduce(&ProcResult, &Result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
void ParallelDVmulCalculation(double pProcLeft, double*& pProcRight, double*& ProcResult, int RowNum) {
	int i;
	for (i = 0; i < RowNum; i++) {
		ProcResult[i] = pProcLeft * pProcRight[i];
	}
}
void ParallelFirstIter(double*& pProcRows, double*& pVector, double*& pProcG, double*& pPrevProcResult, double*& tmpVec,
	int*& pReceiveNum, int*& pReceiveInd, int*& pSendInd, int ProcRank, int RowNum, int Size) {
	double* pPrevProcResultFull = new double[Size];
	MPI_Allgatherv(pPrevProcResult, pReceiveNum[ProcRank], MPI_DOUBLE, pPrevProcResultFull,
		pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
	ParallelMVmulCalculation(pProcRows, pPrevProcResultFull, tmpVec, Size, RowNum);
	ParallelVVsubCalculation(tmpVec, pVector + pSendInd[ProcRank] / Size, pProcG, Size, RowNum);
}
void ParallelSecondIter(double*& pProcG, double*& pPrevProcG, double*& pPrevProcResult, double*& pProcD, double*& pPrevProcD, double*& tmpVec,
	int*& pReceiveNum, int*& pReceiveInd, int*& pSendInd, int ProcRank, int RowNum, int Size) {
	double* parts = new double[2];
	ParallelVVmulCalculation(pProcG, pProcG, parts[0], RowNum);
	double* resparts = new double[2];
	ParallelVVmulCalculation(pPrevProcG, pPrevProcG, parts[1], RowNum);
	ParallelDVmulCalculation((double)(parts[0] / parts[1]), pPrevProcD, tmpVec, RowNum);

	ParallelVVsubCalculation(tmpVec, pProcG, pProcD, Size, RowNum);

}
void ParallelThirdIter(double*& pProcRows, double*& pProcG, double*& pProcD, double*& tmpVec,
	int*& pReceiveNum, int*& pReceiveInd, int*& pSendInd, double& s, int ProcRank, int RowNum, int Size) {
	double above, resabove;
	double* pProcDFull = new double[Size];
	ParallelVVmulCalculation(pProcD, pProcG, above, RowNum);
	MPI_Allgatherv(pProcD, pReceiveNum[ProcRank], MPI_DOUBLE, pProcDFull,
		pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
	ParallelMVmulCalculation(pProcRows, pProcDFull, tmpVec, Size, RowNum);
	double beyond;
	ParallelVVmulCalculation(pProcD, tmpVec, beyond, RowNum);
	s = -(double)(above / beyond);

}
void FourthIter(double*& pProcResult, double*& pPrevProcResult, double*& pProcD, double*& tmpVec,
	int*& pReceiveNum, int*& pReceiveInd, int*& pSendInd, double& s, int ProcRank, int RowNum, int Size) {
	ParallelDVmulCalculation(s, pProcD, tmpVec, RowNum);
	ParallelVVsumCalculation(tmpVec, pPrevProcResult, pProcResult, RowNum);

}


void ProcessTermination(double* pMatrixGA, double* pMatrixGR, double* pVectorGA, double* pVectorGR, double* pResultGA,
	double* pResultGR, double* pResultGRSeq, double* pPrevProcResult, double* pProcResult, double* pPrevResult, double* pPrevProcD,
	double* pProcD, double* pPrevD, double* pD, double* pPrevProcG, double* pProcG, double* pPrevG, double* pG, int ProcRank) {
	delete[] pPrevProcResult; delete[] pProcResult; delete[] pPrevResult; delete[] pPrevProcD;
	delete[] pProcD;  delete[] pPrevProcG; delete[] pProcG;
	if (ProcRank == 0) {
		delete[] pMatrixGR; delete[] pMatrixGA; delete[] pVectorGR; delete[] pVectorGA; delete[] pResultGA;
		delete[] pResultGR;
		delete[] pD; delete[] pG;
	}
}
void main(int argc, char* argv[]) {
	setlocale(LC_ALL, "Russian");
	double* pMatrixGA = 0; double* pMatrixGR = 0;
	double* pVectorGA = 0; double* pVectorGR = 0;
	double* pResultGA = 0; double* pResultGR = 0; double* pResultGRSeq = 0;
	double* pPrevProcResult = 0; double* pProcResult; double* pPrevResult = 0;
	double* pPrevProcD = 0; double* pProcD = 0; double* pPrevD = 0; double* pD = 0;
	double* pPrevProcG = 0; double* pProcG = 0;	double* pPrevG = 0; double* pG = 0;
	int* pSendNum; // Количество отправленных процессу элементов
	int* pSendInd; // Индекс первого среди них
	int* pReceiveNum; // Количество элементов, которые будет отправлять данный процесс
	int* pReceiveInd; // Индекс первого среди них
	double eps;//Погрешность
	double* tmpVec = 0;
	double s;//Значение, вычисляемое на шаге 3 каждой итерации
	int Size; // Размер матрицы и стобца свободных членов
	double* pProcRows;//Строки, выделенных данном процессу
	int RowNum;//Их количество
	int ProcRank, ProcNum;
	double Start, Finish, Duration, DurationPar;

	//------ Параллельная версия метода сопряженных градиентов ------//
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	ProcessInitialization(pMatrixGA, pVectorGA, pResultGA, pMatrixGR, pVectorGR, pResultGR, pProcRows, pProcResult, Size, RowNum, ProcRank, ProcNum,
		pPrevProcResult, pPrevProcD, pProcD, pPrevProcG, pProcG, tmpVec, eps);
	DataDistribution(pMatrixGR, pProcRows, pVectorGR, pPrevProcResult, pPrevProcD, pPrevProcG, pSendNum, pSendInd, Size, RowNum, ProcRank, ProcNum);
	ReceiveInfoCalculation(pReceiveNum, pReceiveInd, Size, ProcNum);
	int flag = 1;
	Start = MPI_Wtime();
	while(flag){
		double sumpogr;
		double sum = 0;
		ParallelFirstIter(pProcRows, pVectorGR, pProcG, pPrevProcResult, tmpVec, pReceiveNum, pReceiveInd, pSendInd, ProcRank, RowNum, Size);
		ParallelSecondIter(pProcG, pPrevProcG, pPrevProcResult, pProcD, pPrevProcD, tmpVec, pReceiveNum, pReceiveInd, pSendInd, ProcRank, RowNum, Size);
		ParallelThirdIter(pProcRows, pProcG, pProcD, tmpVec, pReceiveNum, pReceiveInd, pSendInd, s, ProcRank, RowNum, Size);
		FourthIter(pProcResult, pPrevProcResult, pProcD, tmpVec, pReceiveNum, pReceiveInd, pSendInd, s, ProcRank, RowNum, Size);
		for (int i = 0; i < RowNum; i++)
			sum += fabs(pProcResult[i] - pPrevProcResult[i]);
		MPI_Reduce(&sum, &sumpogr,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		if (ProcRank == 0) {
			if (sumpogr < eps) {
				flag = 0;
			}
		}
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		delete pPrevProcG; delete pPrevProcD; delete pPrevProcResult;
		pPrevProcG = pProcG; pPrevProcD = pProcD; pPrevProcResult = pProcResult;
		pProcG = new double[RowNum]; pProcD = new double[RowNum]; pProcResult = new double[RowNum];
	}
	MPI_Gatherv(pPrevProcResult, RowNum, MPI_DOUBLE, pResultGR,
		pReceiveNum, pReceiveInd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	Finish = MPI_Wtime();
	DurationPar = Finish - Start;
	//------------------------------------//
	if (ProcRank == 0) {
		
		Start = MPI_Wtime();
		SerialResultCalculation(pMatrixGA, pVectorGA, pResultGA, Size);
		Finish = MPI_Wtime();
		Duration = Finish - Start;
		cout << "\n\tThe result of the sequential Gauss method:\n";
		cout << "\tOperating time:" << Duration << " sec.\n";

		ProcessInitializationGradient(pVectorGR, pPrevResult, pResultGRSeq, pPrevD, pD, pPrevG, pG, Size);
		Start = MPI_Wtime();
		SoprGradSeq(pMatrixGR, pVectorGR, pPrevResult, pResultGRSeq, pPrevD, pD, pPrevG, pG, s, eps, Size);
		Finish = MPI_Wtime();
		Duration = Finish - Start;
		cout << "\n\tThe result of the sequential conjugate gradient method:\n";
		
		cout << "\tOperating time:" << Duration << " sec.\n";
		//--------------------------------------------------------------------------------------------//
		cout << "\n\tThe result of the parallel conjugate gradient method:\n";
		
		cout << "\tOperating time:" << DurationPar << " sec.\n";
		bool IsCorrect = 1;

		for (int i = 0; i < Size; i++) {
			if (abs(pResultGA[i] - pResultGR[i]) > eps)
				IsCorrect = 0;
		}
		if (IsCorrect) cout << "\n\tThe algorithm worked correctly.\n";
		else cout << "\n\tThere are errors in the algorithm.\n";
	}
	ProcessTermination(pMatrixGA, pMatrixGR, pVectorGA, pVectorGR, pResultGA, pResultGR, pResultGRSeq, pPrevProcResult, pProcResult, pPrevResult,
		pPrevProcD, pProcD, pPrevD, pD, pPrevProcG, pProcG, pPrevG, pG, ProcRank);
	MPI_Finalize();
}