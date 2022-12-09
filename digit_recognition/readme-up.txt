// HMM_3.cpp : Defines the entry point for the console application.
//

//

#include "stdafx.h"
#include "math.h"
#include "float.h"
#include "string.h"
#include "limits.h"
#include "windows.h"
const int N = 5;
const int M = 32;
const int MAXT = 160;
int T = 0;
#define p 12
long double A[N + 1][N + 1] = {0};
long double Bini[N + 1][M + 1] = {0};
long double Aini[N + 1][N + 1] = {0};
long double B[N + 1][M + 1] = {0};
long double Abar[N + 1][N + 1] = {0};
long double Bbar[N + 1][M + 1] = {0};
long double Amatrix[N + 1][N + 1] = {0}, Bmatrix[N + 1][M + 1] = {0};

long double tokhuraWt[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
long double const threshold = 1e-30;
int frames = 1;
int O[MAXT + 1] = {0};
long double PI[N + 1] = {0};
long double PImatrix[N + 1] = {0};
long double PIbar[N + 1] = {0};

long double alpha[MAXT + 1][N + 1] = {0};
long double beta[MAXT + 1][N + 1] = {0};
long double gamma[MAXT + 1][N + 1] = {0};
long double delta[MAXT + 1][N + 1] = {0};
long double psi[MAXT + 1][N + 1] = {0};
long double Xi[MAXT + 1][N + 1][N + 1] = {0};
int Q[MAXT + 1] = {0};
long double pstar = 0;
long double codebook[32][12];

/// @brief Function to read Initial Model
void readvalues()
{
	FILE *ptr = fopen("HMM_AIJ_FINAL.txt", "r"); // reading AIJ values

	while (!feof(ptr))
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				fscanf(ptr, "%lf", &Aini[i][j]);
			}
		}
	}
	fclose(ptr);
	FILE *fptr = fopen("HMM_BJK_FINAL.txt", "r"); // reading BIJ file
	while (!feof(fptr))
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= M; j++)
			{
				fscanf(fptr, "%lf", &Bini[i][j]);
			}
			printf("\n");
		}
	}

	fclose(fptr);

	FILE *tr = fopen("PI.txt", "r"); // reading PI values
	while (!feof(tr))
	{
		for (int i = 1; i <= N; i++)
		{
			fscanf(tr, "%lf", &PI[i]);
		}
	}
	fclose(tr);
}

/// @brief Function to assign initial model to Matrix A and B
void read_initial_model()
{

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			A[i][j] = Aini[i][j];
		}
	}
	/*	for(int i=1;i<=N;i++){
			for(int j=1;j<=N;j++){
			   printf("%.32e ",A[i][j]);
			}
			printf("\n");
		}*/
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			B[i][j] = Bini[i][j];
		}
	}
	/*	for(int i=1;i<=N;i++){
			for(int j=1;j<=M;j++){
			   printf("%.32e ",B[i][j]);
			}
			printf("\n");
		}*/
}

/// @brief Function to read the Codebook
void readcodebook()
{
	FILE *fptr = fopen("codebook-chi.txt", "r");
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 12; j++)
		{
			fscanf(fptr, "%lf", &codebook[i][j]);
		}
	}
}


/// @brief Function to Calculate tokhura distance
/// @param currentregion 
/// @param currentuniverse 
/// @return tokhuradistane
long double tokhuradistance(long double currentregion[], long double currentuniverse[])
{
	long double ans = 0.0;
	for (int i = 1; i <= p; i++)
	{
		//  printf("%lf ",currentregion[i]);
	}
	for (int i = 0; i < p; i++)
	{
		//   printf("%lf ",currentuniverse[i]);
	}
	for (int i = 0; i < p; i++)
	{
		long double d = currentregion[i + 1] - currentuniverse[i]; /// calculating tokhura distance b/e universe code vector and current cluster
		ans += (tokhuraWt[i] * d * d);
	}
	// printf("\n%lf",ans);
	return ans;
}

/// @brief Function to Calculate Ci and generating observation sequence using tokhura distance
/// @param R 
void durbins(long double R[])
{
	long double alpha[p + 1][p + 1];
	long double K[p + 1];
	long double sum;
	long double E[p + 1];
	E[0] = R[0];
	FILE *ptr;

	for (int i = 1; i <= p; i++)
	{
		sum = 0;
		for (int j = 1; j <= i - 1; j++)
		{
			sum += alpha[j][i - 1] * R[i - j];
		}
		K[i] = (R[i] - sum) / E[i - 1];
		alpha[i][i] = K[i];
		for (int j = 1; j <= i - 1; j++)
		{
			alpha[j][i] = alpha[j][i - 1] - K[i] * alpha[i - j][i - 1]; // calculating Ai's
		}
		E[i] = (1 - K[i] * K[i]) * E[i - 1];
	}

	long double c[p + 1];
	long double sumc;
	c[1] = alpha[1][p];
	for (int i = 2; i <= p; i++)
	{
		sumc = 0.0;
		for (int k = 1; k <= i - 1; k++)
		{

			sumc += ((double)k / i) * c[k] * alpha[i - k][p];
		}
		c[i] = alpha[i][p] + sumc; // calculating ci's
	}
	for (int m = 1; m <= p; m++)
	{
		c[m] *= 1 + ((long double)p / 2) * sin(((long double)22 / 7 * m) / p); // applyng raised sine window
	}

	long double mindistance = 111111;
	int index = 1;

	long double currentcodebook[12];
	for (int i = 0; i < 32; i++)
	{
		for (int k = 0; k < 12; k++)
		{

			currentcodebook[k] = codebook[i][k];
			// printf("%lf ",currentcodebook[k]);
		}
		long double distance = tokhuradistance(c, currentcodebook);
		//  printf("fds =\n%lf ",distance);
		if (distance < mindistance)
		{
			mindistance = distance; // findingthe minimun distance
			index = i + 1;
		}
	}
	// printf("%d",index);
	O[frames++] = index;
	// ci++;
	// fprintf(ptr, "\n");
	// fclose(ptr);
}


/// @brief Fucntion to calculate Soltion 1
/// @return probility P(O/lambda)
long double Forwardprocedure()                    //change the return type to void during training
{ // change to void during testing
	/*for(int i=1;i<=N;i++){
		for(int j=1;j<=N;j++){
		  printf("%.32e ",A[i][j]);
		}
		printf("\n");
	}
	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++){
		  printf("%.32e ",B[i][j]);
		}
		printf("\n");
	}*/
	// step-1 :initialization
	for (int i = 1; i <= N; i++)
	{
		alpha[1][i] = PI[i] * B[i][O[1]];
	}

	/*	for(int i=1;i<=N;i++){
			printf("\n%.32e",alpha[1][i]);
		}*/
	// step-2 Induction
	//	printf("%d\n",T);
	//	printf("%d\n",O[T]);
	for (int t = 1; t <= T - 1; t++)
	{
		for (int j = 1; j <= N; j++)
		{
			long double sum = 0;
			for (int i = 1; i <= N; i++)
			{
				sum += alpha[t][i] * A[i][j];
			}
			alpha[t + 1][j] = sum * B[j][O[t + 1]];
		}
	}

	long double probability = 0;
	/*	for(int i=1;i<=T;i++){
			for(int j=1;j<=N;j++){
		   printf("%.32e ",alpha[i][j]);
			}
			printf("\n");
		}*/
	for (int i = 1; i <= N; i++)
	{
		probability += alpha[T][i];
	}
		printf("\nProbability of (O/lambda) = %.32e\n",probability);     //Probability of P(O/lambda)
	//	fprintf(fptr,"\nProbability of (O/lambda) = %.32e\n",probability);
	// fprintf(fpte,"\n");
	// step-1 initialization
	for (int i = 1; i <= N; i++)
	{
		beta[T][i] = 1;
	}

	// step-2 Induction
	for (int t = T - 1; t > 0; t--)
	{
		for (int i = 1; i <= N; i++)
		{

			long double sum = 0;
			for (int j = 1; j <= N; j++)
			{
				sum += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
			}
			beta[t][i] = sum;
		}
	}
	/*for(int i=1;i<=T;i++){
		for(int j=1;j<=N;j++){
	   printf("%.32e",beta[i][j]);
		}
		printf("\n");
	}*/
	return probability;
}

/// @brief Function to calulate gamma values
void gamma_values()
{

	for (int t = 1; t <= T; t++)
	{
		for (int j = 1; j <= N; j++)
		{
			long double sum = 0;
			for (int i = 1; i <= N; i++)
			{
				sum += alpha[t][i] * beta[t][i];
			}
			gamma[t][j] = alpha[t][j] * beta[t][j] / sum; // caluclating gamma values
		}
	}
	int stateseq[MAXT + 1] = {0};
	for (int t = 1; t <= T; t++)
	{
		long double max = DBL_MIN;
		int index = 0;
		for (int i = 1; i <= N; i++)
		{
			if (gamma[t][i] > max)
			{
				max = gamma[t][i];
				index = i;
			}
		}
		stateseq[t] = index; // finding most likely states
	}
}

/// @brief Function to calculate Solution 2 
void viterbi()
{
	// Viterbi Algorithm
	// step-1 Initialization
	for (int i = 1; i <= N; i++)
	{
		delta[1][i] = PI[i] * B[i][O[1]];
		psi[1][i] = 0;
	}

	// step-2 Recursion
	for (int t = 2; t <= T; t++)
	{
		for (int j = 1; j <= N; j++)
		{
			long double max = 0, ti = 0;
			int index = 0;
			for (int i = 1; i <= N; i++)
			{

				ti = delta[t - 1][i] * A[i][j];
				if (ti > max)
				{
					max = ti;
					index = i;
				}
			}
			delta[t][j] = max * B[j][O[t]];
			psi[t][j] = index;
		}
	}

	// step-3 Termination
	long double max = 0;
	for (int i = 1; i <= N; i++)
	{
		if (delta[T][i] > max)
		{
			max = delta[T][i];
			Q[T] = i;
		}
		pstar = max;
	}
	for (int t = T - 1; t > 0; t--)
	{
		Q[t] = psi[t + 1][Q[t + 1]];
	}
	printf("\n pstar= %.32\n", pstar);
}

/// @brief Function to calculate Xi values
void Xi_values()
{
	long double summation[MAXT + 1];
	for (int t = 1; t <= T - 1; t++)
	{

		summation[t] = 0;
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				summation[t] += alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
			}
		}

		for (int i = 1; i <= N; i++)
		{
			long double x;
			for (int j = 1; j <= N; j++)
			{
				x = alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
				Xi[t][i][j] = x / summation[t];
			}
		}
	}
}


/// @brief Function to restimate Model again (Solution to 3)
void restimating_values()
{
	// reevaluating PI
	for (int i = 1; i <= N; i++)
	{
		PIbar[i] = gamma[1][i];
	}

	// Re evaluating A
	for (int i = 1; i <= N; i++)
	{
		long double max = 0, sum = 0;
		int ind_i = 0, ind_j = 0;

		for (int j = 1; j <= N; j++)
		{
			long double t1 = 0, t2 = 0;
			for (int t = 1; t <= T - 1; t++)
			{
				t1 += Xi[t][i][j];
				t2 += gamma[t][i];
			}

			Abar[i][j] = t1 / t2;
			sum += Abar[i][j];
			if (Abar[i][j] > max)
			{
				max = Abar[i][j];
				ind_i = i;
				ind_j = j;
			}
		}
		if (sum > 1)
			Abar[ind_i][ind_j] -= (sum - 1);
		if (sum < 1)
			Abar[ind_i][ind_j] += (1 - sum);
	}

	// Reevaluating B

	long double sum1 = 0, sum2 = 0;
	for (int j = 1; j <= N; j++)
	{
		int count = 0;
		long double max = 0;
		int ind_j = 0, ind_k = 0;

		for (int k = 1; k <= M; k++)
		{
			sum1 = 0, sum2 = 0;
			for (int t = 1; t <= T; t++)
			{
				sum1 = sum1 + gamma[t][j];
				// if(Q[t]==j){
				if (O[t] == k)
				{
					sum2 = sum2 + gamma[t][j];
					//}
				}
			}
			Bbar[j][k] = sum2 / sum1;

			// finding max
			if (Bbar[j][k] > max)
			{
				max = Bbar[j][k];
				ind_j = j;
				ind_k = k;
			}

			// updating new bij with threshold value if it is zero
			if (Bbar[j][k] == 0)
			{
				Bbar[j][k] = threshold;
				count++;
			}
		}
		Bbar[ind_j][ind_k] = max - count * threshold;
	}
}


/// @brief Function to update model Abar and Bbar
void updating_model()
{

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			A[i][j] = Abar[i][j];
			//  Amatrix[i][j]+=Abar[i][j];
		}
	}
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			B[i][j] = Bbar[i][j];
			//    Bmatrix[i][j]+=Bbar[i][j];
		}
	}
	for (int i = 1; i <= N; i++)
	{
		//	PImatrix[i]+=PIbar[i];
		//	PI[i]=PIbar[i];
	}
}


/// @brief Function to store A matrix values
/// @param digit 
void write_A_values(int digit)
{
	char fn[5000];
	char fno[5000];
	sprintf(fno, "%d", digit);
	strcpy(fn, "");
	strcat(fn, "outputA");
	strcat(fn, fno);
	strcat(fn, ".txt");

	printf("%s", fn);
	FILE *fptr = fopen(fn, "a");
	// FILE* fptr=fopen(fn,"a");
	if (fptr == NULL)
	{
		printf("A me loach");
	}
	printf("\nPrinting updated A matrix's\n");
	//	fprintf(fptr,"Printing updated Amatrix's");
	//	fprintf(fptr,"\n");
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			//  printf("%.32e ",Abar[i][j]);
			// printing updated AIJ
			Amatrix[i][j] += Abar[i][j];
			fprintf(fptr, "%.32e ", Abar[i][j]);
		}
		fprintf(fptr, "\n");
		printf("\n");
	}
	// fprintf(fptr,"\n");
	fclose(fptr);
}

/// @brief Functiont to store B matrix values
/// @param digit 
void write_B_values(int digit)
{
	char fn[500];
	char fno[500];
	sprintf(fno, "%d", digit);
	strcpy(fn, "");
	strcat(fn, "outputB");
	strcat(fn, fno);
	strcat(fn, ".txt");

	printf("%s", fn);
	FILE *fptr = fopen(fn, "a");
	if (fptr == NULL)
	{
		printf("B me loach");
	}

	// 	FILE* fptr=fopen("outputB.txt","a");
	printf("\nPrinting updated B matrix's\n");
	// fprintf(fptr,"Printing updated B matrix's");
	//	fprintf(fptr,"\n");
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			//   printf("%.32e ",Bbar[i][j]);
			// printing updated BJK
			Bmatrix[i][j] += Bbar[i][j];
			fprintf(fptr, "%.32e ", Bbar[i][j]);
		}
		fprintf(fptr, "\n");
		printf("\n");
	}
	// fprintf(fptr,"\n");
	fclose(fptr);
}

/// @brief Function to store PI matrix values
/// @param digit 
void write_PI_values(int digit)
{
	char fn[5000];
	char fno[5000];
	sprintf(fno, "%d", digit);
	strcpy(fn, "");
	strcat(fn, "outputPI");
	strcat(fn, fno);
	strcat(fn, ".txt");

	printf("%s", fn);
	FILE *fptr = fopen(fn, "a");
	if (fptr == NULL)
	{
		printf("PI me loach");
	}

	// FILE* fptr=fopen("outputPI.txt","a");
	for (int i = 1; i <= N; i++)
	{
		// printf("%.32e ",PIbar[i]);
		// printing updated PI
		PImatrix[i] += PIbar[i];
		fprintf(fptr, "%.32e ", PIbar[i]);
	}
	fprintf(fptr, "\n");
	fclose(fptr);
}
/// @brief Function to store Final A matrix value of specified digit
/// @param digit 
void write_A_values_final(int digit)
{
	char fn[5000];
	char fno[5000];
	sprintf(fno, "%d", digit);
	strcpy(fn, "");
	strcat(fn, "outputAFinal");
	strcat(fn, fno);
	strcat(fn, ".txt");

	printf("%s", fn);
	FILE *fptr = fopen(fn, "a");
	printf("\nPrinting updated A matrix's Final\n");
	//	fprintf(fptr,"Printing updated Amatrix's");
	//	fprintf(fptr,"\n");
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			//  printf("%.32e ",Amatrix[i][j]);                            //printing updated AIJ
			fprintf(fptr, "%.32e ", Amatrix[i][j]);
		}
		fprintf(fptr, "\n");
		printf("\n");
	}
	// fprintf(fptr,"\n");
	fclose(fptr);
}
/// @brief Function to store Final B matrix value of specified digit
/// @param digit 
void write_B_values_final(int digit)
{
	char fn[500];
	char fno[500];
	sprintf(fno, "%d", digit);
	strcpy(fn, "");
	strcat(fn, "outputBFinal");
	strcat(fn, fno);
	strcat(fn, ".txt");

	printf("%s", fn);
	FILE *fptr = fopen(fn, "a");

	// 	FILE* fptr=fopen("outputB.txt","a");
	printf("\nPrinting updated B matrix's Final\n");
	// fprintf(fptr,"Printing updated B matrix's");
	//	fprintf(fptr,"\n");
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			//  printf("%.32e ",Bmatrix[i][j]);                              //printing updated BJK
			fprintf(fptr, "%.32e ", Bmatrix[i][j]);
		}
		fprintf(fptr, "\n");
		printf("\n");
	}
	// fprintf(fptr,"\n");
	fclose(fptr);
}
/// @brief Function to store Final PI matrix value of specified digit
/// @param digit 
void write_PI_values_final(int digit)
{
	char fn[100];
	char fno[100];
	sprintf(fno, "%d", digit);
	strcpy(fn, "");
	strcat(fn, "outputPIFinal");
	strcat(fn, fno);
	strcat(fn, ".txt");

	printf("%s", fn);
	FILE *fptr = fopen(fn, "a");

	// FILE* fptr=fopen("outputPI.txt","a");
	for (int i = 1; i <= N; i++)
	{
		printf("%.32e ", PImatrix[i]); // printing updated PI
		fprintf(fptr, "%.32e ", PImatrix[i]);
	}
	fprintf(fptr, "\n");
	fclose(fptr);
}
/// @brief Function to calculate CI and observation values
/// @param filename 
void training_values(char *filename)
{
	FILE *ptr, *fp, *read; // Declaring  the file pointer to read txt
	frames = 1;
	ptr = fopen("silence.txt", "r"); // Reading silence file to remove DcShift
	long double maxvalue = 1;		 // Declaring variable for maxvalue , sum,count for no of data in file
	long double sum = 0;
	long double *CISarray;

	long double count = 0;

	while (!feof(ptr))
	{ // Reading file till EOF
		long double lineno;

		fscanf(ptr, "%lf", &lineno);

		sum += lineno; // Sum of all data
		count++;
	}

	long double avg = sum / count; // calculating the mean of all data

	fclose(ptr);
	read = fopen(filename, "r");
	// reading input file - change  filename to yes.txt for yes input file
	if (read == NULL)
	{
		printf("cant read the fie");
	}
	long double counter = 0;
	char ignore[5000];
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	while (!feof(read))
	{
		long double lineno;

		fscanf(read, "%lf", &lineno);
		counter++;
		if (abs(maxvalue) < abs(lineno))
			maxvalue = lineno; // Calculating the max value in the data file
	}
	fseek(read, 0, SEEK_SET);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);

	long double arr[50000];
	int i = 0;
	long double energyarray[100];
	while (!feof(read))
	{
		long double lineno;
		fscanf(read, "%lf", &lineno);
		lineno = lineno - avg;				   // DCshift
		lineno = ((lineno * 5000) / maxvalue); // Normalizing the data
		arr[i] = lineno;
		i++;
	}
	int framecounter = 0;
	long double energy = 0.0;
	int value = 0;
	int window = 0;
	int step = 1;
	for (int k = 0; k <= counter; k++) // iterating for all data in text File and calculating energy of each 320 frames in file
	{
		long double w = 0.54 - 0.46 * cos((2 * 3.1428 * window) / 319); // calculating hamming window
		arr[k] = arr[k] * w;

		if (framecounter == 320)
		{

			// energyarray[value] = energy / 320.0;    //storing the energyofeach frame in an array
			// energy = 0.0;
			framecounter = 1;
			window = 0;
			long double R[p + 1];
			long double sums;
			for (int i = 0; i <= p; i++)
			{
				// R[i]=0.0;
				sums = 0.0;
				for (int j = k - 320; j < k - i; j++)
				{
					sums += arr[j] * arr[j + i];
				}
				R[i] = sums; // calculating Ri's
			}
			durbins(R);
			k = 80 * step;
			step++;
		}
		if (frames >= 161)
			break;
		// energy += arr[k] * arr[k];
		framecounter++;
		window++;
	}

	// printf("\nenergyindex is=%d",energyindex);
	T = frames - 1;
	// printf("\nframes is =%d\n",frames);
	 for(int i=1;i<=frames-1;i++){

	 printf("%d ",O[i]);
	 }
	////framecounter++;
	// printf("\n");
}
/*void testing_values(char *filename)
{
	FILE *ptr, *fp, *read; // Declaring  the file pointer to read txt
	frames = 1;
	ptr = fopen("silence.txt", "r"); // Reading silence file to remove DcShift
	long double maxvalue = 1;		 // Declaring variable for maxvalue , sum,count for no of data in file
	long double sum = 0;
	long double *CISarray;

	long double count = 0;

	while (!feof(ptr))
	{ // Reading file till EOF
		long double lineno;

		fscanf(ptr, "%lf", &lineno);

		sum += lineno; // Sum of all data
		count++;
	}

	long double avg = sum / count; // calculating the mean of all data

	fclose(ptr);
	read = fopen(filename, "r");
	// reading input file - change  filename to yes.txt for yes input file
	if (read == NULL)
	{
		printf("cant read the fie");
	}
	long double counter = 0;
	char ignore[5000];
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	while (!feof(read))
	{
		long double lineno;

		fscanf(read, "%lf", &lineno);
		counter++;
		if (abs(maxvalue) < abs(lineno))
			maxvalue = lineno; // Calculating the max value in the data file
	}
	fseek(read, 0, SEEK_SET);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);
	fgets(ignore, sizeof(ignore), read);

	long double arr[50000];
	int i = 0;
	long double energyarray[100];
	while (!feof(read))
	{
		long double lineno;
		fscanf(read, "%lf", &lineno);
		lineno = lineno - avg;				   // DCshift
		lineno = ((lineno * 5000) / maxvalue); // Normalizing the data
		arr[i] = lineno;
		i++;
	}
	int framecounter = 0;
	long double energy = 0.0;
	int value = 0;
	int window = 0;
	int step = 1;
	long double energy = 0.0;
	
	
	for (int k = 0; k <= counter; k++) // iterating for all data in text File and calculating energy of each 320 frames in file
	{
		long double w = 0.54 - 0.46 * cos((2 * 3.1428 * window) / 319); // calculating hamming window
		arr[k] = arr[k] * w;

		if (framecounter == 320)
		{

			// energyarray[value] = energy / 320.0;    //storing the energyofeach frame in an array
			// energy = 0.0;
			framecounter = 1;
			window = 0;
			long double R[p + 1];
			long double sums;
			for (int i = 0; i <= p; i++)
			{
				// R[i]=0.0;
				sums = 0.0;
				for (int j = k - 320; j < k - i; j++)
				{
					sums += arr[j] * arr[j + i];
				}
				R[i] = sums; // calculating Ri's
			}
			durbins(R);
			k = 80 * step;
			step++;
		}
		if (frames >= 161)
			break;
		// energy += arr[k] * arr[k];
		framecounter++;
		window++;
	}

	// printf("\nenergyindex is=%d",energyindex);
	T = frames - 1;
	// printf("\nframes is =%d\n",frames);
	 for(int i=1;i<=frames-1;i++){

	 printf("%d ",O[i]);
	 }
	////framecounter++;
	// printf("\n");
}*/



/// @brief Function to average model of specified digit
void average_model()
{

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			Amatrix[i][j] /= 25;
		}
	}
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			Bmatrix[i][j] /= 25;
		}
	}
	for (int i = 1; i <= N; i++)
	{
		PImatrix[i] /= 25;
	}
}


/// @brief Function to erase all models and update initial model
/// @param digit 
void updating_initial_and_erasing_all_models(int digit)
{
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			Aini[i][j] = Amatrix[i][j];
		}
	}

	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			Bini[i][j] = Bmatrix[i][j];
		}
	}
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			Amatrix[i][j] = 0;
		}
	}
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			A[i][j] = 0;
		}
	}
	for (int i = 1; i <= N; i++)
	{
		PImatrix[i] = 0;
	}
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			Bmatrix[i][j] = 0;
		}
	}
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j <= M; j++)
		{
			B[i][j] = 0;
		}
	}

	char fn[100];
	char fno[100];
	sprintf(fno, "%d", digit);
	strcpy(fn, "");
	strcat(fn, "outputB");
	strcat(fn, fno);
	strcat(fn, ".txt");
	FILE *ptr;
	fclose(fopen(fn, "w"));
	char fn1[100];
	char fno1[100];
	sprintf(fno1, "%d", digit);
	strcpy(fn1, "");
	strcat(fn1, "outputA");
	strcat(fn1, fno1);
	strcat(fn1, ".txt");

	fclose(fopen(fn1, "w"));
	char fn2[50];
	char fno2[50];
	sprintf(fno2, "%d", digit);
	strcpy(fn2, "");
	strcat(fn2, "outputPI");
	strcat(fn2, fno2);
	strcat(fn2, ".txt");

	fclose(fopen(fn2, "w"));
}


/// @brief Function to read generated models of each digit
/// @param digit 
void read_models(int digit)
{
	char fn[5000];
	char fno[5000];
	sprintf(fno, "%d", digit);
	strcpy(fn, "");
	strcat(fn, "outputAFinal");
	strcat(fn, fno);
	strcat(fn, ".txt");

	//  printf("%s",fn);
	FILE *fptr = fopen(fn, "r");
	// reading AIJ values

	while (!feof(fptr))
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				fscanf(fptr, "%lf", &A[i][j]);
			}
		}
	}
	fclose(fptr);
	// reading BIJ file
	char fn1[5000];
	char fno1[5000];
	sprintf(fno1, "%d", digit);
	strcpy(fn1, "");
	strcat(fn1, "outputBFinal");
	strcat(fn1, fno1);
	strcat(fn1, ".txt");

	// printf("%s",fn1);
	fptr = fopen(fn1, "r");
	while (!feof(fptr))
	{
		for (int i = 1; i <= N; i++)
		{
			for (int j = 1; j <= M; j++)
			{
				fscanf(fptr, "%lf", &B[i][j]);
			}
			// printf("\n");
		}
	}
	fclose(fptr);
	FILE *tr = fopen("PI.txt", "r"); // reading PI values
	while (!feof(tr))
	{
		for (int i = 1; i <= N; i++)
		{
			fscanf(tr, "%lf", &PI[i]);
		}
	}
	fclose(tr);
}

/// @brief Function to test the data i.e. each utterance of different digit
void testing()
{
	char fno[5000];
	char fno1[5000];
	char fn[5000];
	readcodebook();
	int probCount = 0;
	for (int j = 0; j <= 9; j++)
	{
		for (int i = 26; i <= 30; i++)
		{
			/*sprintf(fno, "%d", i);
			sprintf(fno1, "%d", j);
			strcpy(fn, "");
			strcat(fn, "new-data/");
			strcat(fn, fno1);
			strcat(fn, "/224101056_");
			strcat(fn, "E_");
			strcat(fn, fno1);
			strcat(fn, "_");
			strcat(fn, fno);
			strcat(fn, ".txt");*/
			sprintf(fno, "%d", i);
				sprintf(fno1, "%d", j);
				strcpy(fn, "");
			//	strcat(fn, "new-data/");
				strcat(fn, fno1);
				strcat(fn, "/224101063_");
				strcat(fn, fno1);
				strcat(fn, "_");
				
				strcat(fn, "E_");
				strcat(fn, fno);
				strcat(fn, ".txt");
				printf("%s", fn);
				
			printf("digit = %d utterance no=%d ", j, i);
			training_values(fn);
			long double maxprobability = INT_MIN;
			long double probability = 0.0;
			int ind = 0;
			for (int i = 0; i <= 9; i++)
			{
				read_models(i);
				probability = Forwardprocedure();
				//	printf("%.32e",probability);
				if (probability > maxprobability)
				{
					maxprobability = probability;
					ind = i;
				}
			}
			if (ind == j)
				probCount++;
			printf("\nFile is of digit=%d\n", ind);
		}
	}

	printf("\n\n count = %d and Accuracy = %f\n", probCount, probCount/0.5);
}


/// @brief Function to train the model
void training()
{

	char fno[5000];
	char fno1[5000];
	char fn[5000];

	readcodebook();

	int count = 1;
	for (int d = 0; d <= 9; d++)
	{
		//  if(d==6)
		readvalues();
		count = 1;
	//	while (count <= 2)
		//{
			for (int i = 1; i <= 25; i++)
			{
				sprintf(fno, "%d", i);
				sprintf(fno1, "%d", d);
				strcpy(fn, "");
			//	strcat(fn, "new-data/");
				strcat(fn, fno1);
				strcat(fn, "/224101063_");
				strcat(fn, fno1);
				strcat(fn, "_");
				
				strcat(fn, "E_");
				strcat(fn, fno);
				strcat(fn, ".txt");
				printf("%s", fn);
				read_initial_model();

				training_values(fn);
				int iterations = 1;
				while (iterations <= 100)
				{
					iterations++;
					Forwardprocedure();
					gamma_values();
					viterbi();
					Xi_values();
					restimating_values();
					updating_model();
				}
				// write_PI_values(d);
				write_A_values(d);
				write_B_values(d);
			}
			average_model();
			//count++;
			//if (count == 3)
			//{
				write_A_values_final(d);
				write_B_values_final(d);
				//  write_PI_values_final(d);
			//}
			updating_initial_and_erasing_all_models(d);
		}
	//}
}

//@brief Functin to Live test the model
void Livetesting(){
	char inputFile[100];
	char File[100];
	readcodebook();
	system("Recording_Module.exe 2 livetesting.wav livetesting.txt");
	sprintf(inputFile, "livetesting.txt");
//	printf("%s", inputFile);
	int counter=0;
	FILE* read=fopen(inputFile,"r");
	FILE *fptr=fopen("livetesting1.txt","w");
	while (!feof(read) && counter<4500)
	{
		long double lineno;

		fscanf(read, "%lf", &lineno);
		counter++;
		 // Calculating the max value in the data file
	}
	counter=0;
	while(!feof(read) && counter<13000){
		int data;
		fscanf(read,"%d",&data);
		fprintf(fptr,"%d\n",data);
		counter++;
	 
	}
	fclose(fptr);
	fclose(read);
	sprintf(File,"livetesting1.txt");
	printf("%s",File);
	training_values(inputFile);
			long double maxprobability = INT_MIN;
			long double probability = 0.0;
			int ind = 0;
			for (int i = 0; i <= 9; i++)
			{
				read_models(i);
				probability = Forwardprocedure();
				//	printf("%.32e",probability);
				if (probability > maxprobability)
				{
					maxprobability = probability;
					ind = i;
				}
			}
			printf("\nFile is of digit=%d\n", ind);
		}
	



int _tmain(int argc, _TCHAR *argv[])

{

	// training();     //Function to train the model
	printf("\n");
	int ch;
	while(1){
		printf("\n Press 1 for manual testing\n 2 For live testing\n 3 Exit\n");
		scanf("%d",&ch);
	 switch(ch){
	 case 1: testing();
		      break;
	 case 2: Livetesting();
		      break;
	 case 3: exit(EXIT_FAILURE);
		      break;
	 }
	}

  getchar();
	return 0; 
}
