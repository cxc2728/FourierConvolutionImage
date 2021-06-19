//  Project Title: corr2021
#define _CRT_SECURE_NO_WARNINGS

#include < iostream >
#include < fstream >
#include < string >
#include < io.h >
#include < dos.h >
#include < conio.h >
#include < stdlib.h >
#include < sstream >
#include < stdio.h >
#include < iomanip >
#include < istream >
#include < math.h >

using namespace std;

class corrCov2021 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data {

		double** Image_1; // pointer to the matrix entry 

		double** Image_2; // pointer to the matrix entry

		double** covariance;

		double** correlation;

	}*pointer; // pointer to the matrices

public:

	corrCov2021(int x, int y) : n1(x), n2(y) { };// constructor 

	void allocateData();

	~corrCov2021() { } // destructor

};

void corrCov2021::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	pointer = new data;

	pointer->Image_1 = new double* [this->n2];

	pointer->Image_2 = new double* [this->n2];

	pointer->covariance = new double* [this->n2];

	pointer->correlation = new double* [this->n2];


	for (int v = 0; v < this->n2; v++) { // (1)

		pointer->Image_1[v] = new double[this->n1];

		pointer->Image_2[v] = new double[this->n1];

		pointer->covariance[v] = new double[this->n1];

		pointer->correlation[v] = new double[this->n1];

	} // (1) allocate struct 'data' (end)


	  // (2) initialize (begin)
	for (int v = 0; v < this->n2; v++) { // (a)

		for (int f = 0; f < this->n1; f++) { // (b)

			pointer->Image_1[v][f] = (double)0.0;

			pointer->Image_2[v][f] = (double)0.0;

			pointer->covariance[v][f] = (double)0.0;

			pointer->correlation[v][f] = (double)0.0;

		} //(b)

	} //(a)
   // (2) initialize (end)

} // allocate data


int main(int argc, char* argv[]) {

	if (argc < 6) {
		std::cout << endl;
		std::cout << "Please type the 2 images file names" << endl;
		std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
		std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
		std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
		std::cout << "Please enter the size of the square window (integer in [2, 10])" << endl;
		std::cout << endl;
		exit(0);
	}

	else { // run the program (begin)

		char imageFileName_1[300];
		char imageFileName_2[300];

		sprintf(imageFileName_1, "%s", argv[1]);
		sprintf(imageFileName_2, "%s", argv[2]);

		int n1 = atoi(argv[3]);
		int n2 = atoi(argv[4]);

		int windowSize = atoi(argv[5]);

		if (windowSize < 2 || windowSize > 10)
		{

			std::cout << "The size of the square window must be in [2, 10] " << endl;
			exit(0);

		}

		std::cout << "The first image file name is: " << imageFileName_1 << endl;
		std::cout << "The second image file name is: " << imageFileName_2 << endl;
		std::cout << "The number of pixels along the X direction is: " << atoi(argv[3]) << endl;
		std::cout << "The number of pixels along the Y direction is: " << atoi(argv[4]) << endl;
		std::cout << "The size of the square window is: " << atoi(argv[5]) << endl;

		corrCov2021 computeCorrCov(n1, n2);

		computeCorrCov.allocateData();

		/// read image file (begin)
		FILE* pf;

		if ((pf = fopen(imageFileName_1, "rb+")) == NULL)
		{

			std::cout << "Cannot open file: " << imageFileName_1 << endl;
			exit(0);

		}
		else { // else

			double number;

			for (int i1 = 0; i1 < n2; i1++) {// x dim

				for (int i2 = 0; i2 < n1; i2++) { // y dim

					fread(&number, sizeof(double), 1, pf);

					computeCorrCov.pointer->Image_1[i1][i2] = (double)number;

				} // y dim

			}  // x dim 


			fclose(pf);


		} // else 
		/// read image file (end)

		std::cout << "First Image data loaded" << endl;

		/// read image file (begin)
		if ((pf = fopen(imageFileName_2, "rb+")) == NULL)
		{

			std::cout << "Cannot open file: " << imageFileName_2 << endl;
			exit(0);

		}
		else { // else

			double number;

			for (int i1 = 0; i1 < n2; i1++) {// x dim

				for (int i2 = 0; i2 < n1; i2++) { // y dim

					fread(&number, sizeof(double), 1, pf);

					computeCorrCov.pointer->Image_2[i1][i2] = (double)number;

				} // y dim

			}  // x dim 


			fclose(pf);


		} // else 
		/// read image file (end)

		std::cout << "Second Image data loaded" << endl;


		/// compute covariance & correlation between images (begin)
		double average_1 = 0, average_2 = 0;
		double stdev_1 = 0, stdev_2 = 0;
		long int counter = 0;
		double var_1 = 0, var_2 = 0;
		int i1, i2;

		for (i1 = 0; i1 < n2 - windowSize; i1++) {// x dim

			for (i2 = 0; i2 < n1 - windowSize; i2++) { // y dim

				average_1 = 0;
				average_2 = 0;

				counter = 0;

				stdev_1 = 0;
				stdev_2 = 0;

				for (int wx = i1; wx < i1 + windowSize; wx++) {// x dim

					for (int wy = i2; wy < i2 + windowSize; wy++) { // y dim

						counter++;

						average_1 += (double)computeCorrCov.pointer->Image_1[wx][wy];
						average_2 += (double)computeCorrCov.pointer->Image_2[wx][wy];

					
					} // y dim

				}  // x dim 

				average_1 /= (double)counter;
				average_2 /= (double)counter;

				for (int wx = i1; wx < i1 + windowSize; wx++) {// x dim

					for (int wy = i2; wy < i2 + windowSize; wy++) { // y dim

						
						stdev_1 += ((double)computeCorrCov.pointer->Image_1[wx][wy] - (double)average_1) *
						   	       ((double)computeCorrCov.pointer->Image_1[wx][wy] - (double)average_1);

						stdev_2 += ((double)computeCorrCov.pointer->Image_2[wx][wy] - (double)average_2) *
							       ((double)computeCorrCov.pointer->Image_2[wx][wy] - (double)average_2);

					} // y dim

				}  // x dim 

				var_1 = (double)sqrt(((double)stdev_1));
				var_2 = (double)sqrt(((double)stdev_2));

				for (int kx = i1; kx < i1 + windowSize; kx++) {// x dim

					for (int ky = i2; ky < i2 + windowSize; ky++) { // y dim

						computeCorrCov.pointer->covariance[i1][i2] +=

							((double)computeCorrCov.pointer->Image_1[kx][ky] - (double)average_1) *
							((double)computeCorrCov.pointer->Image_2[kx][ky] - (double)average_2);

						
					} // y dim

				}  // x dim 


				if (((double)var_1 * var_2) != 0.0)
				{

					computeCorrCov.pointer->correlation[i1][i2] =

						((double)computeCorrCov.pointer->covariance[i1][i2] / 
							
						((double)var_1 * var_2));

				}
				
			} // y dim

		}  // x dim 
		/// compute covariance & correlation between images (end)

		FILE* savedata;
		char outputFile[128];

		sprintf(outputFile, "%s", "covarianceImage.img");

		if ((savedata = fopen(outputFile, "wb+")) == NULL)
		{

			std::cout << "Cannot open output file, Now Exit..." << endl;

		}
		else { // (save)

			for (int i1 = 0; i1 < n2; i1++) {// x dim

				for (int i2 = 0; i2 < n1; i2++) { // y dim

					double save = (double)computeCorrCov.pointer->covariance[i1][i2];

					fwrite(&save, sizeof(double), 1, savedata);

				} // y dim

			}  // x dim 

			fclose(savedata);

		} // (save)


		sprintf(outputFile, "%s", "correlationImage.img");

		if ((savedata = fopen(outputFile, "wb+")) == NULL)
		{

			std::cout << "Cannot open output file, Now Exit..." << endl;

		}
		else { // (save)

			for (int i1 = 0; i1 < n2; i1++) {// x dim

				for (int i2 = 0; i2 < n1; i2++) { // y dim

					double save = (double)computeCorrCov.pointer->correlation[i1][i2];

					fwrite(&save, sizeof(double), 1, savedata);

				} // y dim

			}  // x dim 

			fclose(savedata);

		} // (save)


		std::cout << "Covariance & Correlation Images saved" << endl;
		std::cout << "End of Computation..." << endl;
		std::cout << endl;

		delete computeCorrCov.pointer;
		computeCorrCov.~corrCov2021();
	}  // run the program (end)

} // end of main