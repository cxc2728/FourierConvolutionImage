//  Project Title: BPF2D2018v2 (Band Pass Filter)
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

#define ft_SCALE_in 1

void OnZTransform(char imageFilename[], int rcxres, int rcyres, double m_Real, double m_Imaginary);
void OnInverseZTransformTransferFunction2021(char filename[], char BPF[], int rcxres, int rcyres, double m_Real, double m_Imaginary);

class BPF2D2018 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data {

		double** Image; // pointer to the matrix entry 

		double** BPF_2D; // pointer to the matrix entry

	}*pointer; // pointer to the matrices

public:

	BPF2D2018(int x, int y) : n1(x), n2(y) { };// constructor 

	void allocateData();

	void save();

	~BPF2D2018() { } // destructor

};

void BPF2D2018::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	pointer = new data;

	pointer->Image = new double* [this->n2];

	pointer->BPF_2D = new double* [this->n2];


	for (int v = 0; v < this->n2; v++) { // (1)

		pointer->Image[v] = new double[this->n1];

		pointer->BPF_2D[v] = new double[this->n1];

	} // (1) allocate struct 'data' (end)


	  // (2) initialize (begin)
	for (int v = 0; v < this->n2; v++) { // (a)

		for (int f = 0; f < this->n1; f++) { // (b)

			pointer->Image[v][f] = (double)0.0;

			pointer->BPF_2D[v][f] = (double)1.0;

		} //(b)

	} //(a)
   // (2) initialize (end)

} // allocate data


void BPF2D2018::save() { // saveImages

	FILE* savedata;
	char outputFile[128];

	sprintf(outputFile, "%s", "Image.img");

	if ((savedata = fopen(outputFile, "wb+")) == NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	}
	else { // (save)


		for (int v = 0; v < this->n2; v++) { // (a)

			for (int f = 0; f < this->n1; f++)

				fwrite(&pointer->Image[v][f], sizeof(double), 1, savedata);

		} // (a)

		fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s", "BPF_2D.img");
	if ((savedata = fopen(outputFile, "wb+")) == NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	}
	else { // (save)


		for (int v = 0; v < this->n2; v++) { // (a)

			for (int f = 0; f < this->n1; f++)

				fwrite(&pointer->BPF_2D[v][f], sizeof(double), 1, savedata);

		} // (a)

		fclose(savedata);

	} // (save)

} // saveImages


int main(int argc, char* argv[]) {

	char outputFile[128] = "BPF2D2018.log";

	FILE* savedata;

	if (argc < 7) {
		std::cout << endl;
		std::cout << "Please type the image file name" << endl;
		std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
		std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
		std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
		std::cout << "Please enter the cut-off frequency of the traditional BPF (double) > 0 and < 1" << endl;
		std::cout << "Please enter the value of the real part of the complex number (double) possibly in ]0, 1]" << endl;
		std::cout << "Please enter the value of the imaginary part of the complex number (double) possibly in [0, 1]" << endl;
		std::cout << endl;
		exit(0);
	}

	else { // run the program (begin)


		if ((savedata = fopen(outputFile, "w")) == NULL)
		{

			std::cout << "Cannot open output file, Now Exit..." << endl;

		}
		else { // processing (begin)

			int n1 = atoi(argv[2]);
			int n2 = atoi(argv[3]);

			char imageFileName[128];

			sprintf(imageFileName, "%s", argv[1]);

			double cutOFF_frequency = atof(argv[4]);

			double m_Real = atof(argv[5]);
			double m_Imaginary = atof(argv[6]);

			std::cout << endl;
			std::cout << "The image file name is: " << imageFileName << endl;
			std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
			std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
			std::cout << "The cutOFF frequency of the traditional BPF is: " << atof(argv[4]) << endl;
			std::cout << "The value of m_Real is: " << m_Real << endl;
			std::cout << "The value of m_Imaginary is: " << m_Imaginary << endl;


			fprintf(savedata, "%s%s\n", "The image file name is: ", imageFileName);
			fprintf(savedata, "%s%d\n", "The number of pixels along the X direction is: ", n1);
			fprintf(savedata, "%s%d\n", "The number of pixels along the Y direction is: ", n2);
			fprintf(savedata, "%s%lf\n", "The cutOFF frequency of the traditional BPF is: ", cutOFF_frequency);
			fprintf(savedata, "%s%f\n", "The value of m_Real is: ", m_Real);
			fprintf(savedata, "%s%f\n", "The value of m_Imaginary is: ", m_Imaginary);

			double pi = 3.141592;
			double Zphase = (double)2.0 * pi * atan2((double)m_Imaginary, (double)m_Real) / (n1 * n2);
			double magnitude = (double)sqrt((double)m_Real * m_Real + (double)m_Imaginary * m_Imaginary);

			std::cout << "The magnitude of the Complex number is: " << magnitude << endl;
			std::cout << "The phase of the Complex number is: " << Zphase << endl;

			BPF2D2018 BPF2D(n1, n2);

			BPF2D.allocateData();

			/// read image file (begin)
			FILE* pf;

			if ((pf = fopen(imageFileName, "rb+")) == NULL)
			{

				std::cout << "Cannot open file: " << imageFileName << endl;
				fprintf(savedata, "%s%s\n", "Cannot open file: ", imageFileName);
				exit(0);

			}
			else { // else

				double number;

				for (int i1 = 0; i1 < n2; i1++) {// x dim

					for (int i2 = 0; i2 < n1; i2++) { // y dim

						fread(&number, sizeof(double), 1, pf);

						BPF2D.pointer->Image[i1][i2] = (double)number;

					} // y dim

				}  // x dim 


				fclose(pf);


			} // else 
			/// read image file (end)

			std::cout << "Image data loaded" << endl;

			 
	int XNEI = (int)2;
    int n7 = ( (int)floor( (double)n2/2.0) );  
	int n8 = ( (int)floor( (double)n1/2.0) );
	double deltaT = 0.5; // 0.0
	double lowpassFilter = 0.0, highpassFilter = 0.0;
	double alpha = 0.0, RC = 0.5; // 0.0
	double x = 0.0, y = 0.0, w = 0.0;	
	
	
	// build the FILTER function www.wikipedia.org & www.dspguide.com/ch14/5.htm (begin)
	for (int pp =-n7+XNEI; pp < n7-XNEI; pp++) {

            for (int qq =-n8; qq < n8; qq++) {

				x = BPF2D.pointer->Image[pp + n7][qq + n8];

				y = BPF2D.pointer->BPF_2D[pp + n7 - 1][qq + n8];


				// low pass filter
				RC = (double) 1.0 / ((double) 2.0 * pi * cutOFF_frequency);

				deltaT = (double) 0.3; // 0.5

				if ( ((double) RC + deltaT) != 0.0 )

				alpha = ((double) deltaT) / ((double) RC + deltaT);

				else alpha = (double) 1.0;

				lowpassFilter = ( (double) alpha * x ) + ( ((double) 1.0-alpha ) * y );


				// high pass filter
				deltaT = (double) 1.0;

				if ( ((double) RC + deltaT) != 0.0 )

				alpha = ((double) RC) / ((double) RC + deltaT);
				
				else alpha = (double) 1.0;

				w = BPF2D.pointer->Image[pp + n7 - 1][qq + n8];

				highpassFilter = ((double) ((double) alpha * y ) + (((double)alpha) * ((double)x - w )) );
				

				// band pass filter
				// sum low pass and high pass
				BPF2D.pointer->BPF_2D[pp + n7][qq + n8] = ((double) lowpassFilter) + ((double) highpassFilter);
	
			}

     } // build the FILTER function www.wikipedia.org & www.dspguide.com/ch14/5.htm (end)

			std::cout << "Filter calculated" << endl;
			
			BPF2D.save();

			std::cout << "Now calculating the Direct Z Transform of the Signal..." << endl;
			OnZTransform(imageFileName, n1, n2, m_Real, m_Imaginary);

			std::cout << "Now calculating the Direct Z Transform of the BPF Signal..." << endl;
			OnZTransform("BPF_2D.img", n1, n2, m_Real, m_Imaginary);

			std::cout << "Now calculating the Transfer Function using Inverse Z Transform..." << endl;
			OnInverseZTransformTransferFunction2021(imageFileName, "BPF_2D.img", n1, n2, m_Real, m_Imaginary);


			std::cout << "End of Computation..." << endl;
			std::cout << endl;

			fprintf(savedata, "%s\n", "End of Computation...");
			fprintf(savedata, "\n");

			fclose(savedata);
			delete BPF2D.pointer;
			BPF2D.~BPF2D2018();
		} // processing (end)

	} // run the program (end)

	return 0;
} // end of main 

void OnZTransform(char imageFilename[], int rcxres, int rcyres, double m_Real, double m_Imaginary)
{

	int NofXpixels = rcxres;
	int NofYpixels = rcyres;

	int i, j, index;
	int dx, dy;
	int ds, dp;
	int k2, k3, w, t;

	double pi = 3.141592;

	double* ZSpaceR = 0;
	double* ZSpaceI = 0;
	double* Signal = 0;

	FILE* logfile;

	char logfilename[128] = "Z-T.log";

	if ((logfile = fopen(logfilename, "w+")) == NULL)
	{

		printf("%s\n %s\n", "Unable to open log File", "Now Exit");

		exit(0);

	}
	else { // allocate memory 


		if ((ZSpaceR = (double*)calloc(NofXpixels * NofYpixels, sizeof(double))) == NULL)
		{

			fprintf(logfile, "%s\n", "Not enough memory to allocate Real Image data: Exit");

			exit(0);

		}

		if ((ZSpaceI = (double*)calloc(NofXpixels * NofYpixels, sizeof(double))) == NULL)
		{

			fprintf(logfile, "%s\n", "Not enough memory to allocate Real Image data: Exit");

			// FIFO memory deallocation method
			free(ZSpaceR);
			exit(0);

		}

		if ((Signal = (double*)calloc(NofXpixels * NofYpixels, sizeof(double))) == NULL)
		{

			fprintf(logfile, "%s\n", "Not enough memory to allocate Real Image data: Exit");

			// FIFO memory deallocation method
			free(ZSpaceR);
			free(ZSpaceI);
			exit(0);

		}

	} // allocate memory 

	//// read image data and initialize pointers
	double number = 0.0;

	for (i = 0; i < NofYpixels; i++)
	{
		for (j = 0; j < NofXpixels; j++)
		{

			index = ((j * NofYpixels) + i);

			*(ZSpaceR + index) = (double)0.0;

			*(ZSpaceI + index) = (double)0.0;

		}

	}

	FILE* pf;
	char SignalFilename[128];
	double readData;

	sprintf(SignalFilename, "%s", imageFilename);

	if ((pf = fopen(SignalFilename, "rb+")) == NULL)
	{

		fprintf(logfile, "%s\n", "Cannot open file to read Signal");

		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(Signal);

		exit(0);

	}
	else { // read data


		for (i = 0; i < rcyres; i++)
		{ ///read signal data
			for (j = 0; j < rcxres; j++)
			{

				index = ((j * rcyres) + i);

				fread(&readData, sizeof(double), 1, pf);

				*(Signal + index) = (double)readData;

			}
		} ///read signal data

		fprintf(logfile, "%s\n", "Signal Read in DOUBLE (64bits) format");

		fclose(pf);
	} // save data

	// scale Signal (begin) 
	double max = *(Signal);
	double min= *(Signal);

	/// compute max and min of data (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if( *(Signal+index) > (double)max ) 
			
					max = (double) *(Signal+index);
              
				if( *(Signal+index) < (double)min ) 
			
					min = (double) *(Signal+index);
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if ( max == min ) *(Signal+index) = (double)0.0;

				else { 
					
				*(Signal+index) = (double) fabs( (min - (double)*(Signal+index) ) / (min - max) );

				*(Signal+index) *= ft_SCALE_in;

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)


	double phase;
	double exponent, r;
	double complexRZ, complexIZ;

	double magnitude = (double)sqrt((double)m_Real * m_Real + (double)m_Imaginary * m_Imaginary);
	double Zphase = (double)2.0 * pi * atan2((double)m_Imaginary, (double)m_Real) / ((double)NofXpixels * NofYpixels);

	///// Z Transform //////
	for (i = 0; i < NofYpixels; i++)
	{ ///calculate Z-Space data

		for (j = 0; j < NofXpixels; j++)
		{


			dx = ((int)i - NofYpixels / 2);
			dy = ((int)j - NofXpixels / 2);

			k2 = ((int)(dy * NofYpixels) + dx);

			w = ((j * NofYpixels) + i);

			for (int s = 0; s < NofYpixels; s++)
			{ ///calculate Z-Space data 
				for (int p = 0; p < NofXpixels; p++)
				{


					ds = ((int)s - NofYpixels / 2);
					dp = ((int)p - NofXpixels / 2);

					k3 = ((int)(ds * NofXpixels) + dp);

					t = ((p * NofYpixels) + s);


					phase = ((double)(2.0 * pi * k2 * k3) / (NofXpixels * NofYpixels));

					exponent = (double)2.0 * pi * t * (double)Zphase / ((double)pow((double)NofXpixels * NofYpixels, 2.0));

					exponent = (double)fabs((double)exponent);

					r = (double)pow((double)magnitude, -(double)exponent);


					complexRZ = (double)cos((double)phase) + (double)sin((double)phase);

					complexIZ = -(double)sin((double)phase) + (double)cos((double)phase);


					*(ZSpaceR + w) += (double)*(Signal + t) * (double)complexRZ * (double)r;

					*(ZSpaceI + w) += (double)*(Signal + t) * (double)complexIZ * (double)r;
				}

			}///calculate Z-Space data 


		}
	} ///calculate Z-Space data
	///// Z Transform //////

	double savedata = 0.0;
	char Zfilename[128];

	sprintf(Zfilename, "%s%s", "Z-SpaceR-", imageFilename);

	fprintf(logfile, "%s\t%s\n", "Now Saving Z-Space Signal (Real) in File: ", Zfilename);

	if ((pf = fopen(Zfilename, "wb+")) == NULL)
	{

		fprintf(logfile, "%s\n", "Cannot open file to save Z-Space Signal");


		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(Signal);

		exit(0);

	}
	else { // save data


		for (i = 0; i < NofYpixels; i++)
		{ ///save Z-Space data
			for (j = 0; j < NofXpixels; j++)
			{

				index = ((j * NofYpixels) + i);

				savedata = (double)*(ZSpaceR + index);

				fwrite(&savedata, sizeof(double), 1, pf);

			}
		} ///save Z-Space data

		fprintf(logfile, "%s\n", "Z-Space Signal (Real) Saved");

		fclose(pf);
	} // save data



	sprintf(Zfilename, "%s%s", "Z-SpaceI-", imageFilename);

	fprintf(logfile, "%s\t%s\n", "Now Saving Z-Space Signal (Imaginary) in File: ", Zfilename);

	if ((pf = fopen(Zfilename, "wb+")) == NULL)
	{

		fprintf(logfile, "%s\n", "Cannot open file to save Z-Space Signal");

		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(Signal);

		exit(0);

	}
	else { // save data


		for (i = 0; i < NofYpixels; i++)
		{ ///save Z-Space data
			for (j = 0; j < NofXpixels; j++)
			{

				index = ((j * NofYpixels) + i);

				savedata = (double)*(ZSpaceI + index);

				fwrite(&savedata, sizeof(double), 1, pf);

			}
		} ///save Z-Space data

		fprintf(logfile, "%s\n", "Z-Space Signal (Imaginary) Saved");

		fclose(pf);

	} // save data

	sprintf(Zfilename, "%s%s", "Z-SpaceM-", imageFilename);

	fprintf_s(logfile, "%s\t%s\n", "Now Saving Z-Space Magnitude of the Signal in File: ", Zfilename);

	if ((pf = fopen(Zfilename, "wb+")) == NULL)
	{

		fprintf_s(logfile, "%s\n", "Cannot open file to save Z-Space Magnitude of the Signal");

		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(Signal);

		exit(0);

	}
	else { // save data	

	 // save a zero image (begin)
		for (int s = 0; s < NofYpixels; s++)
		{
			for (int p = 0; p < NofXpixels; p++)
			{

				savedata = (double)0.0;

				fwrite(&savedata, sizeof(double), 1, pf);

			}
		} // save a zero image (end)

		fclose(pf);

	}

	if ((pf = fopen(Zfilename, "wb+")) == NULL)
	{

		fprintf_s(logfile, "%s\n", "Cannot open file to save Z-Space Magnitude of the Signal");

		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(Signal);

		exit(0);

	}
	else { // save data

	 // Z-Space Magnitude (begin)
		for (int s = 0; s < (int)NofYpixels; s++)
		{
			for (int p = 0; p < (int)NofXpixels; p++)
			{


				index = ((p * NofYpixels) + s);

				savedata = (double)sqrt((double)*(ZSpaceR + index) * (double)*(ZSpaceR + index) +
					                    (double)*(ZSpaceI + index) * (double)*(ZSpaceI + index));

				fwrite(&savedata, sizeof(double), 1, pf);

			}
		} // Z-Space Magnitude (end)

		fprintf_s(logfile, "%s\n", "Z-Space Magnitude of the Signal Saved");

		fclose(pf);
	} // save data

	printf("%s\n", "Z Processing Completed");
	fprintf_s(logfile, "%s\n", "Z Processing Completed");

	fclose(logfile);


	// FIFO memory deallocation method
	free(ZSpaceR);
	free(ZSpaceI);
	free(Signal);

}

void OnInverseZTransformTransferFunction2021(char filename[], char BPF[], int rcxres, int rcyres, double m_Real, double m_Imaginary)
{

	int NofXpixels = rcxres;
	int NofYpixels = rcyres;

	int i, j, index;
	int dx, dy;
	int ds, dp;
	int k2, k3, w, t;

	double pi = 3.141592;

	double phase;

	//2010
	double emittingSource = 1.4145; // 2021
	double scale = ((double)rcxres * rcyres * emittingSource);
	//2010

	FILE* logfile;
	char logfilename[128] = "INV-ZT.log";

	FILE* image;
	char imageFilename[256];

	double* ZSpaceR = 0;
	double* ZSpaceI = 0;
	double* ZSpaceR_BPF = 0;
	double* ZSpaceI_BPF = 0;
	double* reconSignal = 0;

	if ((logfile = fopen(logfilename, "w+")) == NULL)
	{

		exit(0);

	}
	else { // allocate memory


		printf("%s\n", "Now INV Z Processing...");
		fprintf(logfile, "%s\n", "Now INV Z Processing...");

		if ((ZSpaceR = (double*)calloc(NofXpixels * NofYpixels, sizeof(double))) == NULL)
		{

			fprintf(logfile, "%s\n", "Not enough memory to allocate Real Image data: Exit");

			exit(0);

		}

		if ((ZSpaceI = (double*)calloc(NofXpixels * NofYpixels, sizeof(double))) == NULL)
		{

			fprintf(logfile, "%s\n", "Not enough memory to allocate Real Image data: Exit");

			// FIFO memory deallocation method
			free(ZSpaceR);
			exit(0);

		}


		if ((reconSignal = (double*)calloc(NofXpixels * NofYpixels, sizeof(double))) == NULL)
		{

			fprintf(logfile, "%s\n", "Not enough memory to allocate Imaginary Image data: Exit");

			// FIFO memory deallocation method
			free(ZSpaceR);
			free(ZSpaceI);

			exit(0);

		}


		if ((ZSpaceR_BPF = (double*)calloc(NofXpixels * NofYpixels, sizeof(double))) == NULL)
		{

			fprintf(logfile, "%s\n", "Not enough memory to allocate Real Image data: Exit");

			// FIFO memory deallocation method
			free(ZSpaceR);
			free(ZSpaceI);
			free(reconSignal);

			exit(0);

		}

		if ((ZSpaceI_BPF = (double*)calloc(NofXpixels * NofYpixels, sizeof(double))) == NULL)
		{

			fprintf(logfile, "%s\n", "Not enough memory to allocate Real Image data: Exit");

			// FIFO memory deallocation method
			free(ZSpaceR);
			free(ZSpaceI);
			free(reconSignal);
			free(ZSpaceR_BPF);

			exit(0);

		}

	} // allocate memory

	/// Read Z-Space of Signal (begin)	
	//// read image data and initialize pointers
	sprintf(imageFilename, "%s%s", "Z-SpaceR-", filename);

	if ((image = fopen(imageFilename, "rb+")) == NULL)
	{

		fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);

		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(reconSignal);
		free(ZSpaceR_BPF);
		free(ZSpaceI_BPF);

		exit(0);

	}
	else { // read data and initialize pointers

		double number = 0.0;

		for (i = 0; i < NofYpixels; i++)
		{
			for (j = 0; j < NofXpixels; j++)
			{

				index = ((j * NofYpixels) + i);

				fread(&number, sizeof(double), 1, image);

				*(ZSpaceR + index) = (double)number;


			}

		}

		fclose(image);

	}// read data and initialize pointers


	char imageFilename2[128];

	sprintf(imageFilename2, "%s%s", "Z-SpaceI-", filename);


	if ((image = fopen(imageFilename2, "rb+")) == NULL)
	{

		fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(reconSignal);
		free(ZSpaceR_BPF);
		free(ZSpaceI_BPF);

		exit(0);

	}
	else { // read data and initialize pointers

		double number = 0.0;

		for (i = 0; i < NofYpixels; i++)
		{
			for (j = 0; j < NofXpixels; j++)
			{

				index = ((j * NofYpixels) + i);

				fread(&number, sizeof(double), 1, image);

				*(ZSpaceI + index) = (double)number;

			}

		}

		fclose(image);


		for (i = 0; i < NofYpixels; i++)
		{
			for (j = 0; j < NofXpixels; j++)
			{

				index = ((j * NofYpixels) + i);

				*(reconSignal + index) = (double)0.0;

			}

		}


	}// read data and initialize pointers
	/// Read Z-Space of Signal (end)	


	/// Read Z-Space of BPF Signal (begin)	
	//// read image data and initialize pointers
	sprintf(imageFilename, "%s%s", "Z-SpaceR-", BPF);

	if ((image = fopen(imageFilename, "rb+")) == NULL)
	{

		fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);

		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(reconSignal);
		free(ZSpaceR_BPF);
		free(ZSpaceI_BPF);

		exit(0);

	}
	else { // read data and initialize pointers

		double number = 0.0;

		for (i = 0; i < NofYpixels; i++)
		{
			for (j = 0; j < NofXpixels; j++)
			{

				index = ((j * NofYpixels) + i);

				fread(&number, sizeof(double), 1, image);

				*(ZSpaceR_BPF + index) = (double)number;


			}

		}

		fclose(image);

	}// read data and initialize pointers


	sprintf(imageFilename2, "%s%s", "Z-SpaceI-", BPF);

	if ((image = fopen(imageFilename2, "rb+")) == NULL)
	{

		fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(reconSignal);
		free(ZSpaceR_BPF);
		free(ZSpaceI_BPF);

		exit(0);

	}
	else { // read data and initialize pointers

		double number = 0.0;

		for (i = 0; i < NofYpixels; i++)
		{
			for (j = 0; j < NofXpixels; j++)
			{

				index = ((j * NofYpixels) + i);

				fread(&number, sizeof(double), 1, image);

				*(ZSpaceI_BPF + index) = (double)number;

			}

		}

		fclose(image);

	}// read data and initialize pointers
	/// Read Z-Space of BPF Signal (end)	


	// scale Signal (begin) 
	double max = *(ZSpaceR);
	double min= *(ZSpaceR);

	/// compute max and min of data (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if( *(ZSpaceR+index) > (double)max ) 
			
					max = (double) *(ZSpaceR+index);
              
				if( *(ZSpaceR+index) < (double)min ) 
			
					min = (double) *(ZSpaceR+index);
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if ( max == min ) *(ZSpaceR+index) = (double)0.0;

				else { 
					
				*(ZSpaceR+index) = (double) fabs( (min - (double)*(ZSpaceR+index) ) / (min - max) );

				*(ZSpaceR+index) *= ft_SCALE_in;

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)

	// scale Signal (begin) 
	max = *(ZSpaceI);
	min= *(ZSpaceI);

	/// compute max and min of data (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if( *(ZSpaceI+index) > (double)max ) 
			
					max = (double) *(ZSpaceI+index);
              
				if( *(ZSpaceI+index) < (double)min ) 
			
					min = (double) *(ZSpaceI+index);
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if ( max == min ) *(ZSpaceI+index) = (double)0.0;

				else { 
					
				*(ZSpaceI+index) = (double) fabs( (min - (double)*(ZSpaceI+index) ) / (min - max) );

				*(ZSpaceI+index) *= ft_SCALE_in;

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)


	// scale Signal (begin) 
	max = *(ZSpaceR_BPF);
	min= *(ZSpaceR_BPF);

	/// compute max and min of data (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if( *(ZSpaceR_BPF+index) > (double)max ) 
			
					max = (double) *(ZSpaceR_BPF+index);
              
				if( *(ZSpaceR_BPF+index) < (double)min ) 
			
					min = (double) *(ZSpaceR_BPF+index);
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if ( max == min ) *(ZSpaceR_BPF+index) = (double)0.0;

				else { 
					
				*(ZSpaceR_BPF+index) = (double) fabs( (min - (double)*(ZSpaceR_BPF+index) ) / (min - max) );

				*(ZSpaceR_BPF+index) *= ft_SCALE_in;

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)

	// scale Signal (begin) 
	max = *(ZSpaceI_BPF);
	min= *(ZSpaceI_BPF);

	/// compute max and min of data (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if( *(ZSpaceI_BPF+index) > (double)max ) 
			
					max = (double) *(ZSpaceI_BPF+index);
              
				if( *(ZSpaceI_BPF+index) < (double)min ) 
			
					min = (double) *(ZSpaceI_BPF+index);
		
		} // y dim
        
	}  // x dim
	/// compute max and min of data (end)

	// scale (begin)
	for (i=0; i<rcyres; i++)
	{
		for (j=0; j<rcxres; j++)
		{

				index = ((j*rcyres)+i);

				if ( max == min ) *(ZSpaceI_BPF+index) = (double)0.0;

				else { 
					
				*(ZSpaceI_BPF+index) = (double) fabs( (min - (double)*(ZSpaceI_BPF+index) ) / (min - max) );

				*(ZSpaceI_BPF+index) *= ft_SCALE_in;

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)

	double real = 0.0, imaginary = 0.0;
	double real_BPF = 0.0, imaginary_BPF = 0.0;
	double real_tf = 0.0, imaginary_tf = 0.0;
	double exponent, r;

	double magnitude = (double)sqrt((double)m_Real * m_Real + (double)m_Imaginary * m_Imaginary);
	double Zphase = (double)2.0 * pi * atan2((double)m_Imaginary, (double)m_Real) / ((double)NofXpixels * NofYpixels);

	///// INV Z Transform //////
	for (i = 0; i < NofYpixels; i++)
	{ ///process Z-Space data

		for (j = 0; j < NofXpixels; j++)
		{

			dx = ((int)i - NofYpixels / 2);
			dy = ((int)j - NofXpixels / 2);

			k2 = ((int)(dx * NofXpixels) + dy);

			w = ((j * NofYpixels) + i);

			real = (double)0.0;
			imaginary = (double)0.0;

			real_BPF = (double)0.0;
			imaginary_BPF = (double)0.0;


			for (int s = 0; s < NofYpixels; s++)
			{ ///process Z-Space data

				for (int p = 0; p < NofXpixels; p++)
				{

					ds = ((int)s - NofYpixels / 2);
					dp = ((int)p - NofXpixels / 2);

					k3 = ((int)(dp * NofYpixels) + ds);

					t = ((p * NofYpixels) + s);


					phase = ((double)(2.0 * pi * k2 * k3) / (NofXpixels * NofYpixels));

					exponent = (double)2.0 * pi * t * (double)Zphase / ((double)pow((double)NofXpixels * NofYpixels, 2.0));

					exponent = (double)fabs((double)exponent);

					r = (double)pow((double)magnitude, (double)exponent);


					real += ((double)*(ZSpaceR + t) * (double)cos((double)phase)) -
						    ((double)*(ZSpaceI + t) * (double)sin((double)phase));

					real_BPF += ((double)*(ZSpaceR_BPF + t) * (double)cos((double)phase)) -
						        ((double)*(ZSpaceI_BPF + t) * (double)sin((double)phase));

					real *= (double)r;

					real_BPF *= (double)r;


					imaginary += ((double)*(ZSpaceR + t) * (double)sin((double)phase)) +
						         ((double)*(ZSpaceI + t) * (double)cos((double)phase));

					imaginary_BPF += ((double)*(ZSpaceR_BPF + t) * (double)sin((double)phase)) +
						             ((double)*(ZSpaceI_BPF + t) * (double)cos((double)phase));


					imaginary *= (double)r;

					imaginary_BPF *= (double)r;

				}

			}///process Z-Space data 


			if (((double)real) != 0.0) {

				real_tf = ((double)real_BPF / real);

				if ((_isnan((double)real_tf)) == 0) {}

				else { real_tf = (double)1.0; }

			}
			else { real_tf = (double)1.0; }


			if (((double)imaginary) != 0.0) {

				imaginary_tf = ((double)imaginary_BPF / imaginary);

				if ((_isnan((double)imaginary_tf)) == 0) {}

				else { imaginary_tf = (double)1.0; }

			}
			else { imaginary_tf = (double)1.0; }


			*(reconSignal + w) = (double)sqrt(((double)real_tf * real_tf) + ((double)imaginary_tf * imaginary_tf));

			*(reconSignal + w) /= (double)scale;

		}
	} ///process Z-Space data


	double savedata = 0.0;
	FILE* pf;
	char reconFilename[128];

	sprintf(reconFilename, "%s%s", "reconSignal-", filename);


	fprintf(logfile, "%s\t%s\n", "Now Saving Reconstructed Signal in File: ", reconFilename);

	if ((pf = fopen(reconFilename, "wb+")) == NULL)
	{

		fprintf(logfile, "%s\n", "Cannot open file to save Z-Space Signal");

		// FIFO memory deallocation method
		free(ZSpaceR);
		free(ZSpaceI);
		free(reconSignal);
		free(ZSpaceR_BPF);
		free(ZSpaceI_BPF);

		exit(0);

	}
	else { // save data


		for (i = 0; i < NofYpixels; i++)
		{ ///save Z-Space data
			for (j = 0; j < NofXpixels; j++)
			{

				index = ((j * NofYpixels) + i);

				savedata = (double)*(reconSignal + index);

				fwrite(&savedata, sizeof(double), 1, pf);

			}
		} ///save Z-Space data

		fprintf(logfile, "%s\n", "Reconstructed Signal Saved");

		fclose(pf);
	} // save data


	printf("%s\n", "Inverse Z Processing Completed");
	fprintf(logfile, "%s\n", "Inverse Z Processing Completed");

	fclose(logfile);


	// FIFO memory deallocation method
	free(ZSpaceR);
	free(ZSpaceI);
	free(reconSignal);
	free(ZSpaceR_BPF);
	free(ZSpaceI_BPF);

}

