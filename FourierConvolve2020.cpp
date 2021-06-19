#define _CRT_SECURE_NO_WARNINGS
// name of the project: FourierConvolve2020
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

#define SCALE_in 255
#define ft_SCALE_in 1

void OnFourierTransform(char imageFilename[], int rcxres, int rcyres);
void OnInverseFourierTransformConvolve(char fodfilename[], char icfFilename[], int rcyres, int rcxres);
void runKurtosis(char imageFileName [], int n1, int n2, FILE * savedata);
void runSRE2D(char imageFileName [], int xd, int yd, FILE * savedata);
void runFOD2D(char imageFileName [], int xd, int yd, FILE * savedata);
void runFOD2D_nopad(char imageFileName [], int xd, int yd, FILE * savedata);


class Kurtosis2D2020 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data { 

		double **Signal; // pointer to the matrix entry

		double **Kurtosis;

		double **average;

		double **m4;

		double **m2;

	}*pointer; // pointer to the matrices

public:

	Kurtosis2D2020(int x, int y) : n1(x), n2(y) { };// constructor 
	
	void allocateData();

	void save();

	~Kurtosis2D2020() { } // destructor

};

void Kurtosis2D2020::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;

	 pointer->Signal = new double*[this->n2];

	 pointer->Kurtosis = new double*[this->n2];

	 pointer->average = new double*[this->n2];

	 pointer->m4 = new double*[this->n2];

	 pointer->m2 = new double*[this->n2];


	 for( int v=0; v < this->n2; v++ ) { // (1)

		 pointer->Signal[v] = new double[this->n1];

		 pointer->Kurtosis[v] = new double[this->n1];

		 pointer->average[v] = new double[this->n1];

		 pointer->m4[v] = new double[this->n1];

		 pointer->m2[v] = new double[this->n1];


	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for( int v=0; v < this->n2; v++ ) { // (a)

			for( int f=0; f < this->n1 ; f++ ) { // (b)

			pointer->Signal[v][f] = (double)0.0;

			pointer->Kurtosis[v][f] = (double)0.0;

			pointer->average[v][f] = (double)0.0;

			pointer->m4[v][f] = (double)0.0;

			pointer->m2[v][f] = (double)0.0;

			 } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data

class SRE2D2013 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data {

		double **fMRI; // pointer to the matrix entry 

		double **ICF; // pointer to the matrix entry

	}*pointer; // pointer to the matrices

public:

	SRE2D2013(int x, int y) : n1(x), n2(y) { };// constructor 
	
	void allocateData();

	~SRE2D2013() { } // destructor

};

void SRE2D2013::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;
			
	 pointer->fMRI = new double*[this->n2];

	 pointer->ICF = new double*[this->n2];
	

	 for( int v=0; v < this->n2; v++ ) { // (1)
		 
		 pointer->fMRI[v] = new double[this->n1];

		 pointer->ICF[v] = new double[this->n1];
	
	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for(int v=0; v < this->n2; v++ ) { // (a)

			for( int f=0; f < this->n1 ; f++ ) { // (b)
		 
			pointer->fMRI[v][f] = (double)0.0;

			pointer->ICF[v][f] = (double)0.0;

			 } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


class FOD2D2018 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data { 

		double **Signal; // pointer to the matrix entry

		double **fodSRE2D;

	}*pointer; // pointer to the matrices

public:

	FOD2D2018(int x, int y) : n1(x), n2(y) { };// constructor 
	
	void allocateData();

	~FOD2D2018() { } // destructor

};


void FOD2D2018::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;

	 pointer->Signal = new double*[this->n2];

	 pointer->fodSRE2D = new double*[this->n2];


	 for( int v=0; v < this->n2; v++ ) { // (1)

		 pointer->Signal[v] = new double[this->n1];

		 pointer->fodSRE2D[v] = new double[this->n1];

	  } // (1) allocate struct 'data' (end)


		// (2) initialize (begin)
		for( int v=0; v < this->n2; v++ ) { // (a)

			for( int f=0; f < this->n1 ; f++ ) { // (b)

			pointer->Signal[v][f] = (double)0.0;

			pointer->fodSRE2D[v][f] = (double)0.0;

			 } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


// read the input parameters from the
// console command line
int main ( int argc, char * argv[] ) { 

	// assign to the char string 'outputFile'
	// the value "FourierConvolve2020.log"
	char outputFile[128] = "FourierConvolve2020.log";

	// declare a pointer to a file by the
	// name 'savedata'
	FILE * savedata;

// tell the user of the list of input parameters necessary to tun the program:
if (argc < 4) { std::cout << endl;
				 std::cout << "Please type the MR image file name" << endl;
				 std::cout << "Please make sure that the images format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << endl;
				 exit(0); }

else { // run the program (begin)
	
	// opens the log file which name is 
	// contained in 'outputFile', opens 
	// in write mode
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{
		// alert the user of possible failure in opening the file
		std::cout << "Cannot open output file, Now Exit..." << endl;
		exit(0);

	} else  { // processing (begin)

	// declare an array of char to contain a string
	char imageFileName[128];

	// transfer into the array 'imageFileName'
	// the image file name as per input from
	// the command console. The image file name
	// is 'argv[1]'
	sprintf(imageFileName, "%s", argv[1]);

	// reads from the command console the
	// value of the image size (number of rows
	// and number of columns)
	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);

	// inform the user of the image size 
	// (number of rows and number of columns
	// of the matrix containing the image)
	std::cout << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;

	// save into the log file the image size 
	// (number of rows and number of columns
	// of the matrix containing the image)
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);

	/// processing (begins)
	// run FOD of the SRE2D model polynomial function 
	runFOD2D(imageFileName, n1, n2, savedata);

	// calculate the SRE2D ICF
	runSRE2D(imageFileName, n1, n2, savedata);

	// Fourier transform the FOD2D of the image
	std::cout << "Direct Fourier Transforming the FOD of the Image..." << endl;
	char fodFile[200];
	sprintf(fodFile, "%s%s","fodSRE2D-", imageFileName);
	OnFourierTransform(fodFile, n1, n2);
	
	// Fourier transform the ICF of the image
	std::cout << "Direct Fourier Transforming the ICF of the Image..." << endl;
	char icfFile[200];
	sprintf(icfFile, "%s%s", "ICF_SRE2D-" , argv[1]);
	OnFourierTransform(icfFile, n1, n2);


	OnInverseFourierTransformConvolve(fodFile, icfFile, n1, n2);

	
	// run FOD of the reconstructed signal
	runFOD2D_nopad("reconSignal.img", n1, n2, savedata);

	// run Kurtosis of the reconstructed signal
	runKurtosis("reconSignal.img", n1, n2, savedata);
	
	// run kurtosis of the Signal
	runKurtosis(imageFileName, n1, n2, savedata);

	// run kurtosis of the FOD of the signal
	runKurtosis(fodFile, n1, n2, savedata);

	// run kurtosis of the ICF of the signal
	runKurtosis(icfFile, n1, n2, savedata);
	/// processing (ends)

	// run kurtosis of the FOD of the reconstructed signal
	runKurtosis("fodSRE2D-reconSignal.img", n1, n2, savedata);

	// alert the user of the end of the program
	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	// save to log file
	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	fclose(savedata); 
	} // processing (end)

	} // run the program (end)

	// ANSI C requires the 'main'
	// function returning a value: 
	// zero in this case
	return 0;
} // end of main 


void runSRE2D(char imageFileName [], int xd, int yd, FILE * savedata)
{

	int n1 = xd;
	int n2 = yd;

	double XPixelSize = 1.0;
	double YPixelSize = 1.0;

	double x_misplacement_X = 0.5;
	double y_misplacement_Y = 0.5;

	double theta = 0.0;
	
	printf("%s\n", "Calculation of the SRE2D ICF");
	fprintf(savedata,"%s\n", "Calculation of the SRE2D ICF");
	
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n2);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n1);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The XY rotation angle is: ", theta);

    double misplacement_X = ((double)1.0 - ( cos( (double)theta ) + sin( (double)theta ) ) + x_misplacement_X);
    double misplacement_Y = ((double)1.0 - ( -sin( (double)theta ) + cos( (double)theta ) ) + y_misplacement_Y);

      misplacement_X = ((double)misplacement_X/XPixelSize);
      misplacement_Y = ((double)misplacement_Y/YPixelSize);

	  //////////////////***********//////////////////////
	  // Above formula scales the misplacement to the  //
	  // pixel size the same way the following formula //
	  // would do: (min - misplacement)/(min - max)    //  
	  //////////////////***********//////////////////////

	int PAD = 1;

	SRE2D2013 SRE(n1+PAD*2,n2+PAD*2);

	SRE.allocateData();

	/// read image file (begin)
	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		exit(0);

	} else { // else

	double number;

	
		for( int v=PAD+1; v < n2+PAD+1; v++ ) { // (a)

				for( int f=PAD+1; f < n1+PAD+1; f++ ) { // (b)

				fread(&number,sizeof(double),1,pf);
		
				SRE.pointer->fMRI[v][f] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl;


	// scale Signal (begin) 
	double max = -SRE.pointer->fMRI[n2/2][n1/2];
	double min = SRE.pointer->fMRI[n2/2][n1/2];


	for( int v=PAD+1; v < n2+PAD+1; v++ ) { // (a)

	    for( int f=PAD+1; f < n1+PAD+1; f++ ) { // (b)
		

			if( SRE.pointer->fMRI[v][f] > (double)max ) 
			
					max = (double) SRE.pointer->fMRI[v][f];
              
			if( SRE.pointer->fMRI[v][f] < (double)min ) 
			
					min = (double) SRE.pointer->fMRI[v][f];
	
		}
	}
		
		// scale (begin)
		for( int v=PAD+1; v < n2+PAD+1; v++ ) { // (a)

				for( int f=PAD+1; f < n1+PAD+1; f++ ) { // (b)


				if ( max == min ) SRE.pointer->fMRI[v][f] = (double)0.0;

				else { 
					
				SRE.pointer->fMRI[v][f] = (double) fabs( (min - (double)SRE.pointer->fMRI[v][f] ) / (min - max) );

				SRE.pointer->fMRI[v][f] *= SCALE_in;

				// the absolute vale is needed otherwise the minimum is -0 which is nonsense in computers. The value -0 was discovered
				// printing to the screen. Mathematically, the scaling formula gives values in between [0, SCALE] and not [-0, SCALE].
				// Therefore the nonsense is removed taking the absolute value.

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)
	
	
	// calculate ICF (begin)
	std::cout << "Compute Intensity-Curvature Functional SRE2D" << endl;

	double k=0, s=0;

	for( int i1=PAD+2; i1 < n2+PAD; i1++ ) { // (a)

	    for( int i2=PAD+2; i2 < n1+PAD; i2++ ) { // (b)

       k = (double) SRE.pointer->fMRI[i1][i2] * misplacement_X * misplacement_Y + 
				   ( misplacement_X * misplacement_Y * misplacement_X / 2.0) * 
				   ( SRE.pointer->fMRI[i1+1][i2] - SRE.pointer->fMRI[i1][i2] ) + 
				   ( misplacement_X * misplacement_Y * misplacement_Y / 2.0) * 
				   ( SRE.pointer->fMRI[i1][i2+1] - SRE.pointer->fMRI[i1][i2] ) +
				   ( misplacement_X * misplacement_X * misplacement_Y * misplacement_Y / 4.0) * 
				   ( SRE.pointer->fMRI[i1+1][i2+1] + SRE.pointer->fMRI[i1][i2] - SRE.pointer->fMRI[i1+1][i2] - SRE.pointer->fMRI[i1][i2+1] );
		
        s = (double) SRE.pointer->fMRI[i1][i2] * misplacement_X * misplacement_Y; 
       
		   if ( (double)s == 0.0 && (double)k == 0.0 ) SRE.pointer->ICF[i1][i2] = (double)1.0; // de l'Hopital
           
		   else if ( (double)s == 0.0 && (double)k != 0.0 ) SRE.pointer->ICF[i1][i2] = (double)0.0;
		   
		   else if ( (double)s != 0.0 && (double)k == 0.0 ) SRE.pointer->ICF[i1][i2] = (double)0.0;

		   else  if ( (double)s != 0.0 && (double)k != 0.0 ) SRE.pointer->ICF[i1][i2] = (double)s/k;


		} // y dim
        
	}  // x dim

	std::cout << "Intensity-Curvature Functional Calculated" << endl;
	// calculate ICF(end)


	char outputFile[200];
	FILE * writeimage;

	sprintf(outputFile, "%s%s","ICF_SRE2D-", imageFileName);

	if ((writeimage = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=PAD+1; v < n2+PAD+1; v++ ) { // (a)

		for( int f=PAD+1; f < n1+PAD+1; f++ ) { // b
	
		fwrite(&SRE.pointer->ICF[v][f],sizeof(double),1,writeimage);

		} // b

	} // (a)

	fclose(writeimage);

	} // (save)

	delete SRE.pointer;
	SRE.~SRE2D2013();

}


void runKurtosis(char imageFileName [], int n1, int n2, FILE * savedata)
{

	int maskSize = 2;

	printf("%s%s\n", "Calculation of the Kurtosis of: ", imageFileName);
	fprintf(savedata,"%s\n", "Calculation of the Kurtosis of: ", imageFileName);

	fprintf(savedata,"%s%s\n", "The image file name is: " , imageFileName);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%d\n", "The square mask size is:  ", maskSize);

	Kurtosis2D2020 Kurtosis(n1,n2);

	Kurtosis.allocateData();

	FILE * pf;

	/// read image file (begin)
	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		Kurtosis.~Kurtosis2D2020();
		exit(0);

	} else { // else

	double number;

	for (int i1=0; i1 < n2; i1++) {// x dim
       	
		for (int i2=0; i2 < n1; i2++) { // y dim
		
		// at each iteration of the two for loops
		// the program reads the pixel value from the
		// file containing the image and 
		fread(&number,sizeof(double),1,pf);

		Kurtosis.pointer->Signal[i1][i2] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl; 

		// scale Signal (begin) 
	double max = -Kurtosis.pointer->Signal[0][0];
	double min = Kurtosis.pointer->Signal[0][0];

	for( int v=0; v < n2; v++ ) { // (a)

	    for( int f=0; f < n1; f++ ) { // (b)


			if( Kurtosis.pointer->Signal[v][f] > (double)max ) 
			
					max = (double) Kurtosis.pointer->Signal[v][f];
              
			if( Kurtosis.pointer->Signal[v][f] < (double)min ) 
			
					min = (double) Kurtosis.pointer->Signal[v][f];
	
		}
	}

	// scale (begin)
	for( int v=0; v < n2; v++ ) { // (a)

	    for( int f=0; f < n1; f++ ) { // (b)

				if ( max == min ) Kurtosis.pointer->Signal[v][f] = (double)0.0;

				else { 
					
				Kurtosis.pointer->Signal[v][f] = (double) fabs( (min - (double)Kurtosis.pointer->Signal[v][f] ) / (min - max) );

				Kurtosis.pointer->Signal[v][f] *= SCALE_in;

				// the absolute vale is needed otherwise the minimum is -0 which is nonsense in computers. The value -0 was discovered
				// printing to the screen. Mathematically, the scaling formula gives values in between [0, SCALE] and not [-0, SCALE].
				// Therefore the nonsense is removed taking the absolute value.

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)


	/// Calculate Kurtosis of the image (begins)
	int n7 = ((int)maskSize);  
	int n8 = ((int)maskSize);

	// pad image memory allocation (begins)
	struct pad_data {

		double **pad_Image; // pointer to the matrix entry 

	}*pad_pointer; // pointer to the matrices

	pad_pointer = new pad_data;

	pad_pointer->pad_Image = new double*[n2+n8];

	for( int v=0; v < n2+n8; v++ ) { // (1)
		 
	pad_pointer->pad_Image[v] = new double[n1+n7];

	  } // (1) allocate struct 'pad_data' (end)

	// (2) initialize (begin)
	for( int v=0; v < n2+n8; v++ ) { // (a)

	    for( int f=0; f < n1+n7 ; f++ ) { // (b)
		 
			pad_pointer->pad_Image[v][f] = (double)0.0;

			 } //(b)

		 } //(a)

	for( int v=0; v < n2; v++ ) { // (a)

	    for( int f=0; f < n1; f++ ) { // (b)
		 
			pad_pointer->pad_Image[v+(int)maskSize/2][f+(int)maskSize/2] = (double)Kurtosis.pointer->Signal[v][f];

			 } //(b)

		 } //(a)
	// (2) initialize (end)

	// pad image memory allocation (ends)

	/// calculate average (begins)
	for (int i = 0; i < (int)n2; i++) { // n1

        for (int j = 0; j < (int)n1; j++) { // n2

			double sum = (double)0.0;

	for (int pp = 0; pp < n8; pp++) { // n7

           for (int qq = 0; qq < n7; qq++) { // n8
			
					sum += (double)pad_pointer->pad_Image[i+pp][j+qq];
          
					} // n8 dim
        
				}  // n7 dim 

			Kurtosis.pointer->average[i][j] = ((double)sum / n7*n8); 

		} // n2 dim
        
	}  // n1 dim 
	/// calculate average (ends)

	
	
	/// calculate Kurtosis (begins)
	for (int i = 0; i < (int)n2; i++) { // n1

        for (int j = 0; j < (int)n1; j++) { // n2

			Kurtosis.pointer->m4[i][j] = (double) 0.0;
			
			Kurtosis.pointer->m2[i][j] = (double) 0.0;

	for (int pp = 0; pp < n8; pp++) { // n7

           for (int qq = 0; qq < n7; qq++) { // n8
			
					Kurtosis.pointer->m4[i][j] += (double) pow( ((double) pad_pointer->pad_Image[i+pp][j+qq] - 
						                                         (double) Kurtosis.pointer->average[i][j]), 4.0);

					Kurtosis.pointer->m2[i][j] += (double) pow( ((double) pad_pointer->pad_Image[i+pp][j+qq] - 
						                                         (double) Kurtosis.pointer->average[i][j]), 2.0);
        
					} // n8 dim
        
				}  // n7 dim 

			if ( (double)Kurtosis.pointer->m2[i][j] != (double) 0.0 )
			
				 Kurtosis.pointer->Kurtosis[i][j] = ( ((double)Kurtosis.pointer->m4[i][j]) / 
				                                      ((double)Kurtosis.pointer->m2[i][j] * Kurtosis.pointer->m2[i][j])); 
			
			else Kurtosis.pointer->Kurtosis[i][j] = (double) 0.0;


		} // n2 dim
        
	}  // n1 dim 
	/// calculate Kurtosis (ends


	std::cout << "Kurtosis of the image calculated" << endl;
	// calculate Kurtosis (end)

	FILE * saveimage;
	char PadimageFileName[128];
	sprintf(PadimageFileName, "%s%s", "Pad-", imageFileName);

	// save padded image (begins)
	if ((saveimage = fopen(PadimageFileName,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;
		exit(0);

	} else  { // (save)


	for( int v=0; v < n2+n8; v++ ) { // (a)

		for( int f=0; f < n1+n7; f++ ) 
	
		fwrite(&pad_pointer->pad_Image[v][f],sizeof(double),1,saveimage);

	} // (a)

	fclose(saveimage);

	} // (save)
	// save padded image (ends)

	char outputFile[128];
	sprintf(outputFile, "%s%s", "Kurtosis-", imageFileName);

	if ((saveimage = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;
		exit(0);

	} else  { // (save)


	for( int v=0; v < n2; v++ ) { // (a)

		for( int f=0; f < n1; f++ ) 
	
		fwrite(&Kurtosis.pointer->Kurtosis[v][f],sizeof(double),1,saveimage);

	} // (a)

	fclose(saveimage);

	} // (save)

	
	sprintf(outputFile, "%s%s", "m4-", imageFileName);
	
	if ((saveimage = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < n2; v++ ) { // (a)

		for( int f=0; f < n1; f++ ) 
	
		fwrite(&Kurtosis.pointer->m4[v][f],sizeof(double),1,saveimage);

	} // (a)

	fclose(saveimage);

	} // (save)

	sprintf(outputFile, "%s%s", "m2-", imageFileName);

	if ((saveimage = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < n2; v++ ) { // (a)

		for( int f=0; f < n1; f++ ) 
	
		fwrite(&Kurtosis.pointer->m2[v][f],sizeof(double),1,saveimage);

	} // (a)

	} // (save)

	fclose(saveimage);

	delete Kurtosis.pointer;
	Kurtosis.~Kurtosis2D2020();
	// run Kurtosis (end)
}


void runFOD2D(char imageFileName [], int xd, int yd, FILE * savedata)
{

	int n1 = xd;
	int n2 = yd;

	double XPixelSize = 1.0;
	double YPixelSize = 1.0;

	double x_misplacement_X = 0.5;
	double y_misplacement_Y = 0.5;

	double theta = 0.0;
	
	printf("%s\n", "Calculation of FOD using SRE2D");
	fprintf(savedata,"%s\n", "Calculation of FOD using SRE2D");

	
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n2);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n1);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The XY rotation angle is: ", theta);

    double misplacement_X = ((double)1.0 - ( cos( (double)theta ) + sin( (double)theta ) ) + x_misplacement_X);
    double misplacement_Y = ((double)1.0 - ( -sin( (double)theta ) + cos( (double)theta ) ) + y_misplacement_Y);

      misplacement_X = ((double)misplacement_X/XPixelSize);
      misplacement_Y = ((double)misplacement_Y/YPixelSize);

	  //////////////////***********//////////////////////
	  // Above formula scales the misplacement to the  //
	  // pixel size the same way the following formula //
	  // would do: (min - misplacement)/(min - max)    //  
	  //////////////////***********//////////////////////
	int PAD = 1;

	FOD2D2018 FOD(n1+PAD*2,n2+PAD*2);

	FOD.allocateData();

	/// read image file (begin)
	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		FOD.~FOD2D2018();
		exit(0);

	} else { // else

	double number;

	for( int v=PAD+1; v < n2+PAD+1; v++ ) { // (a)

	    for( int f=PAD+1; f < n1+PAD+1; f++ ) { // (b)

			fread(&number,sizeof(double),1,pf);
		
			FOD.pointer->Signal[v][f] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl;

	// scale Signal (begin) 
	double max = -FOD.pointer->Signal[n2/2][n1/2];
	double min = FOD.pointer->Signal[n2/2][n1/2];

	for( int v=PAD+1; v < n2+PAD+1; v++ ) { // (a)

	    for( int f=PAD+1; f < n1+PAD+1; f++ ) { // (b)

	
			if( FOD.pointer->Signal[v][f] > (double)max ) 
			
					max = (double) FOD.pointer->Signal[v][f];
              
			if( FOD.pointer->Signal[v][f] < (double)min ) 
			
					min = (double) FOD.pointer->Signal[v][f];
	
		}
	}


	for( int v=PAD+1; v < n2+PAD+1; v++ ) { // (a)

	    for( int f=PAD+1; f < n1+PAD+1; f++ ) { // (b)

				if ( max == min ) FOD.pointer->Signal[v][f] = (double)0.0;

				else { 
					
				FOD.pointer->Signal[v][f] = (double) fabs( (min - (double)FOD.pointer->Signal[v][f] ) / (min - max) );

				FOD.pointer->Signal[v][f] *= SCALE_in;

				// the absolute vale is needed otherwise the minimum is -0 which is nonsense in computers. The value -0 was discovered
				// printing to the screen. Mathematically, the scaling formula gives values in between [0, SCALE] and not [-0, SCALE].
				// Therefore the nonsense is removed taking the absolute value.

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)

	double dfxydx = 0, dfxydy = 0, omega = 0;
	
	// calculate FOD of the SRE2D model function (begin)
	for( int v=PAD+2; v < n2+PAD; v++ ) { // (a)

	    for( int f=PAD+2; f < n1+PAD; f++ ) { // (b)


		     dfxydx = ( FOD.pointer->Signal[v-1][f] - FOD.pointer->Signal[v][f] ) + 
			            misplacement_Y * ( FOD.pointer->Signal[v-1][f-1] + FOD.pointer->Signal[v][f] - 
						                   FOD.pointer->Signal[v-1][f] - FOD.pointer->Signal[v][f-1] );

		     dfxydy = ( FOD.pointer->Signal[v][f-1] - FOD.pointer->Signal[v][f] ) + 
			            misplacement_X * ( FOD.pointer->Signal[v-1][f-1] + FOD.pointer->Signal[v][f] - 
						                   FOD.pointer->Signal[v-1][f] - FOD.pointer->Signal[v][f-1] );


			 FOD.pointer->fodSRE2D[v][f] = (double) sqrt ( ((double) dfxydx*dfxydx + dfxydy*dfxydy) );


		} // y dim
        
	}  // x dim

	std::cout << "FOD of the SRE2D model function calculated" << endl;
	// calculate FOD of model function (end)

	FILE * saveimage;
	char outputFile[128];

	sprintf(outputFile, "%s%s","fodSRE2D-", imageFileName);

	if ((saveimage = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=PAD+1; v < n2+PAD+1; v++ ) { // (a)

	    for( int f=PAD+1; f < n1+PAD+1; f++ ) { // (b)
	
		fwrite(&FOD.pointer->fodSRE2D[v][f],sizeof(double),1,saveimage);

		} // (b)

	} // (a)

	fclose(saveimage);

	} // (save)

	delete FOD.pointer;
	FOD.~FOD2D2018();

}

void runFOD2D_nopad(char imageFileName [], int xd, int yd, FILE * savedata)
{

	int n1 = xd;
	int n2 = yd;

	double XPixelSize = 1.0;
	double YPixelSize = 1.0;

	double x_misplacement_X = 0.5;
	double y_misplacement_Y = 0.5;

	double theta = 0.0;
	
	printf("%s\n", "Calculation of FOD using SRE2D");
	fprintf(savedata,"%s\n", "Calculation of FOD using SRE2D");

	
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n2);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n1);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The XY rotation angle is: ", theta);

    double misplacement_X = ((double)1.0 - ( cos( (double)theta ) + sin( (double)theta ) ) + x_misplacement_X);
    double misplacement_Y = ((double)1.0 - ( -sin( (double)theta ) + cos( (double)theta ) ) + y_misplacement_Y);

      misplacement_X = ((double)misplacement_X/XPixelSize);
      misplacement_Y = ((double)misplacement_Y/YPixelSize);

	  //////////////////***********//////////////////////
	  // Above formula scales the misplacement to the  //
	  // pixel size the same way the following formula //
	  // would do: (min - misplacement)/(min - max)    //  
	  //////////////////***********//////////////////////

	FOD2D2018 FOD(n1,n2);

	FOD.allocateData();

	/// read image file (begin)
	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		FOD.~FOD2D2018();
		exit(0);

	} else { // else

	double number;

	for( int v=0; v < n2; v++ ) { // (a)

	    for( int f=0; f < n1; f++ ) { // (b)

			fread(&number,sizeof(double),1,pf);
		
			FOD.pointer->Signal[v][f] = (double)number;
                          
		} // y dim
        
	}  // x dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl;

	// scale Signal (begin) 
	double max = -FOD.pointer->Signal[n2/2][n1/2];
	double min = FOD.pointer->Signal[n2/2][n1/2];

	for( int v=0; v < n2; v++ ) { // (a)

	    for( int f=0; f < n1; f++ ) { // (b)

	
			if( FOD.pointer->Signal[v][f] > (double)max ) 
			
					max = (double) FOD.pointer->Signal[v][f];
              
			if( FOD.pointer->Signal[v][f] < (double)min ) 
			
					min = (double) FOD.pointer->Signal[v][f];
	
		}
	}


	for( int v=0; v < n2; v++ ) { // (a)

	    for( int f=0; f < n1; f++ ) { // (b)

				if ( max == min ) FOD.pointer->Signal[v][f] = (double)0.0;

				else { 
					
				FOD.pointer->Signal[v][f] = (double) fabs( (min - (double)FOD.pointer->Signal[v][f] ) / (min - max) );

				FOD.pointer->Signal[v][f] *= SCALE_in;

				// the absolute vale is needed otherwise the minimum is -0 which is nonsense in computers. The value -0 was discovered
				// printing to the screen. Mathematically, the scaling formula gives values in between [0, SCALE] and not [-0, SCALE].
				// Therefore the nonsense is removed taking the absolute value.

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)

	double dfxydx = 0, dfxydy = 0, omega = 0;
	int PAD = 1, v, f;

	// calculate FOD of the SRE2D model function (begin)
	for( v=PAD; v < n2-1; v++ ) { // (a)

	    for( f=PAD; f < n1-1; f++ ) { // (b)

			
		     dfxydx = ( FOD.pointer->Signal[v+1][f] - FOD.pointer->Signal[v][f] ) + 
			            misplacement_Y * ( FOD.pointer->Signal[v+1][f+1] + FOD.pointer->Signal[v][f] - 
						                   FOD.pointer->Signal[v+1][f] - FOD.pointer->Signal[v][f+1] );

		     dfxydy = ( FOD.pointer->Signal[v][f+1] - FOD.pointer->Signal[v][f] ) + 
			            misplacement_X * ( FOD.pointer->Signal[v+1][f+1] + FOD.pointer->Signal[v][f] - 
						                   FOD.pointer->Signal[v+1][f] - FOD.pointer->Signal[v][f+1] );


			 if ( v == n2-2 || f == n1-2)
			 {

			 dfxydx = ( FOD.pointer->Signal[v-1][f] - FOD.pointer->Signal[v][f] ) + 
			            misplacement_Y * ( FOD.pointer->Signal[v-1][f-1] + FOD.pointer->Signal[v][f] - 
						                   FOD.pointer->Signal[v-1][f] - FOD.pointer->Signal[v][f-1] );

		     dfxydy = ( FOD.pointer->Signal[v][f-1] - FOD.pointer->Signal[v][f] ) + 
			            misplacement_X * ( FOD.pointer->Signal[v-1][f-1] + FOD.pointer->Signal[v][f] - 
						                   FOD.pointer->Signal[v-1][f] - FOD.pointer->Signal[v][f-1] );
			 }

			 FOD.pointer->fodSRE2D[v][f] = (double) sqrt ( ((double) dfxydx*dfxydx + dfxydy*dfxydy) );


		} // y dim
        
	}  // x dim

	std::cout << "FOD of the SRE2D model function calculated" << endl;
	// calculate FOD of model function (end)

	FILE * saveimage;
	char outputFile[128];

	sprintf(outputFile, "%s%s","fodSRE2D-", imageFileName);

	if ((saveimage = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < n2; v++ ) { // (a)

	    for( int f=0; f < n1; f++ ) { // (b)
	
		fwrite(&FOD.pointer->fodSRE2D[v][f],sizeof(double),1,saveimage);

		} // (b)

	} // (a)

	fclose(saveimage);

	} // (save)

	delete FOD.pointer;
	FOD.~FOD2D2018();

}


void OnFourierTransform(char imageFilename[], int rcxres, int rcyres)
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;

	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;

	double * kSpaceR = 0;
	double * kSpaceI = 0;
	double * Signal = 0;

	FILE * logfile;
	
	char logfilename[128]="Fourier-T.log";

  	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

		printf("%s\n %s\n" , "Unable to open log File", "Now Exit");

		exit(0);
	
	} else { // allocate memory 


	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}

	if ((Signal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);
		exit(0);

	}

	} // allocate memory 

	//// read image data and initialize pointers
	double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);
				
				*(kSpaceR+index) = (double) 0.0;

				*(kSpaceI+index) = (double) 0.0;

			}

		}

	FILE * pf;
	char SignalFilename[128];
	double readData;
	
	sprintf(SignalFilename, "%s", imageFilename);

	if ((pf = fopen(SignalFilename,"rb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to read Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // read data


	for (i=0; i<rcyres; i++)
	{ ///read signal data
		for (j=0; j<rcxres; j++)
		{

			index = ((j*rcyres)+i);
          
            fread(&readData,sizeof(double),1,pf);

			*(Signal+index) = (double)readData;

		}
	} ///read signal data

	fprintf(logfile,"%s\n", "Signal Read in DOUBLE (64bits) format");

	fclose (pf);
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

				// the absolute vale is needed otherwise the minimum is -0 which is nonsense in computers. The value -0 was discovered
				// printing to the screen. Mathematically, the scaling formula gives values in between [0, SCALE] and not [-0, SCALE].
				// Therefore the nonsense is removed taking the absolute value.

				}

		} // y dim
        
	}  // x dim
	// scale Signal (end)
	
	double phase, complexR, complexI;
	
	///// Fourier Transform //////
	for (i=0; i<NofYpixels; i++)
	{ ///calculate k-space data

		for (j=0; j<NofXpixels; j++)
		{

	
			dx = ((int) i - NofYpixels/2);
		    dy = ((int) j - NofXpixels/2);

			k2 = ((int)(dy*NofYpixels)+dx); 

			w = ((j*NofYpixels)+i);

			for (int s=0; s<NofYpixels; s++)
			{ ///calculate k-space data 
				for (int p=0; p<NofXpixels; p++)
				{ 
					

		     		ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);
 
				    k3 = ((int)(ds*NofXpixels)+dp); 

					t = ((p*NofYpixels)+s);

					phase = ((double) (2.0 * pi * k2 * k3) / (NofXpixels*NofYpixels));

					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (begin)**/
					complexR = (double) cos( (double)phase ) + (double) sin( (double)phase ); 

					complexI = -(double) sin( (double)phase ) + (double) cos( (double)phase ); 
					//** nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform (end)**/
				
					*(kSpaceR+w) += (double) *(Signal+t) * (double) complexR;

					*(kSpaceI+w) += (double) *(Signal+t) * (double) complexI;

			}

		}///calculate k-space data 

			    
		}
	} ///calculate k-space data

	///// Fourier Transform //////
	double savedata = 0.0;
	char FTfilename[128];

	sprintf(FTfilename, "%s%s", "K-SpaceR-", imageFilename);

    fprintf(logfile, "%s\t%s\n", "Now Saving K-Space Signal (Real) in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");


	 // FIFO memory deallocation method
 	 free(kSpaceR);
 	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(kSpaceR+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "K-Space Signal (Real) Saved");

	fclose (pf);
	} // save data



	sprintf(FTfilename, "%s%s", "K-SpaceI-", imageFilename);

    fprintf(logfile, "%s\t%s\n", "Now Saving K-Space Signal (Imaginary) in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(kSpaceI+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "K-Space Signal (Imaginary) Saved");

	fclose (pf);
	
	} // save data

	sprintf(FTfilename, "%s%s", "K-SpaceM-", imageFilename);

    fprintf_s(logfile, "%s\t%s\n", "Now Saving K-Space Magnitude of the Signal in File: ", FTfilename);

    if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf_s(logfile, "%s\n", "Cannot open file to save K-Space Magnitude of the Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data	

		// save a zero image (begin)
		for (int s=0; s<NofYpixels; s++)
		{ 
			for (int p=0; p<NofXpixels; p++)
			{ 

			savedata = (double)0.0;
          
            fwrite(&savedata,sizeof(double),1,pf);

			}
		} // save a zero image (end)

	fclose(pf);
	
	}
		
	if ((pf = fopen(FTfilename,"wb+"))==NULL)
	{

	 fprintf_s(logfile, "%s\n", "Cannot open file to save K-Space Magnitude of the Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);

	 exit(0);
	
	} else { // save data
		
		// K-Space Magnitude (begin)
		for (int s=0; s<(int)NofYpixels; s++)
		{ 
			for (int p=0; p<(int)NofXpixels; p++)
			{ 
			
		
			index = ((p*NofYpixels)+s);

			savedata = (double) sqrt( (double)*(kSpaceR+index)*(double)*(kSpaceR+index) + 
		   		                      (double)*(kSpaceI+index)*(double)*(kSpaceI+index) );
          
            fwrite(&savedata,sizeof(double),1,pf);
			
		}
	} // K-Space Magnitude (end)

	fprintf_s(logfile,"%s\n", "K-Space Magnitude of the Signal Saved");

	fclose (pf);
	} // save data

	
	printf("%s\n", "FT Processing Completed");
    fprintf_s(logfile,"%s\n", "FT Processing Completed");

	fclose(logfile);

	// FIFO memory deallocation method
	free(kSpaceR);
	free(kSpaceI);
	free(Signal);

}

void OnInverseFourierTransformConvolve(char fodfilename[], char icfFilename[], int rcxres, int rcyres)
{
	
	int NofXpixels = rcxres;
	int NofYpixels = rcyres;
	
	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;	
	double phase;
	double emittingSource = 1.0; 
	double scale = ((double)rcxres*rcyres*emittingSource); 

	FILE * logfile;
	char logfilename[128]="INV-FourierT.log";

	FILE *image;
	char imageFilename[256];

	double * fodkSpaceR = 0;
	double * fodkSpaceI = 0;
	double * icfkSpaceR = 0;
	double * icfkSpaceI = 0;
	double * reconSignal = 0;

	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

	 exit(0);
	
	} else { // allocate memory

		
	printf("%s\n", "Now INV FT Processing...");
    fprintf(logfile,"%s\n", "Now INV FT Processing...");

	if ((fodkSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		exit(0);

	}

	if ((fodkSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		// FIFO memory deallocation method
		free(fodkSpaceR);
		exit(0);

	}

	if ((icfkSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		free(fodkSpaceR);
		free(fodkSpaceI);

		exit(0);

	}

	if ((icfkSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
   
		// FIFO memory deallocation method
		free(fodkSpaceR);
		free(fodkSpaceI);
		free(icfkSpaceR);

		exit(0);

	}

	if ((reconSignal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
	
		fprintf(logfile,"%s\n", "Not enough memory to allocate Imaginary Image data: Exit");
	
		// FIFO memory deallocation method
		free(fodkSpaceR);
		free(fodkSpaceI);
		free(icfkSpaceR);
		free(icfkSpaceI);

		exit(0);

	}
	
	} // allocate memory

	
	//// read image data and initialize pointers
    sprintf(imageFilename, "%s%s", "K-SpaceR-", fodfilename);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	free(fodkSpaceR);
	free(fodkSpaceI);
	free(icfkSpaceR);
	free(icfkSpaceI);
	free(reconSignal);

    exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(fodkSpaceR+index) = (double) number;
		
			}

		}

		fclose(image);

	}// read data and initialize pointers


    char imageFilename2[128];

	sprintf(imageFilename2, "%s%s", "K-SpaceI-", fodfilename);


    if ((image = fopen(imageFilename2,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

    // FIFO memory deallocation method
	free(fodkSpaceR);
	free(fodkSpaceI);
	free(icfkSpaceR);
	free(icfkSpaceI);
	free(reconSignal);
	
	exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(fodkSpaceI+index) = (double) number;

			}

		}

		fclose(image);


		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				*(reconSignal+index) = (double)0.0;
					
			}

		}

	}// read data and initialize pointers


	//// read image data and initialize pointers
    sprintf(imageFilename, "%s%s", "K-SpaceR-", icfFilename);

    if ((image = fopen(imageFilename,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);
	 
	 // FIFO memory deallocation method
	free(fodkSpaceR);
	free(fodkSpaceI);
	free(icfkSpaceR);
	free(icfkSpaceI);
	free(reconSignal);
   
	exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(icfkSpaceR+index) = (double) number;
		
			}

		}

		fclose(image);

	}// read data and initialize pointers


	sprintf(imageFilename2, "%s%s", "K-SpaceI-", icfFilename);
	
    if ((image = fopen(imageFilename2,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename2);

	// FIFO memory deallocation method
	free(fodkSpaceR);
	free(fodkSpaceI);
	free(icfkSpaceR);
	free(icfkSpaceI);
	free(reconSignal);

    exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofYpixels; i++)
		{ 
			for (j=0; j<NofXpixels; j++)
			{

				index = ((j*NofYpixels)+i);

				fread(&number,sizeof(double),1,image);
				
				*(icfkSpaceI+index) = (double) number;

			}

		}

		fclose(image);

	}// read data and initialize pointers


	double real = 0.0, imaginary = 0.0;
    double fodPhase = 0.0, icfPhase = 0.0, phasexy = 0.0;

	///// INV Fourier Transform //////
	for (i=0; i<NofYpixels; i++)
	{ ///process k-space data

		for (j=0; j<NofXpixels; j++)
		{
		
	    	dx = ((int) i - NofYpixels/2);
		    dy = ((int) j - NofXpixels/2);
		
	  	    k2 = ((int)(dx*NofXpixels)+dy);

			w = ((j*NofYpixels)+i);

			real = (double)0.0;
			imaginary = (double)0.0;

			
			for (int s=0; s<NofYpixels; s++)
			{ ///process k-space data

				for (int p=0; p<NofXpixels; p++)
				{ 

					ds = ((int) s - NofYpixels/2);
		            dp = ((int) p - NofXpixels/2);

					k3 = ((int)(dp*NofYpixels)+ds);  
				
					t = ((p*NofYpixels)+s);
					
					phase = ((double) (2.0 * pi * k2 * k3) / (NofXpixels*NofYpixels));

					
					fodPhase = (double)2.0 * pi * atan2( (double) *(fodkSpaceI+t) , (double) *(fodkSpaceR+t)  ) / (NofXpixels*NofYpixels);
				
					icfPhase = (double)2.0 * pi * atan2( (double) *(icfkSpaceI+t) , (double) *(icfkSpaceR+t)  ) / (NofXpixels*NofYpixels);
	
					phasexy = (double)2.0 * pi * atan2( ((double) *(icfkSpaceI+t) - *(fodkSpaceI+t) ), 
					
						                                ((double) *(icfkSpaceR+t) - *(fodkSpaceR+t) ) ) / (NofXpixels*NofYpixels);
								
					fodPhase += (double)phase;

					icfPhase += (double)phase;

					phasexy += (double)phase;
				
					
					real += ((double) *(fodkSpaceR+t) * cos((double)fodPhase)) - ((double) *(fodkSpaceI+t) * sin((double)fodPhase));

					real += ((double) *(icfkSpaceR+t) * cos((double)icfPhase)) - ((double) *(icfkSpaceI+t) * sin((double)icfPhase));
					
					real += ((double)  2.0 * ((double) *(icfkSpaceR+t) -  *(fodkSpaceR+t)) * cos((double)phasexy)) - 
						    ((double)  2.0 * ((double) *(icfkSpaceI+t) -  *(fodkSpaceI+t)) * sin((double)phasexy));

					
					imaginary += ((double) *(fodkSpaceR+t) * sin((double)fodPhase)) + ((double) *(fodkSpaceI+t) * cos((double)fodPhase)); 

					imaginary += ((double) *(icfkSpaceR+t) * sin((double)icfPhase)) + ((double) *(icfkSpaceI+t) * cos((double)icfPhase)); 
					
					imaginary +=  ((double) 2.0 * ((double) *(icfkSpaceR+t) - *(fodkSpaceR+t)) * sin((double)phasexy)) + 
						          ((double) 2.0 * ((double) *(icfkSpaceI+t) - *(fodkSpaceI+t)) * cos((double)phasexy)); 
			
			}

		}///process k-space data 

			*(reconSignal+w) =  (double) sqrt( ((double) real * real)  + ((double) imaginary * imaginary) );

			*(reconSignal+w) /= (double)scale;
		}
	} ///process k-space data


	double savedata = 0.0;
	FILE * pf;
	char reconFilename[128];

	sprintf(reconFilename, "%s", "reconSignal.img");


    fprintf(logfile, "%s\t%s\n", "Now Saving Reconstructed Signal in File: ", reconFilename);

    if ((pf = fopen(reconFilename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(fodkSpaceR);
	 free(fodkSpaceI);
	 free(icfkSpaceR);
	 free(icfkSpaceI);
	 free(reconSignal);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofYpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofXpixels; j++)
		{

			index = ((j*NofYpixels)+i);

			savedata = (double)*(reconSignal+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "Reconstructed Signal Saved");

	fclose (pf);
	} // save data


    printf("%s\n", "Inverse FT Processing Completed");
    fprintf(logfile,"%s\n", "Inverse FT Processing Completed");

	fclose(logfile);
		
	// FIFO memory deallocation method
	free(fodkSpaceR);
	free(fodkSpaceI);
	free(icfkSpaceR);
	free(icfkSpaceI);
	free(reconSignal);

}