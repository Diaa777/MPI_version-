//
//  convolution.c
//
//
//  Created by Josep Lluis Lerida on 11/03/15.
//
// This program calculates the convolution for PPM images.
// The program accepts an PPM image file, a text definition of the kernel matrix and the PPM file for storing the convolution results.
// The program allows to define image partitions for processing large images (>500MB)
// The 2D image is represented by 1D vector for chanel R, G and B. The convolution is applied to each chanel separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

// Structure to store image.
struct imagenppm{
    int altura;
    int ancho;
    char *comentario;
    int maxcolor;
    int P;
    int *R;
    int *G;
    int *B;
};
typedef struct imagenppm* ImagenData;

// Structure to store the kernel.
struct structkernel{
    int kernelX;
    int kernelY;
    float *vkern;
};
typedef struct structkernel* kernelData;

//Functions Definition
ImagenData initimage(char* nombre, FILE **fp, int partitions, int halo);
ImagenData duplicateImageData(ImagenData src, int partitions, int halo);

int readImage(ImagenData Img, FILE **fp, int dim, int halosize, long int *position);
int duplicateImageChunk(ImagenData src, ImagenData dst, int dim);
int initfilestore(ImagenData img, FILE **fp, char* nombre, long *position);
int savingChunk(ImagenData img, FILE **fp, int dim, int offset);
int convolve2D(int* inbuf, int* outbuf, int sizeX, int sizeY, float* kernel, int ksizeX, int ksizeY);
void freeImagestructure(ImagenData *src);

//Open Image file and image struct initialization
ImagenData initimage(char* nombre, FILE **fp,int partitions, int halo){
    char c;
    char comentario[300];
    int i=0,chunk=0;
    ImagenData img=NULL;
    
    /*Opening ppm*/

    if ((*fp=fopen(nombre,"r"))==NULL){
        perror("Error: ");
    }
    else{
        //Memory allocation
        img=(ImagenData) malloc(sizeof(struct imagenppm));

        //Reading the first line: Magical Number "P3"
        fscanf(*fp,"%c%d ",&c,&(img->P));
        
        //Reading the image comment
        while((c=fgetc(*fp))!= '\n'){comentario[i]=c;i++;}
        comentario[i]='\0';
        //Allocating information for the image comment
        img->comentario = calloc(strlen(comentario),sizeof(char));
        strcpy(img->comentario,comentario);
        //Reading image dimensions and color resolution
        fscanf(*fp,"%d %d %d",&img->ancho,&img->altura,&img->maxcolor);
        chunk = img->ancho*img->altura / partitions;
        //We need to read an extra row.
        chunk = chunk + img->ancho * halo;
        if ((img->R=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
        if ((img->G=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
        if ((img->B=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
    }
    return img;
}

//Duplicate the Image struct for the resulting image
ImagenData duplicateImageData(ImagenData src, int partitions, int halo){
    char c;
    char comentario[300];
    unsigned int imageX, imageY;
    int i=0, chunk=0;
    //Struct memory allocation
    ImagenData dst=(ImagenData) malloc(sizeof(struct imagenppm));

    //Copying the magic number
    dst->P=src->P;
    //Copying the string comment
    dst->comentario = calloc(strlen(src->comentario),sizeof(char));
    strcpy(dst->comentario,src->comentario);
    //Copying image dimensions and color resolution
    dst->ancho=src->ancho;
    dst->altura=src->altura;
    dst->maxcolor=src->maxcolor;
    chunk = dst->ancho*dst->altura / partitions;
    //We need to read an extra row.
    chunk = chunk + src->ancho * halo;
    if ((dst->R=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
    if ((dst->G=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
    if ((dst->B=calloc(chunk,sizeof(int))) == NULL) {return NULL;}
    return dst;
}

//Read the corresponding chunk from the source Image
int readImage(ImagenData img, FILE **fp, int dim, int halosize, long *position){
    int i=0, k=0,haloposition=0;
    if (fseek(*fp,*position,SEEK_SET))
        perror("Error: ");
    haloposition = dim-(img->ancho*halosize*2);
    for(i=0;i<dim;i++) {
        // When start reading the halo store the position in the image file
        if (halosize != 0 && i == haloposition) *position=ftell(*fp);
        fscanf(*fp,"%d %d %d ",&img->R[i],&img->G[i],&img->B[i]);
        k++;
    }
//    printf ("Readed = %d pixels, posicio=%lu\n",k,*position);
    return 0;
}

//Duplication of the  just readed source chunk to the destiny image struct chunk
int duplicateImageChunk(ImagenData src, ImagenData dst, int dim){
    int i=0;
    
    for(i=0;i<dim;i++){
        dst->R[i] = src->R[i];
        dst->G[i] = src->G[i];
        dst->B[i] = src->B[i];
    }
//    printf ("Duplicated = %d pixels\n",i);
    return 0;
}

// Open kernel file and reading kernel matrix. The kernel matrix 2D is stored in 1D format.
kernelData leerKernel(char* nombre){
    FILE *fp;
    int i=0;
    kernelData kern=NULL;
    
    /*Opening the kernel file*/
    fp=fopen(nombre,"r");
    if(!fp){
        perror("Error: ");
    }
    else{
        //Memory allocation
        kern=(kernelData) malloc(sizeof(struct structkernel));
        
        //Reading kernel matrix dimensions
        fscanf(fp,"%d,%d,", &kern->kernelX, &kern->kernelY);
        kern->vkern = (float *)malloc(kern->kernelX*kern->kernelY*sizeof(float));
        
        // Reading kernel matrix values
        for (i=0;i<(kern->kernelX*kern->kernelY)-1;i++){
            fscanf(fp,"%f,",&kern->vkern[i]);
        }
        fscanf(fp,"%f",&kern->vkern[i]);
        fclose(fp);
    }
    return kern;
}

// Open the image file with the convolution results
int initfilestore(ImagenData img, FILE **fp, char* nombre, long *position){
    /*Se crea el fichero con la imagen resultante*/
    if ( (*fp=fopen(nombre,"w")) == NULL ){
        perror("Error: ");
        return -1;
    }
    /*Writing Image Header*/
    fprintf(*fp,"P%d\n%s\n%d %d\n%d\n",img->P,img->comentario,img->ancho,img->altura,img->maxcolor);
    *position = ftell(*fp);
    return 0;
}

// Writing the image partition to the resulting file. dim is the exact size to write. offset is the displacement for avoid halos.
int savingChunk(ImagenData img, FILE **fp, int dim, int offset){
    int i,k=0;
    //Writing image partition
    for(i=offset;i<dim+offset;i++){
        fprintf(*fp,"%d %d %d ",img->R[i],img->G[i],img->B[i]);
//        if ((i+1)%6==0) fprintf(*fp,"\n");
        k++;
    }
//    printf ("Writed = %d pixels, dim=%d, offset=%d\n",k,dim, offset);
    return 0;
}

// This function free the space allocated for the image structure.
void freeImagestructure(ImagenData *src){
    
    free((*src)->comentario);
    free((*src)->R);
    free((*src)->G);
    free((*src)->B);
    
    free(*src);
}

///////////////////////////////////////////////////////////////////////////////
// 2D convolution
// 2D data are usually stored in computer memory as contiguous 1D array.
// So, we are using 1D array for 2D data.
// 2D convolution assumes the kernel is center originated, which means, if
// kernel size 3 then, k[-1], k[0], k[1]. The middle of index is always 0.
// The following programming logics are somewhat complicated because of using
// pointer indexing in order to minimize the number of multiplications.
//
//
// signed integer (32bit) version:
///////////////////////////////////////////////////////////////////////////////
int convolve2D(int* in, int* out, int dataSizeX, int dataSizeY,
               float* kernel, int kernelSizeX, int kernelSizeY)
{
    int i, j, m, n;
    int *inPtr, *inPtr2, *outPtr;
    float *kPtr;
    int kCenterX, kCenterY;
    int rowMin, rowMax;                             // to check boundary of input array
    int colMin, colMax;                             //
    float sum;                                      // temp accumulation buffer
    
    // check validity of params
    if(!in || !out || !kernel) return -1;
    if(dataSizeX <= 0 || kernelSizeX <= 0) return -1;
    
    // find center position of kernel (half of kernel size)
    kCenterX = (int)kernelSizeX / 2;
    kCenterY = (int)kernelSizeY / 2;
    
    // init working  pointers
    inPtr = inPtr2 = &in[dataSizeX * kCenterY + kCenterX];  // note that  it is shifted (kCenterX, kCenterY),
    outPtr = out;
    kPtr = kernel;
    
    // start convolution
	
	i= 0;
//	stopPoint=(dataSizeY*(pid+1)-1);

	
    for(i ; i < dataSizeY; ++i)                   // number of rows
    {
        // compute the range of convolution, the current row of kernel should be between these
        rowMax = i + kCenterY;
        rowMin = i - dataSizeY + kCenterY;
        
        for(j = 0; j < dataSizeX; ++j)              // number of columns
        {
            // compute the range of convolution, the current column of kernel should be between these
            colMax = j + kCenterX;
            colMin = j - dataSizeX + kCenterX;
            
            sum = 0;                                // set to 0 before accumulate
            
            // flip the kernel and traverse all the kernel values
            // multiply each kernel value with underlying input data
            for(m = 0; m < kernelSizeY; ++m)        // kernel rows
            {
                // check if the index is out of bound of input array
                if(m <= rowMax && m > rowMin)
                {
                    for(n = 0; n < kernelSizeX; ++n)
                    {
                        // check the boundary of array
                        if(n <= colMax && n > colMin)
                            sum += *(inPtr - n) * *kPtr;
                        
                        ++kPtr;                     // next kernel
                    }
                }
                else
                    kPtr += kernelSizeX;            // out of bound, move to next row of kernel
                
                inPtr -= dataSizeX;                 // move input data 1 raw up
            }
            
            // convert integer number
            if(sum >= 0) *outPtr = (int)(sum + 0.5f);
//            else *outPtr = (int)(sum - 0.5f)*(-1);
            // For using with image editors like GIMP or others...
            else *outPtr = (int)(sum - 0.5f);
            // For using with a text editor that read ppm images like libreoffice or others...
//            else *outPtr = 0;
            
            kPtr = kernel;                          // reset kernel to (0,0)
            inPtr = ++inPtr2;                       // next input
            ++outPtr;                               // next output
        }
    }
    
    return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION
//////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    int i=0,j=0,k=0;
//    int headstored=0, imagestored=0, stored;
    
/*    if(argc != 5)
    {
        printf("Usage: %s <image-file> <kernel-file> <result-file> <partitions>\n", argv[0]);
        
        printf("\n\nError, Missing parameters:\n");
        printf("format: ./serialconvolution image_file kernel_file result_file\n");
        printf("- image_file : source image path (*.ppm)\n");
        printf("- kernel_file: kernel path (text file with 1D kernel matrix)\n");
        printf("- result_file: result image path (*.ppm)\n");
        printf("- partitions : Image partitions\n\n");
        return -1;
    }
    

*/
	
	

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // READING IMAGE HEADERS, KERNEL Matrix, DUPLICATE IMAGE DATA, OPEN RESULTING IMAGE FILE
    //////////////////////////////////////////////////////////////////////////////////////////////////
    int imagesize, partitions, partsize, chunksize, halo, halosize;
    long position=0;
    double start, tstart=0, tend=0, tread=0, tcopy=0, tconv=0, tstore=0, treadk=0;
	
	int rank,nproc;
	int tagm=11;
	int tags=12;
		
	MPI_Status status;
	
	MPI_Request send_request, recv_request;	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);



    struct timeval tim;
    FILE *fpsrc=NULL,*fpdst=NULL;
    ImagenData source=NULL, output=NULL, output1=NULL;

    // Store number of partitions
    partitions = atoi(argv[4]);
    ////////////////////////////////////////
    //Reading kernel matrix
    gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
    tstart = start;
    kernelData kern=NULL;
    if ( (kern = leerKernel(argv[2]))==NULL) {
        //        free(source);
        //        free(output);
        return -1;
    }
    //The matrix kernel define the halo size to use with the image. The halo is zero when the image is not partitioned.
    if (partitions==1) halo=0;
    else halo = (kern->kernelY/2)*2;
    gettimeofday(&tim, NULL);
    treadk = treadk + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);



    // Master Code : 


	


	if(argc != 5)
    {
        printf("Usage: %s <image-file> <kernel-file> <result-file> <partitions>\n", argv[0]);

        printf("\n\nError, Missing parameters:\n");
        printf("format: ./serialconvolution image_file kernel_file result_file\n");
        printf("- image_file : source image path (*.ppm)\n");
        printf("- kernel_file: kernel path (text file with 1D kernel matrix)\n");
        printf("- result_file: result image path (*.ppm)\n");
        printf("- partitions : Image partitions\n\n");
        return -1;
    }




    ////////////////////////////////////////
    //Reading Image Header. Image properties: Magical number, comment, size and color resolution.
    gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
    //Memory allocation based on number of partitions and halo size.
    if ( (source = initimage(argv[1], &fpsrc, partitions, halo)) == NULL) {
        return -1;
    }
    gettimeofday(&tim, NULL);
    tread = tread + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
    
    //Duplicate the image struct.
    gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
    if ( (output = duplicateImageData(source, partitions, halo)) == NULL) {
        return -1;
    }
    gettimeofday(&tim, NULL);
    tcopy = tcopy + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
    
    ////////////////////////////////////////
    //Initialize Image Storing file. Open the file and store the image header.
    gettimeofday(&tim, NULL);
    start = tim.tv_sec+(tim.tv_usec/1000000.0);
    if (initfilestore(output, &fpdst, argv[3], &position)!=0) {
        perror("Error: ");
        //        free(source);
        //        free(output);
        return -1;
    }
    gettimeofday(&tim, NULL);
    tstore = tstore + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
	

//master



	if (rank==0){ 
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // CHUNK READING
    //////////////////////////////////////////////////////////////////////////////////////////////////
    int c=0, offset=0;
	long chunksize1=0;
	int newOffset=0;
	int *PTRR=NULL, *PTRG=NULL, *PTRB=NULL;
	newOffset=(source->altura/nproc)*partitions;
	
	imagesize=source->altura*source->ancho;
	chunksize1=newOffset*source->ancho;
	partsize  = (source->altura*source->ancho)/partitions;	

	PTRR=source->R;
	PTRG=source->G;
	PTRB=source->B;
   // imagesize = source->altura*source->ancho;
   //partsize  = (source->altura*source->ancho)/partitions;
//    printf("%s ocupa %dx%d=%d pixels. Partitions=%d, halo=%d, partsize=%d pixels\n", argv[1], source->altura, source->ancho, imagesize, partitions, halo, partsize);
    while (c < partitions) {
        ////////////////////////////////////////////////////////////////////////////////
        //Reading Next chunk.
        gettimeofday(&tim, NULL);
        start = tim.tv_sec+(tim.tv_usec/1000000.0);
        if (c==0) {
            halosize  = halo/2;
            chunksize = partsize + (source->ancho*halosize);
            offset   = 0;
        }
        else if(c<partitions-1) {
            halosize  = halo;
            chunksize = partsize + (source->ancho*halosize);
            offset    = (source->ancho*halo/2);
        }
        else {
            halosize  = halo/2;
            chunksize = partsize + (source->ancho*halosize);
            offset    = (source->ancho*halo/2);
        }
        //DEBUG 
//        printf("\nRound = %d, position = %ld, partsize= %d, chunksize=%d pixels\n", c, position, partsize, chunksize);
        
        if (readImage(source, &fpsrc, chunksize, halo/2, &position)) {
            return -1;
        }
        gettimeofday(&tim, NULL);
        tread = tread + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
        
        //Duplicate the image chunk
        gettimeofday(&tim, NULL);
        start = tim.tv_sec+(tim.tv_usec/1000000.0);
        if ( duplicateImageChunk(source, output, chunksize) ) {
            return -1;
        }
        //DEBUG
//        for (i=0;i<chunksize;i++)
//            if (source->R[i]!=output->R[i] || source->G[i]!=output->G[i] || source->B[i]!=output->B[i]) printf("At position i=%d %d!=%d,%d!=%d,%d!=%d\n",i,source->R[i],output->R[i], source->G[i],output->G[i],source->B[i],output->B[i]);

        gettimeofday(&tim, NULL);
        tcopy = tcopy + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
	int value=0;
	int i=1;

//	printf("i am the master and my rank equal %d and my chunk start from %d\n", rank, source->R);


			
	//int endi=nproc+1;
	source->R=source->R+chunksize1;
        source->G=source->G+chunksize1;
        source->B=source->B+chunksize1;

	for ( i; i<nproc; i++){
	MPI_Send(source->R, chunksize1, MPI_INT, i, tagm, MPI_COMM_WORLD);
	MPI_Send(source->G, chunksize1, MPI_INT, i, tagm, MPI_COMM_WORLD);
	MPI_Send(source->B, chunksize1, MPI_INT, i, tagm, MPI_COMM_WORLD);
	source->R=source->R+chunksize1;
	source->G=source->G+chunksize1;
	source->B=source->B+chunksize1;
	//value=value+newOffset;
	} // end for 
        

	   c++;
    }// end while






	printf("i am the master and my rank equal %d\n my chunk equal %d\n", rank, chunksize1);

        convolve2D(PTRR, output->R, source->ancho, (newOffset/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(PTRG, output->G, source->ancho, (newOffset/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(PTRB, output->B, source->ancho, (newOffset/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);

        gettimeofday(&tim, NULL);
        tconv = tconv + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);

         printf("%.6lf seconds elapsed for make the convolution.\n", tconv);


	
	int j=1;
	int endj=nproc+1;

	 output->R=output->R+chunksize1;
         output->G=output->G+chunksize1;
         output->B=output->B+chunksize1;

	for (j; j<nproc; j++){
	MPI_Recv(output->R, chunksize1, MPI_INT, j, tags,  MPI_COMM_WORLD, &status);
	MPI_Recv(output->G, chunksize1, MPI_INT, j, tags, MPI_COMM_WORLD, &status);
	MPI_Recv(output->B, chunksize1, MPI_INT, j, tags, MPI_COMM_WORLD, &status);

	//Storing resulting image partition.
        gettimeofday(&tim, NULL);
        start = tim.tv_sec+(tim.tv_usec/1000000.0);
        if (savingChunk(output, &fpdst, chunksize1, offset)) {
            perror("Error: ");
            //        free(source);
            //        free(output);
            return -1;
        }
        gettimeofday(&tim, NULL);
        tstore = tstore + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
	//freeImagestructure(&output);
	
	//printf("cunk received from node %d and start job from %d\n " , j, output->R);
	 output->R=output->R+chunksize1;
	 output->G=output->G+chunksize1;
	 output->B=output->B+chunksize1;

	 printf("%.6lf seconds elapsed for writing the resulting image.\n", tstore);	

	}// end for


	
	//-------------------------------------------------------
    //printf("%.6lf seconds elapsed for writing the resulting image.\n", tstore);

	//freeImagestructure(&source);
        //freeImagestructure(&output);
	fclose(fpsrc);
        fclose(fpdst);

}// end Master code .


	// Slaves code

	if (rank!=0){ 
	
		
	int offset=0;
        int newOffset;
	long chunksize1;
        newOffset=(source->altura/nproc)*partitions;
        chunksize1=newOffset*source->ancho;

	//memory allocation to save the chunk
	
	output1=(ImagenData)malloc(sizeof(struct imagenppm));

	 output1->R=calloc(chunksize1,sizeof(int));
	 output1->G=calloc(chunksize1,sizeof(int));
	 output1->B=calloc(chunksize1,sizeof(int));




	//Duplicate the image chunk
	 

/*	 gettimeofday(&tim, NULL);
	  start = tim.tv_sec+(tim.tv_usec/1000000.0);
	  if ( duplicateImageChunk(source, output1, chunksize1) ) {
	 return -1;
         }
*/


	MPI_Recv(output1->R, chunksize1, MPI_INT, 0, tagm, MPI_COMM_WORLD, &status);
	MPI_Recv(output1->G, chunksize1, MPI_INT, 0, tagm, MPI_COMM_WORLD, &status);
	MPI_Recv(output1->B, chunksize1, MPI_INT, 0, tagm, MPI_COMM_WORLD, &status);

        ///////////////////////////////////////////////////////////////

        // CHUNK CONVOLUTION
        //////////////////////////////////////////////////////////////////////////////////////////////////
        gettimeofday(&tim, NULL);
        start = tim.tv_sec+(tim.tv_usec/1000000.0);
        
        //convolve2D(source->R, output->R, source->ancho, (source->altura/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
        //convolve2D(source->G, output->G, source->ancho, (source->altura/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
        //convolve2D(source->B, output->B, source->ancho, (source->altura/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);




	convolve2D(output1->R, output->R, source->ancho, (newOffset/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(output1->G, output->G, source->ancho, (newOffset/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
        convolve2D(output1->B, output->B, source->ancho, (newOffset/partitions)+halosize, kern->vkern, kern->kernelX, kern->kernelY);
        
        gettimeofday(&tim, NULL);
        tconv = tconv + (tim.tv_sec+(tim.tv_usec/1000000.0) - start);
     
        //////////////////////////////////////////////////////////////////////////////////////////////////
        
       
	MPI_Send(output->R, chunksize1, MPI_INT, 0, tags, MPI_COMM_WORLD);
	MPI_Send(output->G, chunksize1, MPI_INT, 0, tags, MPI_COMM_WORLD);
	MPI_Send(output->B, chunksize1, MPI_INT, 0, tags, MPI_COMM_WORLD);
    

    fclose(fpsrc);
    fclose(fpdst);
    
//    freeImagestructure(&source);
//    freeImagestructure(&output);
    
    gettimeofday(&tim, NULL);
    tend = tim.tv_sec+(tim.tv_usec/1000000.0);
    printf("i am worker number %d and my job start from pixel %d\n", rank, *output->R);
    printf("my chunk of pixels equal %d\n", chunksize1);
   // printf("Imatge: %s\n", argv[1]);
    printf("ISizeX : %d\n", source->ancho);
    //printf("ISizeY : %d\n", source->altura);
    printf("IsizeY : %d\n", newOffset);
    printf("kSizeX : %d\n", kern->kernelX);
    printf("kSizeY : %d\n", kern->kernelY);
    //printf("%.6lf seconds elapsed for Reading image file.\n", tread);
    printf("%.6lf seconds elapsed for copying image structure.\n", tcopy);
    printf("%.6lf seconds elapsed for Reading kernel matrix.\n", treadk);
    printf("%.6lf seconds elapsed for make the convolution.\n", tconv);

    printf("%.6lf seconds elapsed\n", tend-tstart);
    
    //freeImagestructure(&source);
     freeImagestructure(&output1);
     freeImagestructure(&output);
    


	}// end slaves code 
    
	//freeImagestructure(&source);
        //freeImagestructure(&output);

	//fclose(fpsrc);
	//fclose(fpdst);

	MPI_Finalize();

	return 0; 

}

