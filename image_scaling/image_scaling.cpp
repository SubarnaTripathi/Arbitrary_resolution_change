// image_scaling.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdio.h>

#define PHASE_BITS 4
#define NB_PHASES  (1 << PHASE_BITS)
#define NB_TAPS    4
#define FCENTER    1  /* index of the center of the filter */
#define FILTER_BITS   8
#define LINE_BUF_HEIGHT (NB_TAPS * 4)
#define POS_FRAC_BITS 16
#define POS_FRAC      (1 << POS_FRAC_BITS)

#define FILTER_SHIFT 15
#define FELEM int16_t
#define FELEM_MAX   32767 
#define FELEM_MIN   -32768 

#define M_PI   3.14159265358979323846 
#define FFMAX(a,b) ((a) > (b) ? (a) : (b))
#define FFMIN(a,b) ((a) > (b) ? (b) : (a))
#define MY_ABS(A) ((A)<(0) ? (-(A)):(A))

typedef unsigned char uint8_t;
typedef short int16_t;

/* input */
#define XSIZE 1920
#define YSIZE 1280
uint8_t img[XSIZE * YSIZE];
uint8_t img1[XSIZE * YSIZE];
uint8_t img2[XSIZE * YSIZE];

struct ImgReSampleContext
{
	int iwidth, iheight, owidth, oheight;
	int h_incr, v_incr;
	int16_t h_filters[NB_PHASES][NB_TAPS]; /* horizontal filters */
	int16_t v_filters[NB_PHASES][NB_TAPS]; /* vertical filters */
	uint8_t *line_buf;
};

typedef struct MyPicture {
     uint8_t *data[4];
     int linesize[4];       
 } MyPicture;

static inline int get_phase(int pos)
{
     return ((pos) >> (POS_FRAC_BITS - PHASE_BITS)) & ((1 << PHASE_BITS) - 1);
}

static void h_resample_fast(uint8_t *dst, int dst_width, const uint8_t *src,
                             int src_width, int src_start, int src_incr,
                             int16_t *filters)
 {
     int src_pos, phase, sum, i;
     const uint8_t *s;
     int16_t *filter;
 
     src_pos = src_start;
     for(i=0;i<dst_width;i++) {
         s = src + (src_pos >> POS_FRAC_BITS);
         phase = get_phase(src_pos);
         filter = filters + phase * NB_TAPS;
 #if NB_TAPS == 4
         sum = s[0] * filter[0] +
             s[1] * filter[1] +
             s[2] * filter[2] +
             s[3] * filter[3];
 #else
         {
             int j;
             sum = 0;
             for(j=0;j<NB_TAPS;j++)
                 sum += s[j] * filter[j];
         }
 #endif
         sum = sum >> FILTER_BITS;
         if (sum < 0)
             sum = 0;
         else if (sum > 255)
             sum = 255;
         dst[0] = sum;
         src_pos += src_incr;
         dst++;
     }
 }
static void v_resample(uint8_t *dst, int dst_width, const uint8_t *src,
                        int wrap, int16_t *filter)
 {
     int sum, i;
     const uint8_t *s;

     s = src;
     for(i=0;i<dst_width;i++) {
 #if NB_TAPS == 4
         sum = s[0 * wrap] * filter[0] +
             s[1 * wrap] * filter[1] +
             s[2 * wrap] * filter[2] +
             s[3 * wrap] * filter[3];
 #else
        {
            int j;
             uint8_t *s1 = s;
 
            sum = 0;
             for(j=0;j<NB_TAPS;j++) {
                 sum += s1[0] * filter[j];
                 s1 += wrap;
             }
         }
 #endif
         sum = sum >> FILTER_BITS;
         if (sum < 0)
             sum = 0;
         else if (sum > 255)
             sum = 255;
         dst[0] = sum;
         dst++;
         s++;
     }
 }
/* slow version to handle limit cases. Does not need optimisation */
static void h_resample_slow(uint8_t *dst, int dst_width,
                            const uint8_t *src, int src_width,
                            int src_start, int src_incr, int16_t *filters)
 {
     int src_pos, phase, sum, j, v, i;
     const uint8_t *s, *src_end;
     int16_t *filter;
 
     src_end = src + src_width;
     src_pos = src_start;
     for(i=0;i<dst_width;i++) {
         s = src + (src_pos >> POS_FRAC_BITS);
         phase = get_phase(src_pos);
         filter = filters + phase * NB_TAPS;
         sum = 0;
         for(j=0;j<NB_TAPS;j++) {
             if (s < src)
                 v = src[0];
             else if (s >= src_end)
                 v = src_end[-1];
             else
                 v = s[0];
             sum += v * filter[j];
             s++;
         }
         sum = sum >> FILTER_BITS;
        if (sum < 0)
            sum = 0;
        else if (sum > 255)
            sum = 255;
         dst[0] = sum;
         src_pos += src_incr;
         dst++;
     }
 }
 
static void h_resample(uint8_t *dst, int dst_width, const uint8_t *src,
                        int src_width, int src_start, int src_incr,
                        int16_t *filters)
{
    int n, src_end;

     if (src_start < 0) 
	 {
         n = (0 - src_start + src_incr - 1) / src_incr;
         h_resample_slow(dst, n, src, src_width, src_start, src_incr, filters);
         dst += n;
         dst_width -= n;
         src_start += n * src_incr;
     }
     src_end = src_start + dst_width * src_incr;
     if (src_end > ((src_width - NB_TAPS) << POS_FRAC_BITS)) 
	      n = (((src_width - NB_TAPS + 1) << POS_FRAC_BITS) - 1 - src_start) / src_incr;
     else 
	      n = dst_width;
     h_resample_fast(dst, n, src, src_width, src_start, src_incr, filters);
     if (n < dst_width) 
	 {
         dst += n;
         dst_width -= n;
         src_start += n * src_incr;
         h_resample_slow(dst, dst_width,src, src_width, src_start, src_incr, filters);
	 }
 }
static void component_resample(ImgReSampleContext *s, 
                                uint8_t *output, int owrap, int owidth, int oheight,
                                uint8_t *input, int iwrap, int iwidth, int iheight)
 {
     int src_y, src_y1, last_src_y, ring_y, phase_y, y1, y;
     uint8_t *new_line, *src_line;
 
     last_src_y = - FCENTER - 1;
     /* position of the bottom of the filter in the source image */
     src_y = (last_src_y + NB_TAPS) * POS_FRAC; 
     ring_y = NB_TAPS; /* position in ring buffer */
     for(y=0;y<oheight;y++) 
	 {
         /* apply horizontal filter on new lines from input if needed */
         src_y1 = src_y >> POS_FRAC_BITS;
         while (last_src_y < src_y1) {
             if (++ring_y >= LINE_BUF_HEIGHT + NB_TAPS)
                 ring_y = NB_TAPS;
             last_src_y++;
             /* handle limit conditions : replicate line (slightly
               inefficient because we filter multiple times) */
             y1 = last_src_y;
             if (y1 < 0) 
			 {
                 y1 = 0;
             } 
			 else if (y1 >= iheight) 
			 {
                 y1 = iheight - 1;
             }
             src_line = input + y1 * iwrap;
             new_line = s->line_buf + ring_y * owidth;
             /* apply filter and handle limit cases correctly */
             h_resample(new_line, owidth,src_line, iwidth, - FCENTER * POS_FRAC, s->h_incr,&s->h_filters[0][0]);
             /* handle ring buffer wraping */
             if (ring_y >= LINE_BUF_HEIGHT) 
			 {
                memcpy(s->line_buf + (ring_y - LINE_BUF_HEIGHT) * owidth, new_line, owidth);
             }
         }
         /* apply vertical filter */
         phase_y = get_phase(src_y);
         v_resample(output, owidth,s->line_buf + (ring_y - NB_TAPS + 1) * owidth, owidth,&s->v_filters[phase_y][0]);
         src_y += s->v_incr;
         output += owrap;
     }
 }
//double bessel(double x)
// {
//	 /*0th order modified bessel function of the first kind.*/
//    double v=1;
//    double t=1;
//    int i;
//     
//    for(i=1; i<50; i++)
//	{
//         t *= i;
//         v += pow(x*x/4, i)/(t*t);
//    }
//    return v;
// }
static inline int clip(int a, int amin, int amax)
 {
     if (a < amin)
         return amin;
     else if (a > amax)
         return amax;
     else
        return a;
 }

static inline long int lrintf(double x)
 {
	 /* XXX: incorrect, but make it compile */
     return (int)(x + (x < 0 ? -0.5 : 0.5));
 }


void build_filter(FELEM *filter, double factor, int tap_count, int phase_count, int scale, int type)
 {
     int ph, i, v;
     double x, y, w, tab[10];
     const int center= (tap_count-1)/2;
 
    /* if upsampling, only need to interpolate, no filter */
     if (factor > 1.0)
         factor = 1.0;
 
     for(ph=0;ph<phase_count;ph++) 
	 {
         double norm = 0;
         double e= 0;
         for(i=0;i<tap_count;i++) 
		 {
            switch(type)
			 {
			default:
            case 0:
				{
					const float d= -0.5; /*first order derivative = -0.5*/
					x = MY_ABS(((double)(i - center) - (double)ph / phase_count) * factor);
					if(x<1.0) 
						y= 1 - 3*x*x + 2*x*x*x + d*(-x*x + x*x*x);
					else      
						y=  d*(-4 + 8*x - 5*x*x + x*x*x);
					break;
				}
             /*case 1:
                 w = 2.0*x / (factor*tap_count) + M_PI;
                 y *= 0.3635819 - 0.4891775 * cos(w) + 0.1365995 * cos(2*w) - 0.0106411 * cos(3*w);
                 break;
             case 2:
                 w = 2.0*x / (factor*tap_count*M_PI);
                 y *= bessel(16*sqrt(FFMAX(1-w*w, 0)));
                break;*/
			}
            tab[i] = y;
            norm += y;
         }
 
        /* normalize so that an uniform color remains the same */
         for(i=0;i<tap_count;i++) 
		 {
             v = (int)(clip(lrintf(tab[i] * scale / norm + e), FELEM_MIN, FELEM_MAX));
             filter[ph * tap_count + i] = v;
             e += tab[i] * scale / norm - v;
         }
     }
 }

void *av_malloc(unsigned int size)
 {
     void *ptr;
      /* lets disallow possible ambiguous cases */
     if(size > INT_MAX)
         return NULL;
     ptr = malloc(size);
     return ptr;
 }

void *av_mallocz(unsigned int size)
{
     void *ptr;
     
     ptr = av_malloc(size);
     if (!ptr)
         return NULL;
     memset(ptr, 0, size);
     return ptr;
}
void av_free(void *ptr)
{
     /* XXX: this test should not be needed on most libcs */
     if (ptr)
         free(ptr);
}

 
ImgReSampleContext *img_resample_full_init(int owidth, int oheight,int iwidth, int iheight)
 {
     ImgReSampleContext *s;
 
     s = (ImgReSampleContext*)av_mallocz(sizeof(ImgReSampleContext));
     if (!s)
         return NULL;
     if((unsigned)owidth >= UINT_MAX / (LINE_BUF_HEIGHT + NB_TAPS))
         return NULL;
     s->line_buf = (uint8_t*)av_mallocz(owidth * (LINE_BUF_HEIGHT + NB_TAPS));
     if (!s->line_buf)
	 {
         av_free(s);
		 return NULL;
	 }
     
     s->owidth = owidth;
     s->oheight = oheight;
     s->iwidth = iwidth;
     s->iheight = iheight;
          
     s->h_incr = (iwidth * POS_FRAC) / owidth;
     s->v_incr = (iheight* POS_FRAC) / oheight; 
 
	 /** last argument should be 0, 1 or 2 ***/
     build_filter(&s->h_filters[0][0], (float) owidth  / (float) (iwidth), NB_TAPS, NB_PHASES, 1<<FILTER_BITS, 0);
     build_filter(&s->v_filters[0][0], (float) oheight / (float) (iheight ), NB_TAPS, NB_PHASES, 1<<FILTER_BITS, 0);

     return s;
 }

ImgReSampleContext *img_resample_init(int owidth, int oheight, int iwidth, int iheight)
 {
    return img_resample_full_init(owidth, oheight, iwidth, iheight);
 }

void img_resample_close(ImgReSampleContext *s)
 {
     av_free(s->line_buf);
     av_free(s);
 }

void img_resample(ImgReSampleContext *s, MyPicture *output, const MyPicture *input)
 {
     int i, shift;
     uint8_t* optr;

     for (i=0;i<3;i++) {
         shift = (i == 0) ? 0 : 1;
 
         optr = output->data[i];
 
         component_resample(s, optr, output->linesize[i], 
                 s->owidth >> shift, s->oheight >> shift,
                 input->data[i],
                 input->linesize[i], s->iwidth >> shift,
                            s->iheight >> shift);
     }
 }

void polyphase(int xsize_s, int ysize_s, int xsize, int ysize, FILE *ft)
 {
    ImgReSampleContext *s;  
            
	/***** for y ****/
	s = img_resample_full_init(xsize, ysize, xsize_s, ysize_s);
	//dump_filter(&s->h_filters[0][0]);
	component_resample(s, img2, s->owidth, s->owidth, s->oheight,img, s->iwidth, s->iwidth, s->iheight );
    img_resample_close(s);
    fwrite(img2,1,xsize*ysize,ft);
	memcpy(img1, img2,xsize*ysize);
	
	/**** for u ****/
	s = img_resample_full_init(xsize/2, ysize/2, xsize_s/2, ysize_s/2);
    component_resample(s, img2, xsize/2, xsize/2, ysize/2,img+(xsize_s*ysize_s) , xsize_s/2, xsize_s/2, ysize_s/2 );
    img_resample_close(s);
    fwrite(img2,1,(xsize/2)*(ysize/2),ft);
	memcpy(img1 + (xsize*ysize), img2,(xsize*ysize)/4);
		

	/**** for v ****/
	s = img_resample_full_init(xsize/2, ysize/2, xsize_s/2, ysize_s/2);
	component_resample(s, img2, xsize/2, xsize/2, ysize/2,img+(int)(xsize_s*ysize_s*1.25), xsize_s/2, xsize_s/2, ysize_s/2 );
	img_resample_close(s);
	fwrite(img2,1,(xsize/2)*(ysize/2),ft);
	memcpy(img1 + (int)(xsize*ysize*1.25), img2,(xsize*ysize)/4);

	printf("\n\n\n one frame(polyphase) write complete\n\n\n ");
}

static void anti_aliasing_filter_Y(uint8_t* input, uint8_t* output, int Height, int Width, int filter_type, int size, double sigma)
{
	int i,j,k,l;  
	int middle = (size - size%2)/2;
	double *gaussian_filter=NULL;
	double **average_filter=NULL;
	double sum_filter;
	double PI = 3.14159265358979; 
	double value;

	/*initialize filters*/
	if(filter_type==1)
	{     
		average_filter = (double**)calloc(size,sizeof(double*));   
		if(!average_filter)
		{
			fprintf(stderr,"gaussian_filter buffer allocation error\n");
			exit(-1);
		}

		(average_filter)[0] = (double* )calloc(size*size,sizeof(double ));    
		if(!(average_filter)[0])
		{
		fprintf(stderr,"gaussian_filter buffer allocation error\n");
		exit(-1);
		}
      
		for(i=1 ; i<size ; i++)
			(average_filter)[i] =  (average_filter)[i-1] + size;

		for (i=0;i<size;i++)
			for (j=0;j<size;j++)
				average_filter[i][j]=1;	  

		sum_filter=size*size;	
	}

	if (filter_type==2)
    {      
		gaussian_filter = (double*)calloc(size,sizeof(double));  
		if(!gaussian_filter)
		{
			fprintf(stderr,"gaussian_filter buffer allocation error\n");
			exit(-1);
		}

		sum_filter=0;  
		for (i = 0; i<size;i++)  
			sum_filter+=gaussian_filter[i]= (1/(double)(sqrt(2*PI)*sigma))* (exp(- (double)((i-middle)*(i-middle)/(double)(2*sigma*sigma))));
	}
	  	  

  
	/*no filter*/
	if (filter_type == 0)
	{    
		/*replication of the border*/
		for(i=0;i<Height;i++)        
			for(j=0;j<Width;j++) 	    	     
				output[i*Width+j] = input[i*Width+j];
	}
  
	/*average filter*/
	if (filter_type == 1)
	{      
		/*replication of the border*/
		for(i=0;i<Height;i++)        
		for(j=0;j<Width;j++) 	    	     
		output[i*Width+j] = input[i*Width+j];
	          
		for(i=middle;i<Height-middle;i++) 
		{ 
			for(j=middle;j<Width-middle;j++) 
			{  
				value=0;
				for (k=0;k<size;k++)
					for(l=0;l<size;l++)	      
						value += average_filter[k][l]*input[(i-middle+k)*(Width-middle)+j-middle+l] ;	  
				output[i*(Width-middle)+j]=(unsigned char)(value/sum_filter);     
			}
		}
    }
  
	/*gaussian filter (separated into horizontal and vertical filtering)*/
	if (filter_type==2)
	{      
		/*replication of the border*/
		for(i=0;i<Height;i++)     
		for(j=0;j<Width;j++) 	 
			output[i*Width+j] = input[i*Width+j];
  
		for(i=middle;i<Height-middle;i++)    
		for(j=0;j<Width;j++) 
		{  
			value=0;
			for (k=0;k<size;k++)			
			value += gaussian_filter[k]*input[(i-middle+k)*Width+j] ;	
		  
			output[i*Width+j]=(unsigned char)(value/sum_filter); 
	    
		}   

		/*copy the filter picture into the original for the second filtering*/
		for(i=0;i<Height;i++)     
			for(j=0;j<Width;j++) 		
				input[i*Width+j] = output[i*Width+j];
    
		for(i=0;i<Height;i++) 
			for(j=middle;j<Width-middle;j++) 
			{  
				value=0;
				for (k=0;k<size;k++)			
					value += gaussian_filter[k]*input[i*(Width-middle)+j-middle+k] ;
				    //value += gaussian_filter[k]*input[i*Width+j-middle+k] ;	
	  
				/*filtered_frame*/
				output[i*(Width-middle)+j]=(unsigned char)(value/sum_filter);
				//output[i*Width+j]=(unsigned char)(value/sum_filter);
   			}
    }

	if (filter_type==1)
	{          
		free (average_filter[0]);      
		free (average_filter);
	}

	if (filter_type==2)
       free(gaussian_filter);
}

/*
************************
*
* idem for chroma
*
************************
*/
static void anti_aliasing_filter_UV(uint8_t* input, uint8_t *output, int Height, int Width, int filter_type, int size, double sigma)
{
  int i,j,k,l;  
  int middle = (size - size%2)/2;
  double *gaussian_filter;
  double **average_filter;
  double sum_filter;
  double PI = 3.14159265358979; 
  double value;

  /*initialize filters  */
  if(filter_type==1)
  {     
		average_filter = (double**)calloc(size,sizeof(double*));   
		if(!average_filter)
		{
			fprintf(stderr,"gaussian_filter buffer allocation error\n");
			exit(-1);
		}
		(average_filter)[0] = (double* )calloc(size*size,sizeof(double ));    
		if(!(average_filter)[0])
		{
			fprintf(stderr,"gaussian_filter buffer allocation error\n");
			exit(-1);
		}

      for(i=1 ; i<size ; i++)
		(average_filter)[i] =  (average_filter)[i-1] + size;

      for (i=0;i<size;i++)
		for (j=0;j<size;j++)
		average_filter[i][j]=1;	  

      sum_filter=size*size;	
    }

	if (filter_type==2)
	{      
		gaussian_filter = (double*)calloc(size,sizeof(double));  
		if(!gaussian_filter)
		{
			fprintf(stderr,"gaussian_filter buffer allocation error\n");
			exit(-1);
		}

      
		sum_filter=0; 
		for (i = 0; i<size;i++)  
			sum_filter+=gaussian_filter[i]= (1/(double)(sqrt(2*PI)*sigma))* (exp(- (double)((i-middle)*(i-middle)/(double)(2*sigma*sigma))));
	}
  
  
	/*no filter*/
	if (filter_type == 0)
	{    
		/*replication of the border*/
		for(i=0;i<Height;i++)        
			for(j=0;j<Width;j++) 	    	     
				output[i*Width+j]=input[i*Width+j];  
	}
  
	/*average filter*/
	if (filter_type == 1)
	{      
		/*replication of the border*/
		for(i=0;i<Height;i++)        
			for(j=0;j<Width;j++) 	    	     
				output[i*Width+j]=input[i*Width+j]; 
	          
		for(i=middle;i<Height-middle;i++) 
		{ 
			for(j=middle;j<Width-middle;j++) 
			{  
				value=0;
				for (k=0;k<size;k++)
					for(l=0;l<size;l++)	      
						value += average_filter[k][l]*input[(i-middle+k)*(Width-middle)+j-middle+l] ;	      
		  
				output[i*(Width-middle)+j]=(unsigned char)(value/sum_filter);     
			}
		}
    }
  
	/*gaussian filter (separated into horizontal and vertical filtering)*/
	if (filter_type==2)
	{      
		/*replication of the border*/
		for(i=0;i<Height;i++)     
			for(j=0;j<Width;j++) 	 
				output[i*Width+j]=input[i*Width+j]; 
  
		for(i=middle;i<Height-middle;i++)    
			for(j=0;j<Width;j++) 
			{  
				value=0;
				for (k=0;k<size;k++)			
					value += gaussian_filter[k]*input[(i-middle+k)*Width+j] ;			
	  
				output[i*Width+j]=(unsigned char)(value/sum_filter); 
			}   

		/*copy the filter picture into the original for the second filtering*/
		for(i=0;i<Height;i++)     
			for(j=0;j<Width;j++) 		
				input[i*Width+j]=output[i*Width+j]; 
	    
		for(i=0;i<Height;i++) 
			for(j=middle;j<Width-middle;j++) 
			{  
				value=0;
				for (k=0;k<size;k++)			
					value += gaussian_filter[k]*input[i*(Width-middle)+j-middle+k] ;			
			  
				output[i*(Width-middle)+j]=(unsigned char)(value/sum_filter); 
			}
	}

	if (filter_type==1)
	{          
		free (average_filter[0]);      
		free (average_filter);
	}
  
	if (filter_type==2)
       free(gaussian_filter);
}

void nearest_neighbour(int input_width, int input_height, int output_width, int output_height, FILE *ft)
{
	/*******nearest neighbour ****/
	double nXFactor = (double)input_width/(double)output_width;
	double nYFactor = (double)input_height/(double)output_height;
	int i,j,u,v;

	int size;
	double sigma;

	size = 3;
	sigma = 1;
  
	if ((input_height==720)&&(input_width==1280))/*HD input*/
    {
		if ((output_height==576)&&(input_width==720))/*output SD*/
		{
			size = 3;
			sigma = 0.7;
		}
		else if ((output_height==288)&&(output_width==352))/*output CIF*/
		{
			size = 5;
			sigma = 2;
		}
		else if ((output_height==144)&&(output_width==176))/*output QCIF*/
		{
			size = 11;
			sigma = 3.4;
		}
    }

	if ((input_height==576)&&(input_width==720))/*SD input*/
	{
		if ((output_height==288)&&(output_width==352))/*output CIF*/
		{
			size = 3;
			sigma = 1;
		}
		else if ((output_height==144)&&(output_width==176))/*output QCIF*/
		{
			size = 9;
			sigma = 1.8;
		}
    }

	if ((input_height==288)&&(input_width==352))/*SD input*/
    {     
      if ((output_height==144)&&(output_width==176))/*output QCIF*/
		{
			size = 3;
			sigma = 1;
		}
    }

	if ((output_height >= input_height)&&(output_width >= input_width))/*in the case of up-sampling don't filter*/
    	anti_aliasing_filter_Y(img, img1, input_height, input_width, 0, 0, 0);
	else /*for down sampling apply an anti aliasing filter*/
    	anti_aliasing_filter_Y(img, img1, input_height, input_width, 2, size, sigma); 

	//memcpy(img,img1, input_height*input_width);
	
    /****y ****/
	for (i = 0; i < output_height; i++)
	{
		for (j = 0 ; j < output_width; j++)
		{
			u = (int)(i*nYFactor);
			v = (int)(j*nXFactor);
			//img1[i*output_width + j] = img[u*input_width + v]; 
			img2[i*output_width + j] = img1[u*input_width + v];
		}
	}
	fwrite(img2,1, output_width*output_height,ft);
	memcpy(img1,img2, input_height*input_width);
	

	/****u ****/
	if ((output_height >= input_height)&&(output_width >= input_width))/*in the case of up-sampling don't filter*/
    	anti_aliasing_filter_UV(img + (input_height*input_width), img1, input_height/2, input_width/2, 0, 0, 0);
	else /*for down sampling apply an anti aliasing filter*/
    	anti_aliasing_filter_UV(img + (input_height*input_width), img1, input_height/2, input_width/2, 2, size, sigma); 
	//memcpy(img + (input_height*input_width),img1, (int)(input_height*input_width*0.25));

	for (i = 0; i < output_height/2; i++)
	{
		for (j = 0 ; j < output_width/2; j++)
		{
			u = (int)(i*nYFactor);
			v = (int)(j*nXFactor);
			//img1[i*output_width/2 + j] = img[(input_width*input_height)+u*input_width/2 + v]; 
			img2[i*output_width/2 + j] = img1[u*input_width/2 + v];
		}
	}
	fwrite(img2,1, output_width/2*output_height/2,ft);
	memcpy(img1 + (input_height*input_width),img2, (input_height*input_width)/4);
		
	/****v ****/
	if ((output_height >= input_height)&&(output_width >= input_width))/*in the case of up-sampling don't filter*/
    	anti_aliasing_filter_UV(img + (int)(input_height*input_width*1.25), img1, input_height/2, input_width/2, 0, 0, 0);
	else /*for down sampling apply an anti aliasing filter*/
    	anti_aliasing_filter_UV(img + (int)(input_height*input_width*1.25), img1, input_height/2, input_width/2, 2, size, sigma);
	//memcpy(img + (int)(input_height*input_height*1.25),img1, (int)(input_height*input_width*0.25));

	for (i = 0; i < output_height/2; i++)
	{
		for (j = 0 ; j < output_width/2; j++)
		{
			u = (int)(i*nYFactor);
			v = (int)(j*nXFactor);
			//img1[i*output_width/2 + j] = img[(int)(input_width*input_height*1.25)+u*input_width/2 + v];           
			img2[i*output_width/2 + j] = img1[u*input_width/2 + v];           
		}
	}
	fwrite(img2,1, output_width/2*output_height/2,ft);
	memcpy(img1 + (int)(input_height*input_width*1.25),img2, (input_height*input_width)/4);

	printf("\n one frame write(nearest neighbour) complete \n");
}


void bilinear(int input_width, int input_height, int output_width, int output_height, FILE *ft)
{
	/*****bilinear***/
	double nXFactor = (double)input_width/(double)output_width;
	double nYFactor = (double)input_height/(double)output_height;
	double fraction_x, fraction_y, one_minus_x, one_minus_y;
    int ceil_x, ceil_y, floor_x, floor_y;
    int c1,c2,c3,c4;
	int b1, b2,i,j, temp;

	int size;
	double sigma;

	size = 3;
	sigma = 1;
  
	if ((input_height==720)&&(input_width==1280))/*HD input*/
    {
		if ((output_height==576)&&(input_width==720))/*output SD*/
		{
			size = 3;
			sigma = 0.7;
		}
		else if ((output_height==288)&&(output_width==352))/*output CIF*/
		{
			size = 5;
			sigma = 2;
		}
		else if ((output_height==144)&&(output_width==176))/*output QCIF*/
		{
			size = 11;
			sigma = 3.4;
		}
    }

	if ((input_height==576)&&(input_width==720))/*SD input*/
	{
		if ((output_height==288)&&(output_width==352))/*output CIF*/
		{
			size = 3;
			sigma = 1;
		}
		else if ((output_height==144)&&(output_width==176))/*output QCIF*/
		{
			size = 9;
			sigma = 1.8;
		}
    }

	if ((input_height==288)&&(input_width==352))/*SD input*/
    {     
      if ((output_height==144)&&(output_width==176))/*output QCIF*/
		{
			size = 3;
			sigma = 1;
		}
    }

	int hw = input_height*input_width;
	int hw1 = (int)(hw*1.25);

	if ((output_height >= input_height)&&(output_width >= input_width))/*in the case of up-sampling don't filter*/
    	anti_aliasing_filter_Y(img, img1, input_height, input_width, 0, 0, 0);
	else /*for down sampling apply an anti aliasing filter*/
    	anti_aliasing_filter_Y(img, img1, input_height, input_width, 2, size, sigma); 

	/***** for y *****/
	for (i = 0; i < output_height; ++i)
	{
		for (j = 0; j < output_width; ++j)
		{
			floor_x = (int)(j * nXFactor);
			floor_y = (int)(i * nYFactor);
			ceil_x = floor_x + 1;
			if (ceil_x >= input_width) ceil_x = floor_x;
			ceil_y = floor_y + 1;
			if (ceil_y >= input_height) ceil_y = floor_y;
			fraction_x = j * nXFactor - floor_x;
			fraction_y = i * nYFactor - floor_y;
			one_minus_x = 1.0 - fraction_x;
			one_minus_y = 1.0 - fraction_y;

			c1 = img1[floor_y*input_width + floor_x];
			c2 = img1[floor_y*input_width + ceil_x];
			c3 = img1[ceil_y*input_width + floor_x];
			c4 = img1[ceil_y*input_width + ceil_x];

			b1 = (int)(one_minus_x * c1 + fraction_x * c2);
			b2 = (int)(one_minus_x * c3 + fraction_x * c4);
			temp = (int)(one_minus_y * (double)(b1) + fraction_y * (double)(b2));
			if(temp>255)
				temp = 255;
			if(temp<0)
				temp = 0;

			img2[i*output_width +j] = temp; 
		}
	}
	fwrite(img2,1,(output_width*output_height),ft);
	memcpy(img1, img2, output_height*output_width);


	/****for u ****/
	if ((output_height >= input_height)&&(output_width >= input_width))/*in the case of up-sampling don't filter*/
    	anti_aliasing_filter_UV(img + hw, img1, input_height/2, input_width/2, 0, 0, 0);
	else /*for down sampling apply an anti aliasing filter*/
    	anti_aliasing_filter_UV(img + hw, img1, input_height/2, input_width/2, 2, size, sigma); 
	
	for (i = 0; i < output_height/2; ++i)
	{
		for (j = 0; j < output_width/2; ++j)
		{
			floor_x = (int)(j * nXFactor);
			floor_y = (int)(i * nYFactor);
			ceil_x = floor_x + 1;
			if (ceil_x >= input_width/2) ceil_x = floor_x;
			ceil_y = floor_y + 1;
			if (ceil_y >= input_height/2) ceil_y = floor_y;
			fraction_x = j * nXFactor - floor_x;
			fraction_y = i * nYFactor - floor_y;
			one_minus_x = 1.0 - fraction_x;
			one_minus_y = 1.0 - fraction_y;

			c1 = img1[floor_y*input_width/2 + floor_x];
			c2 = img1[floor_y*input_width/2 + ceil_x];
			c3 = img1[ceil_y*input_width/2 + floor_x];
			c4 = img1[ceil_y*input_width/2 + ceil_x];

			b1 = (int)(one_minus_x * c1 + fraction_x * c2);
			b2 = (int)(one_minus_x * c3 + fraction_x * c4);
			temp = (int)(one_minus_y * (double)(b1) + fraction_y * (double)(b2));
			if(temp>255)
				temp = 255;
			if(temp<0)
				temp = 0;

			img2[i*output_width/2 +j] = temp; 
		}
	}
	fwrite(img2,1,(output_width*output_height)/4,ft);
	memcpy(img1 + (output_height*output_width), img2, (output_height*output_width)/4);
	
	/****for v ****/
	if ((output_height >= input_height)&&(output_width >= input_width))/*in the case of up-sampling don't filter*/
    	anti_aliasing_filter_UV(img + hw1, img1, input_height/2, input_width/2, 0, 0, 0);
	else /*for down sampling apply an anti aliasing filter*/
    	anti_aliasing_filter_UV(img + hw1, img1, input_height/2, input_width/2, 2, size, sigma);
	for (i = 0; i < output_height/2; ++i)
	{
		for (j = 0; j < output_width/2; ++j)
		{
			floor_x = (int)(j * nXFactor);
			floor_y = (int)(i * nYFactor);
			ceil_x = floor_x + 1;
			if (ceil_x >= input_width/2) ceil_x = floor_x;
			ceil_y = floor_y + 1;
			if (ceil_y >= input_height/2) ceil_y = floor_y;
			fraction_x = j * nXFactor - floor_x;
			fraction_y = i * nYFactor - floor_y;
			one_minus_x = 1.0 - fraction_x;
			one_minus_y = 1.0 - fraction_y;

			c1 = img1[floor_y*input_width/2 + floor_x];
			c2 = img1[floor_y*input_width/2 + ceil_x];
			c3 = img1[ceil_y*input_width/2 + floor_x];
			c4 = img1[ceil_y*input_width/2 + ceil_x];

			b1 = (int)(one_minus_x * c1 + fraction_x * c2);
			b2 = (int)(one_minus_x * c3 + fraction_x * c4);
			temp = (int)(one_minus_y * (double)(b1) + fraction_y * (double)(b2));
			if(temp>255)
				temp = 255;
			if(temp<0)
				temp = 0;

			img2[i*output_width/2 +j] = temp; 
		}
	}
	fwrite(img2,1,(output_width*output_height)/4,ft);
	memcpy(img1 + (int)(output_height*output_width*1.25), img2, (output_height*output_width)/4);

	printf("\none frame(bilinear)write complete ");
}

int _tmain(int argc, _TCHAR* argv[])
{
	int choice, frame_num; 
	int output_width, output_height, input_width, input_height;
	char input_file[256],scaled_image_file[256],reconstructed_image_file[256];
	FILE *ft1, *ft2;
	
	/***** take i/p *****/
	/*****for input resolution ****/
	printf(" enter yuv/ Y  file :     ");
	scanf("%s",input_file);

	FILE *fp = fopen(input_file,"rb");
	if(fp == NULL)
	{
		printf("\ncannot open input file\n");
		exit(1);
	}

	printf("\n enter resolution of input picture (should be multiple of 16)  :\n width   =  ");
	scanf("%d", &input_width);
		
	printf(" height  =  ");
	scanf("%d",&input_height);
		
	
	/*****for output resolution ****/
	printf("\n\n enter resolution of target picture (should be multiple of 16)  :\n width   =  ");
	scanf("%d", &output_width);
	output_width = (output_width>>4)<<4;

	printf(" height  =  ");
	scanf("%d", &output_height);
	output_height = (output_height>>4)<<4;

	printf("accepted width_output = %d and height_input = %d\n", output_width, output_height);

	printf("\nenter filter number:\n1.nearest neighbour\n2.bilinear\n3.polyphase\n\n\n  :");
	scanf("%d",&choice);

	printf("\nenter frame number :  ");
	scanf("%d",&frame_num);

	if(choice == 1)
	{
		strcpy(scaled_image_file, "nearest_neighbour_scaled.yuv");
		strcpy(reconstructed_image_file, "nearest_neighbour_reconstructed.yuv");
	}
	else if(choice == 2)
	{
		strcpy(scaled_image_file, "bilinear_scaled.yuv");
		strcpy(reconstructed_image_file, "bilinear_reconstructed.yuv");
	}
	else if(choice == 3)
	{
		strcpy(scaled_image_file, "polyphase_scaled.yuv");
		strcpy(reconstructed_image_file, "polyphase_reconstructed.yuv");

	}
	ft1 = fopen(scaled_image_file, "wb");


	/**** scaling ****/
	/**** creating scaled image ****/
	for(int a = 0; a < frame_num; a++)
	{
		fread(img,1,(int)(input_width*input_height*1.5),fp);
		switch(choice)
		{
		case 1:
			nearest_neighbour(input_width, input_height, output_width, output_height, ft1);
			break;
		case 2:
			bilinear(input_width, input_height, output_width, output_height, ft1);
			break;
		case 3:
			polyphase(input_width, input_height, output_width, output_height, ft1);
			break;
		default:
			printf("\n not a valid filter...\n\n..exiting...");
			break;
		}
	}
	/***** close original file ****/
	fclose(fp);
	fclose(ft1);
	printf("\n\n\n scaled file write complete.....\n\n");

	/**reconstruction ***/
	ft1 = fopen(scaled_image_file, "rb");
	ft2 = fopen(reconstructed_image_file, "wb");
	for(int a = 0; a < frame_num; a++)
	{
		fread(img,1,(int)(input_width*input_height*1.5),ft1);
		switch(choice)
		{
		case 1:
            nearest_neighbour(output_width, output_height,input_width, input_height, ft2);
			break;
		case 2:
			bilinear(output_width, output_height, input_width, input_height, ft2);
			break;
		case 3:
			polyphase(output_width, output_height, input_width, input_height, ft2);
			break;
		default:
			printf("\n not a valid filter...\n\n..exiting...");
			break;
		}
	}
	/**** scaled file read  and reconstruction file write complete ***/
	fclose(ft1);
	fclose(ft2);
	printf("\n\n\n reconstructed file write complete.....\n\n");
	return 0;
}
	


