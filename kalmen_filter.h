#ifndef __KF_H
#define __KF_H		
#include "sys.h"	 
#include "stdlib.h"



struct F 
{
   float angle;
   float gyr;
   float cov[4];//����ģ��Э����
};
struct P 
{
   float angle;
   float gyr;
   float cov[4];//����ģ��Э����
};

extern struct P P1;//�������
float kalmen_filter(float ,float );


#endif




