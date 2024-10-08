#include "kalmen_filter.h"

float dt =0.008;//采样时间
float Q_angle=0.01,Q_gyr=0.005;// 预测模型的方差
float R_angle=2,R_gyr=0.9;//观测模型的方差
struct F F1={0,0,0};  //先验变量
struct P P1={1,0,0.2,0,0,0.2};//后验变量  前两个为后验角度与角速度，后四位为后验协方差矩阵
float k_det=0,k[4];//卡尔曼增益矩阵行列式与卡尔曼增益矩阵
/*
 *  name:   卡尔曼滤波函数
 *  para:   测量的角度，与该轴测量角速度
 *  return: 滤波角度
 *  writer: YXL
 *  function:二阶卡尔曼滤波
 *  time:   2022/1/22
 */
float kalmen_filter(float angle_m ,float gyr_m)
{
    float d_angle,d_gyr;//先验预测量与测量量的差值
    //////////////////预测步/////////////
    ///////期望预测
    F1.angle=P1.angle +dt*P1.gyr;
    F1.gyr=P1.gyr;
    ///////协方差预测
    F1.cov[0]=P1.cov[0]+dt*(P1.cov[1]+P1.cov[2])+Q_angle;
    F1.cov[1]=P1.cov[1]+dt*P1.cov[3];
    F1.cov[2]=P1.cov[2]+dt*P1.cov[3];
    F1.cov[3]=P1.cov[3]+Q_gyr;
    /////////////////更新步////////////////
    //计算 卡尔曼增益矩阵
    k_det=1/((F1.cov[0]+R_angle)*(F1.cov[3]+R_gyr)-F1.cov[1]*F1.cov[2]);
    k[0]=k_det*(F1.cov[0]*(F1.cov[3]+R_gyr)-F1.cov[1]*F1.cov[2]);
    k[1]=k_det*(F1.cov[1]*(F1.cov[0]+R_angle)-F1.cov[0]*F1.cov[1]);
    k[2]=k_det*(F1.cov[2]*(F1.cov[3]+R_gyr)-F1.cov[2]*F1.cov[3]);
    k[3]=k_det*(F1.cov[3]*(F1.cov[0]+R_angle)-F1.cov[1]*F1.cov[2]);
    //更新期望（即生成后验期望）
    d_angle=angle_m-F1.angle;
    d_gyr=gyr_m-F1.gyr;
    P1.angle=F1.angle+k[0]*d_angle+k[1]*d_gyr;
    P1.gyr=F1.gyr+k[2]*d_angle+k[3]*d_gyr;    
    //更新协方差矩阵（即后验协方差矩阵）
    P1.cov[0]=F1.cov[0]*(1-k_det*k[0])-k_det*k[1]*F1.cov[2];
    P1.cov[1]=F1.cov[1]*(1-k_det*k[0])-k_det*k[1]*F1.cov[3];
    P1.cov[2]=F1.cov[2]*(1-k_det*k[3])-k_det*k[2]*F1.cov[0];
    P1.cov[3]=F1.cov[3]*(1-k_det*k[3])-k_det*k[2]*F1.cov[1];
    
    //完成迭代，返回期望角度
    return P1.angle;
}

