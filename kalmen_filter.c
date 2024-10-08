#include "kalmen_filter.h"

float dt =0.008;//����ʱ��
float Q_angle=0.01,Q_gyr=0.005;// Ԥ��ģ�͵ķ���
float R_angle=2,R_gyr=0.9;//�۲�ģ�͵ķ���
struct F F1={0,0,0};  //�������
struct P P1={1,0,0.2,0,0,0.2};//�������  ǰ����Ϊ����Ƕ�����ٶȣ�����λΪ����Э�������
float k_det=0,k[4];//�����������������ʽ�뿨�����������
/*
 *  name:   �������˲�����
 *  para:   �����ĽǶȣ������������ٶ�
 *  return: �˲��Ƕ�
 *  writer: YXL
 *  function:���׿������˲�
 *  time:   2022/1/22
 */
float kalmen_filter(float angle_m ,float gyr_m)
{
    float d_angle,d_gyr;//����Ԥ������������Ĳ�ֵ
    //////////////////Ԥ�ⲽ/////////////
    ///////����Ԥ��
    F1.angle=P1.angle +dt*P1.gyr;
    F1.gyr=P1.gyr;
    ///////Э����Ԥ��
    F1.cov[0]=P1.cov[0]+dt*(P1.cov[1]+P1.cov[2])+Q_angle;
    F1.cov[1]=P1.cov[1]+dt*P1.cov[3];
    F1.cov[2]=P1.cov[2]+dt*P1.cov[3];
    F1.cov[3]=P1.cov[3]+Q_gyr;
    /////////////////���²�////////////////
    //���� �������������
    k_det=1/((F1.cov[0]+R_angle)*(F1.cov[3]+R_gyr)-F1.cov[1]*F1.cov[2]);
    k[0]=k_det*(F1.cov[0]*(F1.cov[3]+R_gyr)-F1.cov[1]*F1.cov[2]);
    k[1]=k_det*(F1.cov[1]*(F1.cov[0]+R_angle)-F1.cov[0]*F1.cov[1]);
    k[2]=k_det*(F1.cov[2]*(F1.cov[3]+R_gyr)-F1.cov[2]*F1.cov[3]);
    k[3]=k_det*(F1.cov[3]*(F1.cov[0]+R_angle)-F1.cov[1]*F1.cov[2]);
    //���������������ɺ���������
    d_angle=angle_m-F1.angle;
    d_gyr=gyr_m-F1.gyr;
    P1.angle=F1.angle+k[0]*d_angle+k[1]*d_gyr;
    P1.gyr=F1.gyr+k[2]*d_angle+k[3]*d_gyr;    
    //����Э������󣨼�����Э�������
    P1.cov[0]=F1.cov[0]*(1-k_det*k[0])-k_det*k[1]*F1.cov[2];
    P1.cov[1]=F1.cov[1]*(1-k_det*k[0])-k_det*k[1]*F1.cov[3];
    P1.cov[2]=F1.cov[2]*(1-k_det*k[3])-k_det*k[2]*F1.cov[0];
    P1.cov[3]=F1.cov[3]*(1-k_det*k[3])-k_det*k[2]*F1.cov[1];
    
    //��ɵ��������������Ƕ�
    return P1.angle;
}

