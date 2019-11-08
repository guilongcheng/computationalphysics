#ifndef GLC
#define GLC
#define PI 3.14159265358979323846264338327950288419716939937510
#define NMAX 21 //拉格朗日插值函数插值点的最大个数
#define NMAXITER 1000 //自适应微分和积分的最大迭代次数
#define MMAX 100 //求解微分方程 微分方程组的最大个数
#define NPARS 20 //参数数组的大小
#define MALLOC(n, type) ((type *) malloc ( (n) * sizeof(type)))
#define SAVEFILE
void savef(char filename[], double *data,int n, int m, int showm);
#endif