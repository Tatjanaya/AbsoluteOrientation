#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define N 6  // 控制点数量
#define PI 3.1415926
#define MAXITERARATION 8 //最大迭代次数


void MatrixTranspose(double* MatrixOrigin, double* MatrixNew, int m, int n)
{
	int i;
	int j;
	for (i = 0; i != n; i++)
	{
		for (j = 0; j != m; j++)
		{
			MatrixNew[i*m + j] = MatrixOrigin[j*n + i];
		}
	}
};

void MatrixInversion(double* Matrix, int m)
{
	int i, j, k;
	for (k = 0; k != m; k++)
	{
		for (i = 0; i != m; i++)
		{
			if (i != k)
				Matrix[i*m + k] = -Matrix[i*m + k] / Matrix[k*m + k];
		}
		Matrix[k*m + k] = 1 / Matrix[k*m + k];

		for (i = 0; i != m; i++)
		{
			if (i != k)
			{
				for (j = 0; j != m; j++)
				{
					if (j != k)
						Matrix[i*m + j] += Matrix[k*m + j] * Matrix[i*m + k];
				}
			}
		}

		for (j = 0; j != m; j++)
		{
			if (j != k)
				Matrix[k*m + j] *= Matrix[k*m + k];
		}
	}
};

void MatrixMultiply(double* MatrixA, double* MatrixB, double* MatrixC, int ARow, int AColumn, int BColumn)
{
	int i;
	int j;
	int k;

	for (i = 0; i != ARow; i++)
	for (j = 0; j<BColumn; j++)
	{
		MatrixC[i*BColumn + j] = 0.0;
		for (k = 0; k<AColumn; k++)
			MatrixC[i*BColumn + j] += MatrixA[i*AColumn + k] * MatrixB[j + k*BColumn];
	}
};

void MatrixAdd(double* MatrixA, double* MatrixB, double* MatrixC, int Row, int Column)
{
	int i;
	int j;
	for (i = 0; i != Row; i++)
	{
		for (j = 0; j != Column; j++)
		{
			MatrixC[i*Column + j] = MatrixA[i*Column + j] + MatrixB[i*Column + j];
		}
	}
};


void MatrixAdd(double* MatrixA, double* MatrixB, double* MatrixC, int length)
{
	int i;
	for (i = 0; i != length; i++)
	{
		MatrixC[i] = MatrixA[i] + MatrixB[i];
	}
};

void MatrixMinus(double* MatrixA, double* MatrixB, double* MatrixC, int length)
{
	int i;
	for (i = 0; i != length; i++)
	{
		MatrixC[i] = MatrixA[i] - MatrixB[i];
	}
};

void MatrixCopy(double* MatrixA, double* MatrixB, int length)
{
	int i;
	for (i = 0; i != length; i++)
	{
		MatrixB[i] = MatrixA[i];
	}
};

//检查改正数是否已经达到精度要求
int CheckPrecision(double* deta)
{
	int ret;
	ret = (fabs(deta[0])<0.5&&fabs(deta[1])<0.5&&fabs(deta[2])<0.5&&
		fabs(deta[4])<0.5 / 60 * PI / 180 && fabs(deta[5])<0.5 / 60 * PI / 180 && fabs(deta[6])<
		0.5 / 60 * PI / 180);
	return ret;
}


struct SourceData    //保存原始数据
{
	double x;
	double y;
	double z;
	double X;
	double Y;
	double Z;
};

int main()
{
	int dh[6], a;
	double x0, y0, z0;
	double X00, Y00, Z00;
	SourceData sd[10];
	char *p = "模型点-控制点-坐标.txt";
	FILE* xp = fopen(p, "r");

	//读入控制点及模型点坐标
	int i = 0;
	
	sd[0].x = -2.994926; sd[0].y = 98.313214; sd[0].z = -165.370335; sd[0].X = 27313.512000; sd[0].Y = 2700167.702000; sd[0].Z = 103.950000;
	sd[1].x = 115.300090; sd[1].y = 106.807568; sd[1].z = -166.986144; sd[1].X = 28500.938000; sd[1].Y = 2700184.416000; sd[1].Z = 97.350000;
	sd[2].x = -10.104023; sd[2].y = -76.494059; sd[2].z = -165.102793; sd[2].X = 27141.968000; sd[2].Y = 2698422.955000; sd[2].Z = 101.994000;
	sd[3].x = 116.937501; sd[3].y = -79.779735; sd[3].z = -162.042707; sd[3].X = 28409.856000; sd[3].Y = 2698319.640000; sd[3].Z = 155.804000;
	sd[4].x = -19.486363; sd[4].y = 13.056943; sd[4].z = -160.562998; sd[4].X = 27102.439000; sd[4].Y = 2699324.440000; sd[4].Z = 163.290000;
	sd[5].x = 90.631173; sd[5].y = 7.206584; sd[5].z = -166.162713; sd[5].X = 28197.742000; sd[5].Y = 2699201.833000; sd[5].Z = 100.000000;
	fclose(xp);

	double lamd = 1.0, phi, omega, kappa, Xs, Ys, Zs, Xtpg = 0, Ytpg = 0, Ztpg = 0, xg = 0, yg = 0, zg = 0;
	double lamd0 = 1.0;
	phi = omega = kappa = 0.0;
	Xs = Ys = Zs = 0.0;

	//计算重心化坐标
	int n = 0;
	for (i = 0; i<6; i++)
	{
		Xtpg += sd[i].X;
		Ytpg += sd[i].Y;
		Ztpg += sd[i].Z;
		n++;
	}
	Xtpg = Xtpg / n;
	Ytpg = Ytpg / n;
	Ztpg = Ztpg / n;
	double Xtp_[6] = { 0 }, Ytp_[6] = { 0 }, Ztp_[6] = { 0 };
	for (i = 0; i<6; i++)
	{
		Xtp_[i] = sd[i].X - Xtpg;
		Ytp_[i] = sd[i].Y - Ytpg;
		Ztp_[i] = sd[i].Z - Ztpg;
	}
	double x_[6] = { 0 }, y_[6] = { 0 }, z_[6] = { 0 };

	/*声明旋转矩阵*/
	double R[9];
	double lX[6], lY[6], lZ[6], L[3];

	double A[3 * 7] = { 0.0 };
	double AT[7 * 3] = { 0.0 };

	double deta[7] = { 1, 1, 1, 1, 1, 1, 1 };

	double Buf1[49] = { 0.0 };
	double Buf2[49] = { 0.0 };//ATA累加

	double Buf3[7] = { 0.0 };
	double Buf4[7] = { 0.0 };//ATL累加

	double Buf5[3 * 6 * 7] = { 0.0 };//存储3*N*6的A矩阵
	double Buf6[3 * 6 * 1] = { 0.0 };//存储3*N*1的L矩阵
	double V[3 * 6 * 1] = { 0.0 };

	double ATA[49] = { 0.0 };
	double ATL[7] = { 0.0 };

	int iCount = 0;//声明迭代次数
	printf("开始迭代计算…………\n\n\n");
	while (!CheckPrecision(deta))
	{
		printf("第%d次迭代：", ++iCount);

		if (iCount == MAXITERARATION)
		{
			printf("迭代次数超限，可能不收敛\n");
			break;
		}

		//每次迭代之前必须清空两个保存累加值的矩阵ATA与ATL
		for (int i = 0; i != 49; i++)
		{
			ATA[i] = 0.0;
			if (i<7)
			{
				ATL[i] = 0.0;
			}
		}

		int n = 0;
		for (i = 0; i<6; i++)
		{
			xg += sd[i].x;
			yg += sd[i].y;
			zg += sd[i].z;
			n++;
		}
		xg = xg / n;
		yg = yg / n;
		zg = zg / n;

		for (i = 0; i<6; i++)
		{
			x_[i] = sd[i].x - xg;
			y_[i] = sd[i].y - yg;
			z_[i] = sd[i].z - zg;
		}


		/////////////////////////////计算常数项
		//计算旋转矩阵R
		R[0] = cos(phi)*cos(kappa) - sin(phi)*sin(omega)*sin(kappa);
		R[1] = -cos(phi)*sin(kappa) - sin(phi)*sin(omega)*cos(kappa);
		R[2] = -sin(phi)*cos(omega);
		R[3] = cos(omega)*sin(kappa);
		R[4] = cos(omega)*cos(kappa);
		R[5] = -sin(omega);
		R[6] = sin(phi)*cos(kappa) + cos(phi)*sin(omega)*sin(kappa);
		R[7] = -sin(phi)*sin(kappa) + cos(phi)*sin(omega)*cos(kappa);
		R[8] = cos(phi)*cos(omega);

		for (i = 0; i<6; i++)
		{
			lX[i] = Xtp_[i] - lamd*(R[0] * x_[i] + R[1] * y_[i] + R[2] * z_[i]) - Xs;
			lY[i] = Ytp_[i] - lamd*(R[3] * x_[i] + R[4] * y_[i] + R[5] * z_[i]) - Ys;
			lZ[i] = Ztp_[i] - lamd*(R[6] * x_[i] + R[7] * y_[i] + R[8] * z_[i]) - Zs;

			A[0] = 1;
			A[1] = 0;
			A[2] = 0;
			A[3] = (R[0] * x_[i] + R[1] * y_[i] + R[2] * z_[i]);;
			A[4] = -lamd*(R[6] * x_[i] + R[7] * y_[i] + R[8] * z_[i]);
			A[5] = -lamd*(R[3] * x_[i] + R[4] * y_[i] + R[5] * z_[i])*sin(phi);
			A[6] = -lamd*(R[3] * x_[i] + R[4] * y_[i] + R[5] * z_[i])*cos(phi)*cos(omega) - lamd*(R[6] * x_[i] + R[7] * y_[i] + R[8] * z_[i])*sin(omega);
			A[7] = 0;
			A[8] = 1;
			A[9] = 0;
			A[10] = (R[3] * x_[i] + R[4] * y_[i] + R[5] * z_[i]);
			A[11] = 0;
			A[12] = lamd*(R[0] * x_[i] + R[1] * y_[i] + R[2] * z_[i])*sin(phi) - lamd*(R[6] * x_[i] + R[7] * y_[i] + R[8] * z_[i])*cos(phi);
			A[13] = lamd*(R[0] * x_[i] + R[1] * y_[i] + R[2] * z_[i])*cos(phi)*cos(omega) + lamd*(R[6] * x_[i] + R[7] * y_[i] + R[8] * z_[i])*sin(phi)*cos(omega);
			A[14] = 0;
			A[15] = 0;
			A[16] = 1;
			A[17] = (R[6] * x_[i] + R[7] * y_[i] + R[8] * z_[i]);
			A[18] = lamd*(R[0] * x_[i] + R[1] * y_[i] + R[2] * z_[i]);
			A[19] = lamd*(R[3] * x_[i] + R[4] * y_[i] + R[5] * z_[i])*cos(phi);
			A[20] = lamd*(R[0] * x_[i] + R[1] * y_[i] + R[2] * z_[i])*sin(omega) - lamd*(R[3] * x_[i] + R[4] * y_[i] + R[5] * z_[i])*sin(phi)*cos(omega);
			//该循环保存A矩阵，最后评定精度
			for (int l = 0; l<21; l++)
			{
				Buf5[21 * i + l] = A[l];
			}

			//所谓逐步法，即要在循环内部将ATA计算出来累加，下面的L矩阵类似
			MatrixTranspose(A, AT, 3, 7);
			MatrixMultiply(AT, A, Buf1, 7, 3, 7);
			MatrixCopy(ATA, Buf2, 49);
			MatrixAdd(Buf1, Buf2, ATA, 49);//为逐步法化后的ATA矩阵累加

			L[0] = lX[i];
			L[1] = lY[i];
			L[2] = lZ[i];

			//保存L矩阵，最后评定精度用
			for (int l = 0; l<3; l++)
			{
				Buf6[3 * i + l] = L[l];
			}

			MatrixMultiply(AT, L, Buf3, 7, 3, 1);
			MatrixCopy(ATL, Buf4, 7);
			MatrixAdd(Buf3, Buf4, ATL, 7);
		}

		//"逐步法"的另一处不同，出循环即可直接计算ATA逆乘ATL
		MatrixInversion(ATA, 7);
		MatrixMultiply(ATA, ATL, deta, 7, 7, 1);

		//deta即为改正数
		Xs += deta[0];
		Ys += deta[1];
		Zs += deta[2];
		lamd = lamd + deta[3];
		phi += deta[4];
		omega += deta[5];
		kappa += deta[6];

		printf("改正数值为：\n");
		for (i = 0; i != 7; i++)
		{
			printf("deta[%d]=%lf\n", i, deta[i]);
		}
	}
	double DX, DY, DZ;
	DX = Xtpg - xg;
	DY = Ytpg - yg;
	DZ = Ztpg - zg;

	printf("正常退出迭代。\n\n\n");

	//精度评定
	double Q[7] = { 0.0 };
	for (int h = 0; h<7; h++)
	{
		Q[h] = ATA[h * 7 + h];
	}

	MatrixMultiply(Buf5, deta, V, 3 * 7, 7, 1);//V=Ax-L
	MatrixMinus(V, Buf6, V, 3 * 6);

	double m0 = 0;//单位权中误差
	double VSum = 0.0;//[vv],即平方和


	for (i = 0; i<3 * 6; i++)
	{
		VSum += V[i] * V[i];
	}
	m0 = sqrt(VSum / (3 * 6 - 7));//中误差m0

	double M[7] = { 0.0 };//保存7个值的中误差
	for (i = 0; i != 7; i++)
	{
		M[i] = m0*sqrt(Q[i]);

	}
	printf("计算结果：\n\n\n");

	printf("六个外方位元素为：\n");
	printf("Xs=%.4lf\n", DX);
	printf("Ys=%.4lf\n", DY);
	printf("Zs=%.4lf\n", DZ);
	printf("lamd=%.4lf\n", lamd);
	printf("phi=%.10lf\n", phi);
	printf("omega=%.10lf\n", omega);
	printf("kappa=%.10lf\n\n", kappa);

	printf("旋转矩阵：\n");

	printf("%lf\t%lf\t%lf\n", R[0], R[1], R[2]);
	printf("%lf\t%lf\t%lf\n", R[3], R[4], R[5]);
	printf("%lf\t%lf\t%lf\n\n", R[6], R[7], R[8]);

	printf("单位权中误差：%lf\n\n", m0);

	printf("六个外方位元素的精度：\n");

	printf("lamd:%lf \n", M[3]);
	printf("phi:%lf 秒\n", M[4]);
	printf("omega:%lf 秒\n", M[5]);
	printf("kappa:%lf 秒\n", M[6]);
	system("pause");
	return 0;
}

