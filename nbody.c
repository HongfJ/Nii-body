#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nbody.h"

#include <time.h>


const double tpi = 2.0 * 3.1415926535897932;
const int nbod = 2; // 天体数量
const int ndim = 6; // 积分维度

double xmass[MAX_BODIES]; // 质量（地球质量单位）
// 新增：观测距离参数（单位：秒差距 pc）
const double observer_distance_pc = 3.0;

// RKF78积分器
int rkf78_no_t(double* ptr_t, double* ptr_dt, double err_tol, double* x, int n);

// 角度转换
double zatan(double sina, double cosa);

// 轨道要素转直角坐标
int ele2rat(double a, double e, double ci, double w, double omega, double cm0, int ii, double xout[6]);

// 直角坐标转轨道要素
int rat2ele(double xin[6], double* a, double* e, double* ci, double* w, double* omega, double* cm, int ii, int ito);

int main()
{
    int i, j;

    // 恒星质量（太阳质量）
    xmass[0] = 1.04;

    clock_t start_total, end_total;
    double total_time_used;
    start_total = clock();

    clock_t start_integrate, end_integrate;
    double integrate_time_used = 0.0; 

    /* 读取初始轨道要素和质量 */
    double p0[MAX_BODIES], e0[MAX_BODIES], ci0[MAX_BODIES], w0[MAX_BODIES];
    double omg0[MAX_BODIES], cm00[MAX_BODIES];
    double raw_mass[MAX_BODIES];     /* 临时存放质量 */

    FILE* fele = fopen("c_ele1.txt", "r");
    if (!fele) { perror("c_ele1.txt"); return 1; }
    // 只读取1颗行星的数据
    for (i = 1; i < nbod; i++)
    {
        fscanf(fele, "%lf %lf %lf %lf %lf %lf %lf",
            &p0[i], &e0[i], &ci0[i], &w0[i], &omg0[i], &cm00[i], &raw_mass[i]);
        xmass[i] = raw_mass[i] * 3e-6;
    }
    fclose(fele);

    FILE* faei = fopen("c_aei.txt", "w");
    FILE* fd = fopen("c_d.txt", "w");
    // 打开赤经赤纬输出文件（单位：微角秒）
    FILE* fra_dec = fopen("ra_dec.txt", "w");
    FILE* fxyz_rv = fopen("c_xyz_rv.txt", "w");
    if (!fra_dec) { perror("ra_dec.txt"); return 1; }
    if (!fd) { perror("c_d.txt"); return 1; }
    if (!fxyz_rv) { perror("c_xyz_rv.txt"); return 1; }
    
    /* 从周期计算半长轴，转换角度单位 */
    double a0[MAX_BODIES], a[MAX_BODIES];
    double e[MAX_BODIES], ci[MAX_BODIES], w[MAX_BODIES], omg[MAX_BODIES], cm0[MAX_BODIES], cm[MAX_BODIES];
    for (i = 1; i < nbod; i++) {
        a0[i] = pow(pow(p0[i] / 365.2422, 2) * xmass[0], 1.0 / 3.0);
        a[i] = a0[i];
        e[i] = e0[i];
        ci[i] = ci0[i] * tpi / 360.0;  // 转换为弧度
        w[i] = w0[i] * tpi / 360.0;
        omg[i] = omg0[i] * tpi / 360.0;
        cm0[i] = cm00[i] * tpi / 360.0;
    }

    /* 输出初始参数 */
    printf("Initial a, e, mass (bodies 2..%d):\n", nbod);
    for (i = 1; i < nbod; i++) {
        printf("% .6f  % .6f  % .6f\n", a[i], e[i], xmass[i]);
    }

    /* 状态向量初始化 */
    double x[MAX_BODIES][6] = { {0} };  // 位置和速度
    double xp[6];                       // 临时存储单个天体的坐标
    double xrk[MAX_BODIES * 6];         // 积分器用的1D数组


    // 初始化直角坐标
    for (i = 1; i < nbod; i++) {
        ele2rat(a[i], e[i], ci[i], w[i], omg[i], cm0[i], i, xp);
        for (j = 0; j < ndim; j++) x[i][j] = xp[j];
    }

    // 质心校准
    double dx[6] = { 0 };
    for (j = 0; j < ndim; j++)
    {
        dx[j] = 0.0;
        for (i = 1; i < nbod; i++)
        {
            dx[j] += x[i][j] * xmass[i];
        }
    }
    double tmass = 0.0;
    for (i = 0; i < nbod; i++) tmass += xmass[i];
    for (j = 0; j < ndim; j++) x[0][j] = -dx[j] / tmass;

    // 调整坐标（相对质心）
    for (i = 1; i < nbod; i++)
    {
        for (j = 0; j < ndim; j++)
        {
            x[i][j] += x[0][j];
        }
    }

    // 转换为1D数组供积分器使用
    for (i = 0; i < nbod; i++)
    {
        for (j = 0; j < ndim; j++)
        {
            xrk[i * ndim + j] = x[i][j];
        }
    }

    // 积分控制参数
    double h = tpi / 365.25 / 100.0;    // 初始步长
    double err = 1e-16;                 // 误差 tolerance
    double t = 0.0;                     // 初始时间
    double tstop = (60 / 365.25) * tpi;       // 总积分时间
    //double tout0 = tstop / 300.0;       // 输出间隔（约100年）
    double tout0 = 0.02 / 365.25 * tpi; // 输出间隔（10天）
    double tout = tout0;                // 下次输出时间

    /*AU_IN_M：1 天文单位（AU）的米数（149597870700.0 m）
    SEC_PER_DAY：1 天的秒数（86400.0 s）
    SEC_PER_YEAR：1 年的秒数（365.25 * SEC_PER_DAY）
    AU_PER_YEAR_TO_MPS：将速度单位从 “AU / 年” 转换为 “米 / 秒” 的因子，计算方式为 AU_IN_M / SEC_PER_YEAR（ 1 AU / 年 = 1 AU 的米数 ÷ 1 年的秒数）*/

    const double AU_IN_M = 149597870700.0;            /* IAU 定义的 1 AU (m) */
    const double SEC_PER_DAY = 86400.0;
    const double SEC_PER_YEAR = 365.25 * SEC_PER_DAY;
    const double AU_PER_YEAR_TO_MPS = AU_IN_M / SEC_PER_YEAR; /* AU/year -> m/s */
    const double AU_PER_PC = 206264.80624709636; /* 1 pc = 206264.806... AU */

    double dist_au = observer_distance_pc * AU_PER_PC; /* 将 pc -> AU */
    double PI = tpi / 2.0;
    double rad2uas = (180.0 * 3600.0 * 1e6) / PI; /* 1 rad -> micro-arcsec */

    // 观测者在 +x 方向
    const double observer_x = dist_au;
    const double observer_y = 0.0;
    const double observer_z = 0.0;

    /* 主积分循环 */
    while (t < tstop) {
        // 更新坐标
        for (i = 0; i < nbod; i++)
        {
            for (j = 0; j < ndim; j++)
            {
                x[i][j] = xrk[i * ndim + j];
            }
        }

        // 转换为轨道要素
        for (i = 1; i < nbod; i++)
        {
            for (j = 0; j < ndim; j++)
            {
                xp[j] = x[i][j] - x[0][j];
            }
            rat2ele(xp, &a[i], &e[i], &ci[i], &w[i], &omg[i], &cm[i], i, 0);
        }

        // 输出数据
        if (t >= tout) {
            // 输出轨道要素
            fprintf(faei,
                "%20.6f"  // t/tpi
                "%20.6f"  // a[1]
                "%20.6f"  // e[1]
                "%20.6f"  // ci[1]
                "%20.6f"  // w[1]
                "%20.6f"  // omg[1]
                "%20.6f"  // cm[1]
                "%20.6f\n", // xmass[1]
                t / tpi,
                a[1],
                e[1],
                ci[1],
                w[1],
                omg[1],
                cm[1],
                xmass[1] / 3e-6
            );

            fflush(faei);

            // 输出距离和单位向量
            fprintf(fd, "%20.6f", t / tpi);
            for (i = 0; i < nbod; i++) {
                double d = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2]);
                fprintf(fd, "%20.6f", d);
                if (d != 0) {
                    fprintf(fd, "%20.6f%20.6f%20.6f", x[i][0] / d, x[i][1] / d, x[i][2] / d);
                }
                else {
                    fprintf(fd, "%20.0f%20.0f%20.0f", 0.0, 0.0, 0.0);
                }
            }
            fprintf(fd, "\n");
            fflush(fd);

            // --- 输出 xyz, vxvyvz (AU/year) 和 RV (m/s) 
            fprintf(fxyz_rv, "%20.6f", (t / tpi) * 365.25); // 第1列：时间

            for (i = 0; i < nbod; i++) {
                // 获取当前天体的位置和速度
                double x_pos = x[i][0], y_pos = x[i][1], z_pos = x[i][2];
                double vx_rad = x[i][3], vy_rad = x[i][4], vz_rad = x[i][5];

                // 把速度从 AU/rad -> AU/year（因为 t 中 2π 对应 1 年）
                double vx_au_per_year = vx_rad * tpi;
                double vy_au_per_year = vy_rad * tpi;
                double vz_au_per_year = vz_rad * tpi;

                /* 观测方向在 +x，视线单位向量约为 (1, 0, 0)。在遥远观测者近似下
                   视向速度约等于 x 方向速度分量（AU/year）。 */
                double rv_au_per_year = - vz_au_per_year;
                double rv_mps = rv_au_per_year * AU_PER_YEAR_TO_MPS;

                /* 输出：x,y,z (AU) ; vx,vy,vz (AU/year) ; RV (m/s) */
                fprintf(fxyz_rv,
                    "%20.6f%20.6f%20.6f"  // x, y, z (AU)
                    "%20.6f%20.6f%20.6f"  // vx, vy, vz (AU/year)
                    "%20.6f",             // RV (m/s)
                    x_pos, y_pos, z_pos,
                    vx_au_per_year, vy_au_per_year, vz_au_per_year,
                    rv_mps
                );
            }

            fprintf(fxyz_rv, "\n");
            fflush(fxyz_rv);

            double t_days = (t / tpi) * 365.25; /* t -> days */
            fprintf(fra_dec, "%.6f", t_days);

            /* 计算主星的小角近似 RA/Dec（相对于观测者，单位 µas） */
            double sx = x[0][0];   /* 主星 X (AU) */
            double sy = x[0][1];   /* 主星 Y (AU) */
            double sz = x[0][2];   /* 主星 Z (AU) */

            double star_ra_rad = sx / dist_au;
            double star_dec_rad = sy / dist_au;

            double star_ra_uas = star_ra_rad * rad2uas;
            double star_dec_uas = star_dec_rad * rad2uas;

            /* 仅写入两列：star_ra_uas, star_dec_uas */
            fprintf(fra_dec, "%20.6f  %20.6f\n", star_ra_uas, star_dec_uas);
            fflush(fra_dec);

            tout += tout0; 
        }
        
        start_integrate = clock();
        /* 自适应RK步长积分 */
        rkf78_no_t(&t, &h, err, xrk, nbod * ndim);

        /* 稳定性检查 */
        for (i = 1; i < nbod; i++) {
            if (fabs(a[i] - a0[i]) > 0.1 * a0[i]) { break; }
        }
        end_integrate = clock();
        integrate_time_used += ((double)(end_integrate - start_integrate)) / CLOCKS_PER_SEC;
    }

    printf("Total integration timem : %f s\n", integrate_time_used);
    
    // 计算并输出运行时间
    /*end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Program run time: %f s\n", cpu_time_used);*/

    // 关闭文件
    fclose(faei);
    fclose(fd);
    fclose(fra_dec);
    fclose(fxyz_rv);
    return 0;
}

