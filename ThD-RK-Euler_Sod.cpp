#include <math.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <string>
#include <fstream>
#include <iterator>
#include <valarray>
using namespace std;
// example:Euler 3-4-L-W   特征分裂
// g++ -O2 -std=c++11 -Wall -Wextra "-Wl,--stack=1073741824" SDMMs-Euler_Exm.cpp -o ref.exe
class cfda
{

private:
    int w;

public:
    double PI = 3.141592653589793;
    // double PI = 3.1415926535897932385;

    int n;                                                         // 网格结点数（有限差分）
    int nb;                                                        // 边界网格结点数（有限差分）
    double t, tt, CFL;                                             // 计算总时 间t,实际计算时间tt;
    int nt, kt, ij, ijj, onet_out = 2;                             // 计算时间步长;
    double nx_begin, nx_end;                                       // 网格长度;
    int n1, ii1, m;                                                // 总的网格结点数（有限差分）
    double *Ltt_n0, *x, *Ltt_n1, *Lttx_n1;                         // 声明变长数组
    double *L_n0, *L_n1, *L_n2;                                    // f=L
    double *Lt_n0, *Lt_n1, *Lt_n2;                                 // L_t
    double *Ltx_n0, *Ltx_n1, *Ltx_n2;                              // L_t
    double *u_n0, *u_n1, *u_n2, *u_n3, *u_n4, *u_n5, *u_n6, *u_n7; // 推进步中间时刻
    double *f_n0, *f_n1, *f_n2, *f_n3, *f_n4, *f_n5, *f_n6, *f_n7; // 推进步中间时刻
    double *u_nn, *f_nn, *Lttx_n0;                                 // 下一时刻（n+1)
    double *TL_n0, *TLt_n0, *TL_n1, *TLt_n1, *TL_n2, *TLt_n2, *Tu1, *Tu2;
    double *flux, *fluxz, *fluxf;  // j+1/2处的通量
    double *df, *u_excat, *uc_max; // du[513 + 6 * 2],df[513 + 6 * 2],
    ////-----------------------------------------
    double dx;
    double dt;
    double *dd; // ss
    double A1, A2, A3, A4, A5, B0, B1, B2, C0, C1, C2, D0, D1, D2, D3, maxuu, TVmax;
    double a21, a31, a32, aa12, aa21, aa31, aa32, aaa31, aaa32, w11, w22, v0, v1, v2, v3, vv0, vv1, vv2, vv3, ww1, ww2, b0, b1, b2, b3, bb1, bb2, bb3, bbb1, bbb2, bbb3;
    double put_obj[10000]; // //onet_out = 10;
    double put_one_n[10000];
    // Euler参数——————————————————————————————————————————————————
    double *u, *c, *p, *a;   // 原始变量  u速度，c声速，p压强，a密度
    double *fz, *ff;         // 守恒通量
    double fza[61], ffa[61]; // 特征变换守恒通量
    double y = 1.4;          // 热力系数
                             ////-----------------------------------------
    void pc_to_uu(double *uuu)
    {
        for (int i = 0; i < n1; i++) // Sod激波管
        {
            uuu[i] = a[i], uuu[n1 + i] = a[i] * u[i], uuu[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
        }
    }

    void pc_to_ff(double *fff)
    {
        for (int i = 0; i < n1; i++) // Sod激波管
        {
            fff[i] = a[i] * u[i], fff[n1 + i] = a[i] * u[i] * u[i] + p[i];
            fff[2 * n1 + i] = (p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0 + p[i]) * u[i];
        }
    }

    void uu_to_ff(double *uuu, double *fff)
    {
        for (int i = 0; i < n1; i++) // Sod激波管
        {
            fff[i] = uuu[i] * u[i], fff[n1 + i] = uuu[i] * u[i] + p[i];
            fff[2 * n1 + i] = (uuu[i] + p[i]) * u[i];
        }
    }

    void uu_to_pc(double *uuu)
    {
        for (int i = 0; i < n1; i++) // Sod激波管
        {
            a[i] = uuu[i], u[i] = uuu[n1 + i] / a[i];
            p[i] = (y - 1.0) * (uuu[n1 * 2 + i] - 0.5 * uuu[n1 + i] * u[i]);
            c[i] = sqrt(y * p[i] / a[i]);
        }
    }

    void intc()
    {
        double p1, a1, u1, c1, p2, a2, u2, c2, ssss;

        dx = (nx_end - nx_begin) / (1.0 * (n - 1)); // Sod激波管   mesh=200  [0-1]  t=0.2    dt = CFL * dx / 2.;   //Sod
        a1 = 1.0, u1 = 0.00, p1 = 1.0;
        a2 = 0.125, u2 = 0.00, p2 = 0.1;
        c1 = sqrt(y * p1 / a1), c2 = sqrt(y * p2 / a2);
        ssss = sqrt(2.0);
        for (int i = 0; i < n1; i++)
        {
            x[i] = nx_begin + (i - nb) * dx;
            if (x[i] < (nx_end - nx_begin) / 2.0)
            {
                u[i] = u1, c[i] = c1, p[i] = p1, a[i] = a1;
                u_n0[i] = a[i], u_n0[n1 + i] = a[i] * u[i], u_n0[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
            }
            else
            {
                u[i] = u2, c[i] = c2, p[i] = p2, a[i] = a2;
                u_n0[i] = a[i], u_n0[n1 + i] = a[i] * u[i], u_n0[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
            }
        }

        // dx = (nx_end - nx_begin) / (1.0 * (n - 1)); // 真空  mesh=200  [0-1]  t=0.15    dt = CFL * dx / 2.;
        // u1 = -2.0, a1 = 1.0, p1 = 0.4;
        // u2 = 2.0, a2 = 1.0, p2 = 0.4;
        // //         u1 = -3.0, a1 = 1.0, p1 = 1.0;
        // // u2 = 3.0, a2 = 1.0, p2 = 1.0;
        // c1 = sqrt(y * p1 / a1), c2 = sqrt(y * p2 / a2);
        // ssss = sqrt(2.0);
        // for (int i = 0; i < n1; i++)
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     if (x[i] < (nx_end - nx_begin) / 2.0)
        //     {
        //         u[i] = u1, c[i] = c1, p[i] = p1, a[i] = a1;
        //         u_n0[i] = a[i], u_n0[n1 + i] = a[i] * u[i], u_n0[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
        //     }
        //     else
        //     {
        //         u[i] = u2, c[i] = c2, p[i] = p2, a[i] = a2;
        //         u_n0[i] = a[i], u_n0[n1 + i] = a[i] * u[i], u_n0[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
        //     }
        // }

        // dx = (nx_end - nx_begin) / (1.0 * (n - 1)); // Shu-Osher   mesh=500  [0-10]  t=1.8     // dt = CFL * dx / 5.0;  //Shu-Osher
        // for (int i = 0; i < n1; i++)
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     if (x[i] <= 1.0)
        //     {
        //         a1 = 3.857134, u1 = 2.629369, p1 = 10.33333, c1 = sqrt(y * p1 / a1);
        //         u[i] = u1, c[i] = c1, p[i] = p1, a[i] = a1;
        //         u_n0[i] = a[i], u_n0[n1 + i] = a[i] * u[i], u_n0[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
        //     }
        //     else
        //     {
        //         a2 = 1. + 0.2 * sin(x[i] * 5.), u2 = 0.00, p2 = 1.0, c2 = sqrt(y * p2 / a2);
        //         u[i] = u2, c[i] = c2, p[i] = p2, a[i] = a2;
        //         u_n0[i] = a[i], u_n0[n1 + i] = a[i] * u[i], u_n0[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
        //     }
        // }

        // dx = (nx_end - nx_begin) / (1.0 * (n - 1));
        // for (int i = 0; i < n1; i++) //  blast wave problem   mesh=500  [0-10]  t=0.38
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     if (x[i] < 1.0)
        //     {
        //         a1 = 1.0, u1 = 0.0, p1 = 1000., c1 = sqrt(y * p1 / a1);
        //         u[i] = u1, c[i] = c1, p[i] = p1, a[i] = a1;
        //         u_n0[i] = a[i], u_n0[n1 + i] = a[i] * u[i], u_n0[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
        //     }
        //     else if (x[i] < 9.0)
        //     {
        //         a2 = 1., u2 = 0.0, p2 = 0.01, c2 = sqrt(y * p2 / a2);
        //         u[i] = u2, c[i] = c2, p[i] = p2, a[i] = a2;
        //         u_n0[i] = a[i], u_n0[n1 + i] = a[i] * u[i], u_n0[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
        //     }
        //     else
        //     {
        //         a2 = 1., u2 = 0.0, p2 = 100.0, c2 = sqrt(y * p2 / a2);
        //         u[i] = u2, c[i] = c2, p[i] = p2, a[i] = a2;
        //         u_n0[i] = a[i], u_n0[n1 + i] = a[i] * u[i], u_n0[2 * n1 + i] = p[i] / (y - 1.0) + a[i] * u[i] * u[i] / 2.0;
        //     }
        // }

        Store_obj(put_obj, u, a, p);
        Write_obj(0);
    }

    double intc_fun(double xxf)
    {
        double yyf;
        // yyf = (sin(2 * PI * xxf) / 2.0 + 0.5) * 1.0; //  / PIbergure
        // yyf = sin(PI * xxf) / 2.0 + 0.25;//bergure
        yyf = sin(PI * xxf); // line
        return yyf;
    }

    double dmax1(double a, double b, double c)
    {
        if (b > a)
        {
            a = b;
        }
        if (c > a)
        {
            a = c;
        }
        return a;
    }
    // //-----------------------------------------
    void border()
    {
        borderfun(u_n0), borderfun(u_n1), borderfun(u_n2), borderfun(u_nn);
        borderfun(f_n0), borderfun(f_n1), borderfun(f_n2), borderfun(f_nn);
        borderfun(L_n0), borderfun(L_n1), borderfun(L_n2), borderfun(Lt_n0);
        borderfun(Lt_n1), borderfun(Lt_n2), borderfun(Ltt_n0);
        borderfun(Ltx_n0), borderfun(Lttx_n0), borderfun(u_n3);
        borderfun(Ltx_n1), borderfun(Ltx_n2), borderfun(f_n3);
    }

    void borderfun(double *ffbc)
    {

        for (int i = 0; i < nb; i++)
        {

            // *(ffbc + i) = *(ffbc + n1 - 2 * nb + i - 1);  //周期边条
            // *(ffbc + n1 - nb + i) = *(ffbc + nb + i + 1); //周期边条

            // *(ffbc + i) = 0;
            //*(ffbc + n1 - nb + i) = 0;
            // f_n0[i] = f_n0[nb];                   //恒定边条
            // f_n0[n1 - 1 - i] = f_n0[n1 - 1 - nb]; //恒定边条

            *(ffbc + i) = *(ffbc + nb);                   // 恒定边条
            *(ffbc + n1 - 1 - i) = *(ffbc + n1 - 1 - nb); // 恒定边条

            // *(ffbc + i) = 0.0;                   //恒定边条
            // *(ffbc + n1 - 1 - i) = 0.0; //恒定边条
            // ffbc++;
        }
    }

    void borderfun3(double *ffbc)
    {

        for (int i = 0; i < nb; i++)
        {

            // *(ffbc + i) = *(ffbc + n1 - 2 * nb + i - 1);  //周期边条
            // *(ffbc + n1 - nb + i) = *(ffbc + nb + i + 1); //周期边条

            // *(ffbc + i) = 0;
            //*(ffbc + n1 - nb + i) = 0;
            // f_n0[i] = f_n0[nb];                   //恒定边条
            // f_n0[n1 - 1 - i] = f_n0[n1 - 1 - nb]; //恒定边条

            *(ffbc + i) = *(ffbc + nb);                                       // 恒定边条
            *(ffbc + n1 - 1 - i) = *(ffbc + n1 - 1 - nb);                     // 恒定边条
            *(ffbc + i + n1) = *(ffbc + nb + n1);                             // 恒定边条
            *(ffbc + n1 - 1 - i + n1) = *(ffbc + n1 - 1 - nb + n1);           // 恒定边条
            *(ffbc + i + n1 + n1) = *(ffbc + nb + n1 + n1);                   // 恒定边条
            *(ffbc + n1 - 1 - i + n1 + n1) = *(ffbc + n1 - 1 - nb + n1 + n1); // 恒定边条
        }
    }
    // //-----------------------------------------
    void carr(double *p1, double *p2, int i)
    {
        while (i-- > 0)
        {
            *p1++ = *p2++;
        }
    }
    // //-----------------------------------------
    void Splitting(double *uuu) //_global
    {
        double lmax = 0.0, lmax0;
        for (int i = 0; i < n1; i++) // comput L(n0)
        {
            lmax0 = dmax1(abs(u[i]), abs(u[i] - c[i]), abs(u[i] + c[i])); // 当地最大特征值
            if (lmax0 > lmax)
            {
                lmax = lmax0;
            }
        }
        lmax = lmax * 1.0;
        for (int i = 0; i < n1; i++) // comput L(n0)
        {
            fz[i] = 0.5 * (a[i] * u[i] + lmax * uuu[i]); // 正通量
            fz[n1 + i] = 0.5 * (a[i] * u[i] * u[i] + p[i] + lmax * uuu[n1 + i]);
            fz[n1 + n1 + i] = 0.5 * ((uuu[n1 + n1 + i] + p[i]) * u[i] + lmax * uuu[n1 + n1 + i]);
            ff[i] = 0.5 * (a[i] * u[i] - lmax * uuu[i]); // 负通量
            ff[n1 + i] = 0.5 * (a[i] * u[i] * u[i] + p[i] - lmax * uuu[n1 + i]);
            ff[n1 + n1 + i] = 0.5 * ((uuu[n1 + n1 + i] + p[i]) * u[i] - lmax * uuu[n1 + n1 + i]);
        }

        borderfun3(ff);
        borderfun3(fz);
    }
    void Splitting_local(double *uuu) //
    {
        double lmax = 0.0, lmax0;
        for (int i = 0; i < n1; i++) // comput L(n0)
        {
            lmax = dmax1(abs(u[i]), abs(u[i] - c[i]), abs(u[i] + c[i])); // 当地最大特征值
            fz[i] = 0.5 * (a[i] * u[i] + lmax * uuu[i]);                 // 正通量
            fz[n1 + i] = 0.5 * (a[i] * u[i] * u[i] + p[i] + lmax * uuu[n1 + i]);
            fz[n1 + n1 + i] = 0.5 * ((uuu[n1 + n1 + i] + p[i]) * u[i] + lmax * uuu[n1 + n1 + i]);
            ff[i] = 0.5 * (a[i] * u[i] - lmax * uuu[i]); // 负通量
            ff[n1 + i] = 0.5 * (a[i] * u[i] * u[i] + p[i] - lmax * uuu[n1 + i]);
            ff[n1 + n1 + i] = 0.5 * ((uuu[n1 + n1 + i] + p[i]) * u[i] - lmax * uuu[n1 + n1 + i]);
        }

        borderfun3(ff);
        borderfun3(fz);
    }

    void Splitting_flux_char(double *uuu)
    {
        double gamma = 1.4, fla[3];
        double a1, u1, H1, a2, u2, H2, r1, r2, r0, u0, H0, c0, t1, s11, s12, s13, s21, s22, s23, s31, s32, s33;
        double t2, D11, D12, D13, D21, D22, D23, D31, D32, D33;
        int k1, mm;
        for (int j = nb - 1; j < n1 - nb; j++) // comput L(n0)
        {
            a1 = a[j];
            u1 = u[j];
            H1 = gamma / (gamma - 1.) * p[j] / a[j] + u[j] * u[j] * 0.5; // ! 总焓
            a2 = a[j + 1];
            u2 = u[j + 1];
            H2 = gamma / (gamma - 1.) * p[j + 1] / a[j + 1] + u[j + 1] * u[j + 1] * 0.5;

            // !        Roe 平均
            r1 = sqrt(a1);
            r2 = sqrt(a2);
            r0 = r1 + r2;
            u0 = (r1 * u1 + r2 * u2) / r0;
            H0 = (r1 * H1 + r2 * H2) / r0;
            c0 = sqrt((gamma - 1.) * (H0 - u0 * u0 * 0.5)); // ! j+1/2点处的声速

            // !  特征矩阵S  (A的右特征矩阵)
            t1 = (gamma - 1.) / c0;
            s11 = u0 * u0 / 2. - c0 * c0 / (gamma - 1.), s12 = -u0, s13 = 1.;
            s21 = -u0 - t1 * u0 * u0 / 2., s22 = 1. + t1 * u0, s23 = -t1;
            s31 = -u0 + t1 * u0 * u0 / 2., s32 = 1. - t1 * u0, s33 = t1;

            // ! 特征矩阵S^-1  (S的逆矩阵，A的左特征矩阵)
            t2 = (gamma - 1.) / (c0 * c0);
            D11 = -t2, D12 = -0.5 / c0, D13 = 0.5 / c0;
            D21 = -t2 * u0, D22 = -(u0 - c0) / (2. * c0), D23 = (u0 + c0) / (2. * c0);
            D31 = -t2 * u0 * u0 / 2., D32 = -(u0 * u0 / 2. + 1. / t2 - u0 * c0) / (2. * c0);
            D33 = (u0 * u0 / 2. + 1. / t2 + u0 * c0) / (2. * c0);

            for (int j0 = -4; j0 < 5; j0++)
            {
                fluxz[j + j0] = s11 * fz[j + j0] + s12 * fz[n1 + j + j0] + s13 * fz[n1 + n1 + j + j0];
                fluxz[n1 + j + j0] = s21 * fz[j + j0] + s22 * fz[n1 + j + j0] + s23 * fz[n1 + n1 + j + j0];
                fluxz[n1 + n1 + j + j0] = s31 * fz[j + j0] + s32 * fz[n1 + j + j0] + s33 * fz[n1 + n1 + j + j0];

                fluxf[j + j0] = s11 * ff[j + j0] + s12 * ff[n1 + j + j0] + s13 * ff[n1 + n1 + j + j0];
                fluxf[n1 + j + j0] = s21 * ff[j + j0] + s22 * ff[n1 + j + j0] + s23 * ff[n1 + n1 + j + j0];
                fluxf[n1 + n1 + j + j0] = s31 * ff[j + j0] + s32 * ff[n1 + j + j0] + s33 * ff[n1 + n1 + j + j0];
            }
            for (int k = 0; k < 3; k++)
            {
                mm = j + k * n1;
                fla[k] = df_dx_char(fluxz, fluxf, mm);
            }

            flux[j] = D11 * fla[0] + D12 * fla[1] + D13 * fla[2];
            flux[n1 + j] = D21 * fla[0] + D22 * fla[1] + D23 * fla[2];
            flux[n1 + n1 + j] = D31 * fla[0] + D32 * fla[1] + D33 * fla[2];

            // lmax = dmax1(abs(u[i]), abs(u[i] - c[i]), abs(u[i] + c[i])); // 当地最大特征值
            // fz[i] = 0.5 * (a[i] * u[i] + lmax * uuu[i]);                 //正通量
            // fz[n1 + i] = 0.5 * (a[i] * u[i] * u[i] + p[i] + lmax * uuu[n1 + i]);
            // fz[n1 + n1 + i] = 0.5 * ((uuu[n1 + n1 + i] + p[i]) * u[i] + lmax * uuu[n1 + n1 + i]);
            // ff[i] = 0.5 * (a[i] * u[i] - lmax * uuu[i]); //负通量
            // ff[n1 + i] = 0.5 * (a[i] * u[i] * u[i] + p[i] - lmax * uuu[n1 + i]);
            // ff[n1 + n1 + i] = 0.5 * ((uuu[n1 + n1 + i] + p[i]) * u[i] - lmax * uuu[n1 + n1 + i]);
        }
        borderfun3(flux);
    }

    //------------------------------------------------
    void Store_obj(double *ss, double *us, double *as, double *ps)
    {
        for (int ik = 0; ik < n1; ik++) //
        {
            ss[ik] = as[ik], ss[n1 + ik] = us[ik], ss[n1 * 2 + ik] = ps[ik];
        }
    }

    void Write_obj(int j)
    {
        if (j < -0.9)
        {
            // string title = "result_t_burgers_excat"; //+ to_string(t);
            string title = to_string(n);
            string Title = title + ".plt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=x,a,u,p " << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (int j = 0; j < 3; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_obj[j * n1 + i];
                }

                ofs << endl;
            }
            ofs.close();
        }

        if (j < 0.5 && j > -0.8)
        {
            string Title = "int.plt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=x,a,u,p " << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (int j = 0; j < 3; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_obj[j * n1 + i];
                }

                ofs << endl;
            }
            ofs.close();
        }

        if (j > 1.5)
        {
            string title = "result_one_t_burgers_every-" + to_string(j);
            string Title = title + ".plt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << "Title=" << title << endl;
            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (j = 0; j < onet_out; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_one_n[i];
                }

                ofs << endl;
            }
            ofs.close();
        }
    }
    ////-----------------------------------------
    void comput_LWqt0(double *fqtt, double *qt, double *uuu)
    {
        int j, k;
        double fq[4][4];
        for (int i = 0; i < n1; i++) // comput L(n0)
        {
            fq[1][1] = 0.0, fq[1][2] = 1.0, fq[1][3] = 0.0;
            fq[2][1] = (y - 3.0) * uuu[n1 + i] * uuu[n1 + i] / 2.0 / uuu[i] / uuu[i];
            fq[2][2] = (3.0 - y) * uuu[n1 + i] / uuu[i];
            fq[2][3] = y - 1.0;
            fq[3][1] = ((y - 1.0) * uuu[n1 + i] * uuu[n1 + i] * uuu[n1 + i] - y * uuu[i] * uuu[n1 + i] * uuu[n1 + n1 + i]);
            fq[3][1] = fq[3][1] / uuu[i] / uuu[i] / uuu[i];
            fq[3][2] = (3 * (1.0 - y) * uuu[n1 + i] * uuu[n1 + i] + 2 * y * uuu[i] * uuu[n1 + n1 + i]) / uuu[i] / uuu[i] / 2.0;
            fq[3][3] = y * uuu[n1 + i] / uuu[i];
            fqtt[i] = -1.0 * (fq[1][1] * qt[i] + fq[1][2] * qt[n1 + i] + fq[1][3] * qt[n1 + n1 + i]);
            fqtt[n1 + i] = -1.0 * (fq[2][1] * qt[i] + fq[2][2] * qt[n1 + i] + fq[2][3] * qt[n1 + n1 + i]);
            fqtt[n1 + n1 + i] = -1.0 * (fq[3][1] * qt[i] + fq[3][2] * qt[n1 + i] + fq[3][3] * qt[n1 + n1 + i]);
        }
        borderfun3(fqtt);
    }
    //  comput_LWqtt0(Lttx_n0, Lt_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
    // void comput_LWqtt0(double *qttt, double *qtt, double *qt, double *uuu)
    // {
    //     int j, k;
    //     double q1, q2, q3, qt1, qt2, qt3, qtt1, qtt2, qtt3;
    //     for (int i = 0; i < n1; i++) // comput L(n0)
    //     {
    //         q1 = uuu[i], q2 = uuu[n1 + i], q3 = uuu[n1 + n1 + i];
    //         qt1 = qt[i], qt2 = qt[n1 + i], qt3 = qt[n1 + n1 + i];
    //         qtt1 = qtt[i], qtt2 = qtt[n1 + i], qtt3 = qtt[n1 + n1 + i];
    //         qttt[i] = -qtt2;
    //         qttt[n1 + i] = (2. * pow(q1, 2) * (pow(qt2, 2) * (y - 3.) - q1 * qtt3 * (y - 1.)) +
    //                         pow(q2, 2) * (2. * pow(qt1, 2) - q1 * qtt1) * (y - 3.) +
    //                         2. * q1 * q2 * (-2. * qt1 * qt2 + q1 * qtt2) * (y - 3.)) /
    //                        (2. * pow(q1, 3));
    //         qttt[n1 + n1 + i] = (2. * pow(q2, 3) * (3. * pow(qt1, 2) - q1 * qtt1) * (y - 1) +
    //                              3 * q1 * pow(q2, 2) * (-4. * qt1 * qt2 + q1 * qtt2) * (y - 1) -
    //                              2. * pow(q1, 2) * (-2. * q3 * qt1 * qt2 + 2. * q1 * qt2 * qt3 + q1 * q3 * qtt2) * y +
    //                              2. * q1 * q2 * (3. * q1 * pow(qt2, 2) * (y - 1) - 2. * q3 * pow(qt1, 2) * y + q1 * (2. * qt1 * qt3 + q3 * qtt1 - q1 * qtt3) * y)) /
    //                             (2. * pow(q1, 4));
    //     }
    //     borderfun3(qttt);
    // }

    void comput_LWqtt0(double *qttt, double *qtt, double *qt, double *uuu)
    {
        int j, k;
        double q1, q2, q3, qt1, qt2, qt3, qtt1, qtt2, qtt3, mq1e2, mq2e2, mqt1e2, mqt2e2, qt1qt2, xxx1, y3, y1;
        for (int i = 0; i < n1; i++) // comput L(n0)
        {
            q1 = uuu[i], q2 = uuu[n1 + i], q3 = uuu[n1 + n1 + i];
            qt1 = qt[i], qt2 = qt[n1 + i], qt3 = qt[n1 + n1 + i];
            qtt1 = qtt[i], qtt2 = qtt[n1 + i], qtt3 = qtt[n1 + n1 + i];
            qttt[i] = -qtt2;
            mq1e2 = 2. * q1 * q1, mq2e2 = q2 * q2, mqt1e2 = 2. * qt1 * qt1, mqt2e2 = qt2 * qt2;
            qt1qt2 = -2. * qt1 * qt2, xxx1 = qt1qt2 + q1 * qtt2, y3 = y - 3., y1 = y - 1.;
            qttt[n1 + i] = (mq1e2 * (mqt2e2 * y3 - q1 * qtt3 * y1) + (mq2e2 * (mqt1e2 - q1 * qtt1) + 2. * q1 * q2 * xxx1) * y3) /
                           (mq1e2 * q1);
            qttt[n1 + n1 + i] = ((2. * q2 * (1.5 * mqt1e2 - q1 * qtt1) + 3 * q1 * (xxx1 + qt1qt2)) * y1 * mq2e2 -
                                 mq1e2 * (q3 * qt1qt2 + 2. * q1 * qt2 * qt3 + q1 * q3 * qtt2) * y +
                                 2. * q1 * q2 * (3. * q1 * mqt2e2 * y1 + (-q3 * mqt1e2 + q1 * (2. * qt1 * qt3 + q3 * qtt1 - q1 * qtt3)) * y)) /
                                (0.5 * mq1e2 * mq1e2);
        }
        borderfun3(qttt);
    }

    void compt_ThD13()
    {
        int j = 0, i, k, ii;

        uu_to_pc(u_n0);
        Splitting(u_n0);
        // Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n0, Lt_n0, L_n0, u_n0); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n0);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n0[ii] = dfdxcent(Lttx_n0, ii);
                u_nn[ii] = u_n0[ii] + dt * L_n0[ii] + dt * dt / 2. * Lt_n0[ii] + dt * dt * dt / 6. * Ltt_n0[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }
    void compt_ThD13_Char()
    {
        int j = 0, i, k, ii;

        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n0, Lt_n0, L_n0, u_n0); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n0);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n0[ii] = dfdxcent(Lttx_n0, ii);
                u_nn[ii] = u_n0[ii] + dt * L_n0[ii] + dt * dt / 2. * Lt_n0[ii] + dt * dt * dt / 6. * Ltt_n0[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void compt_ThD23_Char()
    {
        int j = 0, i, k, ii;

        a21 = 0.500;
        b2 = 0.5;
        bb2 = 0.125;
        bbb1 = 0.020833333333333;
        bbb2 = 1. / 6. * (1. - 3. * a21 * (a21 * b2 + 2. * bb2) - 6. * bbb1);
        b1 = 1. - b2;
        bb1 = 1. / 2. - a21 * b2 - bb2;

        // a21 = 0.614823, b2 = 0.40239, bb2 = 0.082793, b1 = 1 - b2;
        // bb1 = 1. / 2. - a21 * b2 - bb2;
        // bbb1 = 1. / 24. * (4. - 1. / a21 - 8. * a21 * a21 * b2 - 12. * a21 * bb2);
        // bbb2 = 1. / (24. * a21) - 1. / 6. * a21 * (a21 * b2 + 3. * bb2);
        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n0, Lt_n0, L_n0, u_n0); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n0);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n0[ii] = dfdxcent(Lttx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 * dt * dt / 2. * Lt_n0[ii] + a21 * a21 * a21 * dt * dt * dt / 6. * Ltt_n0[ii];
            }
            j++;
        }
        borderfun3(u_n1);

        uu_to_pc(u_n1);
        Splitting(u_n1);
        Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n1, Lt_n1, L_n1, u_n1); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n1);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n1[ii] = dfdxcent(Lttx_n1, ii);
                u_nn[ii] = u_n0[ii] + b1 * dt * L_n0[ii] + bb1 * dt * dt * Lt_n0[ii] + bbb1 * dt * dt * dt * Ltt_n0[ii] +
                           b2 * dt * L_n1[ii] + bb2 * dt * dt * Lt_n1[ii] + bbb2 * dt * dt * dt * Ltt_n1[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void compt_ThD24()
    {
        int j = 0, i, k, ii;
        a21 = 0.614823, b2 = 0.40239, bb2 = 0.082793, b1 = 1 - b2;
        bb1 = 1. / 2. - a21 * b2 - bb2;
        bbb1 = 1. / 24. * (4. - 1. / a21 - 8. * a21 * a21 * b2 - 12. * a21 * bb2);
        bbb2 = 1. / (24. * a21) - 1. / 6. * a21 * (a21 * b2 + 3. * bb2);
        uu_to_pc(u_n0);
        Splitting(u_n0);
        // Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n0, Lt_n0, L_n0, u_n0); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n0);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n0[ii] = dfdxcent(Lttx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 * dt * dt / 2. * Lt_n0[ii] + a21 * a21 * a21 * dt * dt * dt / 6. * Ltt_n0[ii];
            }
            j++;
        }
        borderfun3(u_n1);

        uu_to_pc(u_n1);
        Splitting(u_n1);
        // Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n1, Lt_n1, L_n1, u_n1); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n1);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n1[ii] = dfdxcent(Lttx_n1, ii);
                u_nn[ii] = u_n0[ii] + b1 * dt * L_n0[ii] + bb1 * dt * dt * Lt_n0[ii] + bbb1 * dt * dt * dt * Ltt_n0[ii] +
                           b2 * dt * L_n1[ii] + bb2 * dt * dt * Lt_n1[ii] + bbb2 * dt * dt * dt * Ltt_n1[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }
    void compt_ThD24_Char()
    {
        int j = 0, i, k, ii;
        a21 = 0.614823, b2 = 0.40239, bb2 = 0.082793, b1 = 1 - b2;
        bb1 = 1. / 2. - a21 * b2 - bb2;
        bbb1 = 1. / 24. * (4. - 1. / a21 - 8. * a21 * a21 * b2 - 12. * a21 * bb2);
        bbb2 = 1. / (24. * a21) - 1. / 6. * a21 * (a21 * b2 + 3. * bb2);
        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n0, Lt_n0, L_n0, u_n0); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n0);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n0[ii] = dfdxcent(Lttx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 * dt * dt / 2. * Lt_n0[ii] + a21 * a21 * a21 * dt * dt * dt / 6. * Ltt_n0[ii];
            }
            j++;
        }
        borderfun3(u_n1);

        uu_to_pc(u_n1);
        Splitting(u_n1);
        Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n1, Lt_n1, L_n1, u_n1); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n1);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n1[ii] = dfdxcent(Lttx_n1, ii);
                u_nn[ii] = u_n0[ii] + b1 * dt * L_n0[ii] + bb1 * dt * dt * Lt_n0[ii] + bbb1 * dt * dt * dt * Ltt_n0[ii] +
                           b2 * dt * L_n1[ii] + bb2 * dt * dt * Lt_n1[ii] + bbb2 * dt * dt * dt * Ltt_n1[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void compt_ThD25()
    {
        int j = 0, i, k, ii;
        a21 = 0.6;
        b1 = 1., b2 = 0., bb1 = (2. - 5. * a21 + 10. * pow(a21, 3)) / (20. * pow(a21, 3));
        bb2 = (-2. + 5. * a21) / (20. * pow(a21, 3));
        bbb1 = (3. + 10. * (-1. + a21) * a21) / (60. * pow(a21, 2));
        bbb2 = (3. - 5. * a21) / (60. * pow(a21, 2));
        uu_to_pc(u_n0);
        Splitting(u_n0);
        // Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n0, Lt_n0, L_n0, u_n0); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n0);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n0[ii] = dfdxcent(Lttx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 * dt * dt / 2. * Lt_n0[ii] + a21 * a21 * a21 * dt * dt * dt / 6. * Ltt_n0[ii];
            }
            j++;
        }
        borderfun3(u_n1);

        uu_to_pc(u_n1);
        Splitting(u_n1);
        // Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n1, Lt_n1, L_n1, u_n1); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n1);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n1[ii] = dfdxcent(Lttx_n1, ii);
                u_nn[ii] = u_n0[ii] + b1 * dt * L_n0[ii] + bb1 * dt * dt * Lt_n0[ii] + bbb1 * dt * dt * dt * Ltt_n0[ii] +
                           b2 * dt * L_n1[ii] + bb2 * dt * dt * Lt_n1[ii] + bbb2 * dt * dt * dt * Ltt_n1[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }
    void compt_ThD25_Char()
    {
        int j = 0, i, k, ii;
        a21 = 0.6;
        b1 = 1., b2 = 0., bb1 = (2. - 5. * a21 + 10. * pow(a21, 3)) / (20. * pow(a21, 3));
        bb2 = (-2. + 5. * a21) / (20. * pow(a21, 3));
        bbb1 = (3. + 10. * (-1. + a21) * a21) / (60. * pow(a21, 2));
        bbb2 = (3. - 5. * a21) / (60. * pow(a21, 2));
        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n0, Lt_n0, L_n0, u_n0); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n0);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n0[ii] = dfdxcent(Lttx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 * dt * dt / 2. * Lt_n0[ii] + a21 * a21 * a21 * dt * dt * dt / 6. * Ltt_n0[ii];
            }
            j++;
        }
        borderfun3(u_n1);

        uu_to_pc(u_n1);
        Splitting(u_n1);
        Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n1, Lt_n1, L_n1, u_n1); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n1);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n1[ii] = dfdxcent(Lttx_n1, ii);
                u_nn[ii] = u_n0[ii] + b1 * dt * L_n0[ii] + bb1 * dt * dt * Lt_n0[ii] + bbb1 * dt * dt * dt * Ltt_n0[ii] +
                           b2 * dt * L_n1[ii] + bb2 * dt * dt * Lt_n1[ii] + bbb2 * dt * dt * dt * Ltt_n1[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void compt_ThD36_Char()
    {
        int j = 0, i, k, ii;
        a21 = 0.869270422229783, b1 = 1., b2 = 0., b3 = 0.;
        aaa32 = 0., a32 = 0., bbb2 = 0., bbb3 = 0;
        bb1 = (-2. + 5. * pow(a21, 2.) * (16. + a21 * (-42. + 29. * a21))) / (60. * pow(2. - 3. * a21, 2) * pow(a21, 2.));
        bb2 = 1. / (60. * pow(a21, 2.) * (2. + a21 * (-6. + 5. * a21)));
        bb3 = pow(3. - 5. * a21, 4.) / (60. * pow(2. - 3. * a21, 2.) * (2. + a21 * (-6. + 5. * a21)));
        bbb1 = (1. + 5. * (-1. + a21) * a21) / (60. * a21 * (-2. + 3. * a21));
        a31 = (2. - 3. * a21) / (3. - 5. * a21);
        aa31 = (pow(2. - 3. * a21, 2.) * (-2. + a21 * (6. + a21 * (22. + 15. * a21 * (-6. + 5. * a21))))) / (6. * pow(3. - 5. * a21, 4.) * pow(a21, 2.));
        aa32 = (pow(2. - 3. * a21, 2.) * (2. + a21 * (-6. + 5. * a21))) / (6. * pow(3. - 5. * a21, 4.) * pow(a21, 2.));
        aaa31 = (pow(2. - 3. * a21, 2.) * (-2. + 3. * a21 * (4 + a21 * (-8 + 5. * a21)))) / (6. * pow(3. - 5. * a21, 4.) * a21);

        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n0, Lt_n0, L_n0, u_n0); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n0);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n0[ii] = dfdxcent(Lttx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 * dt * dt / 2. * Lt_n0[ii] + a21 * a21 * a21 * dt * dt * dt / 6. * Ltt_n0[ii];
            }
            j++;
        }
        borderfun3(u_n1);

        uu_to_pc(u_n1);
        Splitting(u_n1);
        Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        comput_LWqtt0(Lttx_n1, Lt_n1, L_n1, u_n1); // (double *qttt, double *qtt, double *qt, double *uuu)
        borderfun3(Lttx_n1);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Ltt_n1[ii] = dfdxcent(Lttx_n1, ii);
                // u_n2[ii] = u_n0[ii] + b1 * dt * L_n0[ii] + bb1 * dt * dt * Lt_n0[ii] + bbb1 * dt * dt * dt * Ltt_n0[ii] +
                //            b2 * dt * L_n1[ii] + bb2 * dt * dt * Lt_n1[ii] + bbb2 * dt * dt * dt * Ltt_n1[ii];
                u_n2[ii] = u_n0[ii] + a31 * dt * L_n0[ii] + aa31 * dt * dt * Lt_n0[ii] + aaa31 * dt * dt * dt * Ltt_n0[ii] +
                           a32 * dt * L_n1[ii] + aa32 * dt * dt * Lt_n1[ii] + aaa32 * dt * dt * dt * Ltt_n1[ii];
            }
            j++;
        }
        borderfun3(u_n2);

        uu_to_pc(u_n2);
        Splitting(u_n2);
        Splitting_flux_char(u_n2);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n2[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n2, L_n2, u_n2); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n2);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n2[ii] = dfdxcent(Ltx_n2, ii);
                // u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }
        // comput_LWqtt0(Lttx_n2, Lt_n2, L_n2, u_n2); // (double *qttt, double *qtt, double *qt, double *uuu)
        // borderfun3(Lttx_n1);

        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // Ltt_n1[ii] = dfdxcent(Lttx_n1, ii);
                u_nn[ii] = u_n0[ii] + b1 * dt * L_n0[ii] + bb1 * dt * dt * Lt_n0[ii] + bbb1 * dt * dt * dt * Ltt_n0[ii] +
                           b2 * dt * L_n1[ii] + bb2 * dt * dt * Lt_n1[ii] + bbb2 * dt * dt * dt * Ltt_n1[ii] +
                           b3 * dt * L_n2[ii] + bb3 * dt * dt * Lt_n2[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }

    ////-----------------------------------------
    void RK33_compt()
    {
        int j, i, k, rr;
        uu_to_pc(u_n0);
        Splitting(u_n0);
        // Splitting_flux_char(u_n0);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                df[rr] = -df_dx(fz, ff, rr);
                u_n1[rr] = u_n0[rr] + dt * df[rr];
            }
            j++;
        }
        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        // Splitting_flux_char(u_n1);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1; // df[rr] = -df_dx(fz, ff, rr);
                                 // df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                df[rr] = -df_dx(fz, ff, rr);
                u_n2[rr] = 0.75 * u_n0[rr] + 0.25 * (u_n1[rr] + dt * df[rr]);
            }
            j++;
        }
        borderfun3(u_n2);
        uu_to_pc(u_n2);
        Splitting(u_n2);
        // Splitting_flux_char(u_n2);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1; // df[rr] = -df_dx(fz, ff, rr);
                                 // df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                df[rr] = -df_dx(fz, ff, rr);
                u_nn[rr] = 1.0 / 3.0 * u_n0[rr] + 2.0 / 3.0 * (u_n2[rr] + dt * df[rr]);
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void RK33_compt_Char()
    {
        int j, i, k, rr;
        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // df[rr] = -df_dx(fz, ff, rr);
                u_n1[rr] = u_n0[rr] + dt * df[rr];
            }
            j++;
        }
        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        Splitting_flux_char(u_n1);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // df[rr] = -df_dx(fz, ff, rr);
                u_n2[rr] = 0.75 * u_n0[rr] + 0.25 * (u_n1[rr] + dt * df[rr]);
            }
            j++;
        }
        borderfun3(u_n2);
        uu_to_pc(u_n2);
        Splitting(u_n2);
        Splitting_flux_char(u_n2);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // df[rr] = -df_dx(fz, ff, rr);
                u_nn[rr] = 1.0 / 3.0 * u_n0[rr] + 2.0 / 3.0 * (u_n2[rr] + dt * df[rr]);
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void SSPRK54_compt_Shu() // Strong stability preserving Runge-Kutta and multistep time discretizations P23
    {
        int j, i, k, rr;
        uu_to_pc(u_n0);
        Splitting(u_n0);
        // Splitting_flux_char(u_n0);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                df[rr] = -df_dx(fz, ff, rr);
                u_n1[rr] = u_n0[rr] + 0.391752226571890 * dt * df[rr];
            }
            j++;
        }

        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        // Splitting_flux_char(u_n1);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                df[rr] = -df_dx(fz, ff, rr);
                u_n2[rr] = 0.444370493651235 * u_n0[rr] + 0.555629506348765 * u_n1[rr] + 0.368410593050371 * dt * df[rr];
            }
            j++;
        }

        borderfun3(u_n2);
        uu_to_pc(u_n2);
        Splitting(u_n2);
        // Splitting_flux_char(u_n2);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                df[rr] = -df_dx(fz, ff, rr);
                u_n3[rr] = 0.620101851488403 * u_n0[rr] + 0.379898148511597 * u_n2[rr] + 0.251891774271694 * dt * df[rr];
            }
            j++;
        }

        borderfun3(u_n3);
        uu_to_pc(u_n3);
        Splitting(u_n3);
        // Splitting_flux_char(u_n3);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // L_n0[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                L_n0[rr] = -df_dx(fz, ff, rr);
                u_n4[rr] = 0.178079954393132 * u_n0[rr] + 0.821920045606868 * u_n3[rr] + 0.544974750228521 * dt * L_n0[rr];
            }
            j++;
        }

        borderfun3(u_n4);
        uu_to_pc(u_n4);
        Splitting(u_n4);
        // Splitting_flux_char(u_n4);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                df[rr] = -df_dx(fz, ff, rr);
                // u_nn[rr] = 0.00683325884039 * u_n0[rr] + 0.51723167208978 * u_n2[rr] + 0.12759831133288 * u_n3[rr] +
                //            0.34833675773694 * u_n4[rr] + 0.08460416338212 * dt * L_n0[rr] + 0.22600748319395 * dt * df[rr];
                u_nn[rr] = 0.517231671970585 * u_n2[rr] + 0.096059710526147 * u_n3[rr] + 0.063692468666290 * dt * L_n0[rr] +
                           0.386708617503269 * u_n4[rr] + 0.226007483236906 * dt * df[rr];
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void SSPRK54_compt_Shu_Char() // Strong stability preserving Runge-Kutta and multistep time discretizations P23
    {
        int j, i, k, rr;
        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // df[rr] = -df_dx(fz, ff, rr);
                u_n1[rr] = u_n0[rr] + 0.391752226571890 * dt * df[rr];
            }
            j++;
        }

        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        Splitting_flux_char(u_n1);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // df[rr] = -df_dx(fz, ff, rr);
                u_n2[rr] = 0.444370493651235 * u_n0[rr] + 0.555629506348765 * u_n1[rr] + 0.368410593050371 * dt * df[rr];
            }
            j++;
        }

        borderfun3(u_n2);
        uu_to_pc(u_n2);
        Splitting(u_n2);
        Splitting_flux_char(u_n2);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // df[rr] = -df_dx(fz, ff, rr);
                u_n3[rr] = 0.620101851488403 * u_n0[rr] + 0.379898148511597 * u_n2[rr] + 0.251891774271694 * dt * df[rr];
            }
            j++;
        }

        borderfun3(u_n3);
        uu_to_pc(u_n3);
        Splitting(u_n3);
        Splitting_flux_char(u_n3);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                L_n0[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // L_n0[rr] = -df_dx(fz, ff, rr);
                u_n4[rr] = 0.178079954393132 * u_n0[rr] + 0.821920045606868 * u_n3[rr] + 0.544974750228521 * dt * L_n0[rr];
            }
            j++;
        }

        borderfun3(u_n4);
        uu_to_pc(u_n4);
        Splitting(u_n4);
        Splitting_flux_char(u_n4);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // df[rr] = -df_dx(fz, ff, rr);
                // u_nn[rr] = 0.00683325884039 * u_n0[rr] + 0.51723167208978 * u_n2[rr] + 0.12759831133288 * u_n3[rr] +
                //            0.34833675773694 * u_n4[rr] + 0.08460416338212 * dt * L_n0[rr] + 0.22600748319395 * dt * df[rr];
                u_nn[rr] = 0.517231671970585 * u_n2[rr] + 0.096059710526147 * u_n3[rr] + 0.063692468666290 * dt * L_n0[rr] +
                           0.386708617503269 * u_n4[rr] + 0.226007483236906 * dt * df[rr];
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void LSRK65_compt2() // Numerical Methods for Ordinary Differential Equations
    {
        int j = 0, i, k, rr;
        uu_to_pc(u_n0);
        Splitting(u_n0);
        // Splitting_flux_char(u_n0);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // L_n0[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                L_n0[rr] = -df_dx(fz, ff, rr);
                u_n1[rr] = u_n0[rr] + 1. / 4. * dt * L_n0[rr];
            }
            j++;
        }
        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        // Splitting_flux_char(u_n1);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // L_n1[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                L_n1[rr] = -df_dx(fz, ff, rr);
                u_n2[rr] = u_n0[rr] + 1. / 8. * dt * L_n0[rr] + 1. / 8. * dt * L_n1[rr];
                // u_n2[rr] = 0.75 * u_n0[rr] + 0.25 * (u_n1[rr] + dt * L_n1[rr]);
            }
            j++;
        }
        borderfun3(u_n2);
        uu_to_pc(u_n2);
        Splitting(u_n2);
        // Splitting_flux_char(u_n2);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // L_n2[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                L_n2[rr] = -df_dx(fz, ff, rr);
                u_n3[rr] = u_n0[rr] + 0. * dt * L_n0[rr] + 0. * dt * L_n1[rr] +
                           0.5 * dt * L_n2[rr];
            }
            j++;
        }
        borderfun3(u_n3);
        uu_to_pc(u_n3);
        Splitting(u_n3);
        // Splitting_flux_char(u_n3);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // Lt_n0[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                Lt_n0[rr] = -df_dx(fz, ff, rr);
                u_n4[rr] = u_n0[rr] + 3. / 16. * dt * L_n0[rr] - 3. / 8. * dt * L_n1[rr] +
                           3. / 8. * dt * L_n2[rr] + 9. / 16. * dt * Lt_n0[rr];
            }
            j++;
        }
        borderfun3(u_n4);
        uu_to_pc(u_n4);
        Splitting(u_n4);
        // Splitting_flux_char(u_n4);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // Lt_n1[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                Lt_n1[rr] = -df_dx(fz, ff, rr);
                u_n5[rr] = u_n0[rr] - 3. / 7. * dt * L_n0[rr] + 8. / 7. * dt * L_n1[rr] + 6. / 7. * dt * L_n2[rr] -
                           12. / 7. * dt * Lt_n0[rr] + 8. / 7. * dt * Lt_n1[rr];
            }
            j++;
        }
        borderfun3(u_n5);
        uu_to_pc(u_n5);
        Splitting(u_n5);
        // Splitting_flux_char(u_n5);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                // df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                df[rr] = -df_dx(fz, ff, rr);
                u_nn[rr] = u_n0[rr] + 7. / 90. * dt * L_n0[rr] + 0. * dt * L_n1[rr] +
                           16. / 45. * dt * L_n2[rr] + 2. / 15. * dt * Lt_n0[rr] +
                           16. / 45. * dt * Lt_n1[rr] + 7. / 90. * dt * df[rr];
                // u_nn[rr] = 1.0 / 3.0 * u_n0[rr] + 2.0 / 3.0 * (u_n2[rr] + dt * L_n2[rr]);
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void LSRK65_compt2_Char() // Numerical Methods for Ordinary Differential Equations
    {
        int j = 0, i, k, rr;
        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                L_n0[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // L_n0[rr] = -df_dx(fz, ff, rr);
                u_n1[rr] = u_n0[rr] + 1. / 4. * dt * L_n0[rr];
            }
            j++;
        }
        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        Splitting_flux_char(u_n1);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                L_n1[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // L_n1[rr] = -df_dx(fz, ff, rr);
                u_n2[rr] = u_n0[rr] + 1. / 8. * dt * L_n0[rr] + 1. / 8. * dt * L_n1[rr];
                // u_n2[rr] = 0.75 * u_n0[rr] + 0.25 * (u_n1[rr] + dt * L_n1[rr]);
            }
            j++;
        }
        borderfun3(u_n2);
        uu_to_pc(u_n2);
        Splitting(u_n2);
        Splitting_flux_char(u_n2);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                L_n2[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // L_n2[rr] = -df_dx(fz, ff, rr);
                u_n3[rr] = u_n0[rr] + 0. * dt * L_n0[rr] + 0. * dt * L_n1[rr] +
                           0.5 * dt * L_n2[rr];
            }
            j++;
        }
        borderfun3(u_n3);
        uu_to_pc(u_n3);
        Splitting(u_n3);
        Splitting_flux_char(u_n3);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                Lt_n0[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // Lt_n0[rr] = -df_dx(fz, ff, rr);
                u_n4[rr] = u_n0[rr] + 3. / 16. * dt * L_n0[rr] - 3. / 8. * dt * L_n1[rr] +
                           3. / 8. * dt * L_n2[rr] + 9. / 16. * dt * Lt_n0[rr];
            }
            j++;
        }
        borderfun3(u_n4);
        uu_to_pc(u_n4);
        Splitting(u_n4);
        Splitting_flux_char(u_n4);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                Lt_n1[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // Lt_n1[rr] = -df_dx(fz, ff, rr);
                u_n5[rr] = u_n0[rr] - 3. / 7. * dt * L_n0[rr] + 8. / 7. * dt * L_n1[rr] + 6. / 7. * dt * L_n2[rr] -
                           12. / 7. * dt * Lt_n0[rr] + 8. / 7. * dt * Lt_n1[rr];
            }
            j++;
        }
        borderfun3(u_n5);
        uu_to_pc(u_n5);
        Splitting(u_n5);
        Splitting_flux_char(u_n5);
        j = nb, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                rr = j + k * n1;
                df[rr] = -(flux[rr] - flux[rr - 1]) / dx;
                // df[rr] = -df_dx(fz, ff, rr);
                u_nn[rr] = u_n0[rr] + 7. / 90. * dt * L_n0[rr] + 0. * dt * L_n1[rr] +
                           16. / 45. * dt * L_n2[rr] + 2. / 15. * dt * Lt_n0[rr] +
                           16. / 45. * dt * Lt_n1[rr] + 7. / 90. * dt * df[rr];
                // u_nn[rr] = 1.0 / 3.0 * u_n0[rr] + 2.0 / 3.0 * (u_n2[rr] + dt * L_n2[rr]);
            }
            j++;
        }
        borderfun3(u_nn);
    }

    // -------------------------
    void compt23_LLttnn()
    {
        int j = 0, i, k, ii;
        a21 = 0.594223212099088, v2 = 0.306027487008159;
        ww1 = 0.0, ww2 = 0.0, w11 = 0.0, w22 = 0.0;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;
        uu_to_pc(u_n0);
        Splitting(u_n0);

        // Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }

        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        // Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); //  Ltx_n1[j] = -L_n1[j] * u_n1[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                u_nn[ii] = u_n0[ii] + dt * (v1 * L_n0[ii] + v2 * L_n1[ii] + w11 * TL_n0[ii] + w22 * TL_n1[ii]) +
                           dt * dt * (vv1 * Lt_n0[ii] + vv2 * Lt_n1[ii] + ww1 * TLt_n0[ii] + ww2 * TLt_n1[ii]);
                // u_nn[ii] = u_n0[ii] + dt * L_n0[ii] + dt * dt / 2. * Lt_n0[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void compt23_LLttnn_Char()
    {
        int j = 0, i, k, ii;
        a21 = 0.594223212099088, v2 = 0.306027487008159;
        ww1 = 0.0, ww2 = 0.0, w11 = 0.0, w22 = 0.0;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;
        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }

        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); //  Ltx_n1[j] = -L_n1[j] * u_n1[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                u_nn[ii] = u_n0[ii] + dt * (v1 * L_n0[ii] + v2 * L_n1[ii] + w11 * TL_n0[ii] + w22 * TL_n1[ii]) +
                           dt * dt * (vv1 * Lt_n0[ii] + vv2 * Lt_n1[ii] + ww1 * TLt_n0[ii] + ww2 * TLt_n1[ii]);
                // u_nn[ii] = u_n0[ii] + dt * L_n0[ii] + dt * dt / 2. * Lt_n0[ii];
            }
            j++;
        }
        borderfun3(u_nn);
    }

    void compt24_LLttnn()
    {
        int j = 0, i, k, ii;

        a21 = 0.5, v1 = 1., v2 = 0., vv1 = 1. / 6., vv2 = 1. / 3.; // C0 = 1. / 3., C1 = 2. / 3.;
        ww1 = 0.0, ww2 = 0.0, w11 = 0.0, w22 = 0.0;
        uu_to_pc(u_n0);
        Splitting(u_n0);
        // Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }

        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        // Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); //  Ltx_n1[j] = -L_n1[j] * u_n1[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                u_nn[ii] = u_n0[ii] + dt * (v1 * L_n0[ii] + v2 * L_n1[ii]) +
                           dt * dt * (vv1 * Lt_n0[ii] + vv2 * Lt_n1[ii]);
            }

            j++;
        }
        borderfun3(u_nn);
    }

    void compt24_LLttnn_Char()
    {
        int j = 0, i, k, ii;

        a21 = 0.5, v1 = 1., v2 = 0., vv1 = 1. / 6., vv2 = 1. / 3.; // C0 = 1. / 3., C1 = 2. / 3.;
        ww1 = 0.0, ww2 = 0.0, w11 = 0.0, w22 = 0.0;
        uu_to_pc(u_n0);
        Splitting(u_n0);
        Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }

        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                // L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); //  Ltx_n1[j] = -L_n1[j] * u_n1[j];
        borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                u_nn[ii] = u_n0[ii] + dt * (v1 * L_n0[ii] + v2 * L_n1[ii]) +
                           dt * dt * (vv1 * Lt_n0[ii] + vv2 * Lt_n1[ii]);
            }

            j++;
        }
        borderfun3(u_nn);
    }

    void compt35_LLttnn()
    {
        int j = 0, i, k, ii;

        a21 = 0.1; // a21 = 0.752;
        aa12 = -1. + 2. * a21;
        vv1 = (1. + 2. * a21 * (-4. + 5. * a21)) / (12. * a21 * (-3. + 5. * a21));
        vv2 = 1. / (12. * a21 * (3. + 10. * (-1. + a21) * a21));
        vv3 = (25. * pow(aa12, 3)) / (12. * (-3. + 5. * a21) * (3. + 10. * (-1 + a21) * a21));
        v1 = 1., v2 = 0., v3 = 0.;
        a31 = (3. - 5. * a21) / (5. - 10. * a21), a32 = 0.;
        aa31 = ((-3. + 5. * a21) * (-3. + 10. * a21) * (1. + 5. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));
        aa32 = ((-3. + 5. * a21) * (3. + 10. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));
        uu_to_pc(u_n0);
        Splitting(u_n0);
        // Splitting_flux_char(u_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n0[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n0[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n0, L_n0, u_n0); // Ltx_n0[j] = -L_n0[j] * u_n0[j];
        // borderfun3(Ltx_n0);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n0[ii] = dfdxcent(Ltx_n0, ii);
                u_n1[ii] = u_n0[ii] + a21 * dt * L_n0[ii] + a21 * a21 / 2.0 * dt * dt * Lt_n0[ii];
            }
            j++;
        }

        borderfun3(u_n1);
        uu_to_pc(u_n1);
        Splitting(u_n1);
        // Splitting_flux_char(u_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n1[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n1[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n1, L_n1, u_n1); //  Ltx_n1[j] = -L_n1[j] * u_n1[j];
        // borderfun3(Ltx_n1);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n1[ii] = dfdxcent(Ltx_n1, ii);
                u_n2[ii] = u_n0[ii] + dt * (a31 * L_n0[ii] + a32 * L_n1[ii]) +
                           dt * dt * (aa31 * Lt_n0[ii] + aa32 * Lt_n1[ii]);
            }
            j++;
        }
        borderfun3(u_n2);
        uu_to_pc(u_n2);
        Splitting(u_n2);
        // Splitting_flux_char(u_n2);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                // L_n2[ii] = -(flux[ii] - flux[ii - 1]) / dx;
                L_n2[ii] = -df_dx(fz, ff, ii);
            }
            j++;
        }
        comput_LWqt0(Ltx_n2, L_n2, u_n2); //  Ltx_n1[j] = -L_n1[j] * u_n1[j];
        // borderfun3(Ltx_n2);
        j = 0, i = n1 - nb * 2;
        while (i-- > 0)
        {
            for (int k = 0; k < 3; k++)
            {
                ii = nb + j + k * n1;
                Lt_n2[ii] = dfdxcent(Ltx_n2, ii);
                u_nn[ii] = u_n0[ii] + dt * (v1 * L_n0[ii] + v2 * L_n1[ii] + v3 * L_n2[ii]) +
                           dt * dt * (vv1 * Lt_n0[ii] + vv2 * Lt_n1[ii] + vv3 * Lt_n2[ii]);
            }
            j++;
        }
        borderfun3(u_nn);
        uu_to_pc(u_nn);
    }
 
    //-----------------------------------------
    int Perimeter()
    {
        cout << " yes ";
        return 2 * t;
    }

    double df_dx(double *fz, double *ff, int ii)
    {
        double y, ss1, ss2;
        int i = ii;
        // 1up
        //  ss1 = (*(fz + i) - *(fz + i - 1)) / dx; //正通量
        //  ss2 = (*(ff + i + 1) - *(ff + i)) / dx;  //负通量
        // 8up
        // ss1 = (-*(ff + i - 3) * 5. + *(ff + i - 2) * 60. - *(ff + i - 1) * 420. - *(ff + i) * 378 + *(ff + i + 1) * 1050. - *(ff + i + 2) * 420. + *(ff + i + 3) * 140. - *(ff + i + 4) * 30. + *(ff + i + 5) * 3.) / 840. / dx;
        // ss2 = -(-*(fz + i + 3) * 5. + *(fz + i + 2) * 60. - *(fz + i + 1) * 420. - *(fz + i) * 378 + *(fz + i - 1) * 1050. - *(fz + i - 2) * 420. + *(fz + i - 3) * 140. - *(fz + i - 4) * 30. + *(fz + i - 5) * 3.) / 840. / dx;

        // 9up
        //  ss1 = (*(ff + i - 4) - *(ff + i - 3) * 12. + *(ff + i - 2) * 72. - *(ff + i - 1) * 336. - *(ff + i) * 100.8 + *(ff + i + 1) * 504. - *(ff + i + 2) * 168. + *(ff + i + 3) * 48. - *(ff + i + 4) * 9. + *(ff + i + 5) * 0.8) / 504. / dx;
        //  ss2 = -(*(fz + i + 4) - *(fz + i + 3) * 12. + *(fz + i + 2) * 72. - *(fz + i + 1) * 336. - *(fz + i) * 100.8 + *(fz + i - 1) * 504. - *(fz + i - 2) * 168. + *(fz + i - 3) * 48. - *(fz + i - 4) * 9. + *(fz + i - 5) * 0.8) / 504. / dx;

        /*  //weno_5-js  正通量
        double ep = 1.E-5, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(fz + i - 2) - 2. * (*(fz + i - 1)) + *(fz + i)), 2) + 0.25 * pow((*(fz + i - 2) - 4. * (*(fz + i - 1)) + 3. * (*(fz + i))), 2);
            sn1 = 13. / 12. * pow((*(fz + i - 1) - 2. * (*(fz + i)) + *(fz + i + 1)), 2) + 0.25 * pow((*(fz + i - 1) - *(fz + i + 1)), 2);
            sn2 = 13. / 12. * pow((*(fz + i) - 2. * (*(fz + i + 1)) + *(fz + i + 2)), 2) + 0.25 * pow((3. * (*(fz + i)) - 4. * (*(fz + i + 1)) + *(fz + i + 2)), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 2. / 6. * (*(fz + i - 2)) - 7. / 6. * (*(fz + i - 1)) + 11. / 6. * (*(fz + i));
            q13 = -1. / 6. * (*(fz + i - 1)) + 5. / 6. * (*(fz + i)) + 2. / 6. * (*(fz + i + 1));
            q23 = 2. / 6. * (*(fz + i)) + 5. / 6. * (*(fz + i + 1)) - 1. / 6. * (*(fz + i + 2));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        ss1 = (vz[0] - vz[1]) / dx;
        // weno_5-js 负通量
        //  double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        //  double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        i = ii, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(ff + i + 1) - 2. * (*(ff + i + 2)) + *(ff + i + 3)), 2) + 0.25 * pow((3. * (*(ff + i + 1)) - 4. * (*(ff + i + 2)) + *(ff + i + 3)), 2);
            sn1 = 13. / 12. * pow((*(ff + i) - 2. * (*(ff + i + 1)) + *(ff + i + 2)), 2) + 0.25 * pow((*(ff + i) - *(ff + i + 2)), 2);
            sn2 = 13. / 12. * pow((*(ff + i - 1) - 2. * (*(ff + i)) + *(ff + i + 1)), 2) + 0.25 * pow((*(ff + i - 1) - 4. * (*(ff + i)) + 3. * (*(ff + i + 1))), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 11. / 6. * (*(ff + i + 1)) - 7. / 6. * (*(ff + i + 2)) + 2. / 6. * (*(ff + i + 3));
            q13 = 2. / 6. * (*(ff + i)) + 5. / 6. * (*(ff + i + 1)) - 1. / 6. * (*(ff + i + 2));
            q23 = -1. / 6. * (*(ff + i - 1)) + 5. / 6. * (*(ff + i)) + 2. / 6. * (*(ff + i + 1));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        ss2 = (vz[0] - vz[1]) / dx;
         */

        /* //weno-z
        double pe = 2.0, epb = 1.E-6, tau, sum_a, aaz[3], wwz[3], qqz[3], b[3];
        double cc0 = 1.0 / 10.0, cc1 = 3.0 / 5.0, cc2 = 3.0 / 10.0, ccz[3], vz[2]; //!理想权重
        ccz[0] = cc0, ccz[1] = cc1, ccz[2] = cc2;
        for (int j = 0; j < 2; j++)
        {
            b[2] = 13. / 12. * pow((*(fz + i) - 2. * (*(fz + i + 1)) + *(fz + i + 2)), 2) + 0.25 * pow((3. * (*(fz + i)) - 4. * (*(fz + i + 1)) + *(fz + i + 2)), 2);
            b[1] = 13. / 12. * pow((*(fz + i - 1) - 2. * (*(fz + i)) + *(fz + i + 1)), 2) + 0.25 * pow((*(fz + i - 1) - *(fz + i + 1)), 2);
            b[0] = 13. / 12. * pow((*(fz + i - 2) - 2. * (*(fz + i - 1)) + *(fz + i)), 2) + 0.25 * pow((*(fz + i - 2) - 4. * (*(fz + i - 1)) + 3. * (*(fz + i))), 2);
            tau = abs(b[2] - b[0]);
            //   lamdx=dx**(2.0/3.0)
            //   tau=abs(b2-b0)
            //   a=c*(1.0+((tau+ept)/(b+epb))**pw+lamdx*(b+epb)/(tau+ept))
            for (int jj = 0; jj < 3; jj++)
            {
                aaz[jj] = ccz[jj] * (1.0 + pow((tau / (b[jj] + epb)), 2));
            }
            sum_a = aaz[0] + aaz[1] + aaz[2];
            wwz[0] = aaz[0] / sum_a, wwz[1] = aaz[1] / sum_a, wwz[2] = aaz[2] / sum_a; //!正式权重
            qqz[0] = 1.0 / 3.0 * (*(fz + i - 2)) - 7.0 / 6.0 * (*(fz + i - 1)) + 11.0 / 6.0 * (*(fz + i));
            qqz[1] = -1.0 / 6.0 * (*(fz + i - 1)) + 5.0 / 6.0 * (*(fz + i)) + 1.0 / 3.0 * (*(fz + i + 1));
            qqz[2] = 1.0 / 3.0 * (*(fz + i)) + 5.0 / 6.0 * (*(fz + i + 1)) - 1.0 / 6.0 * (*(fz + i + 2));

            vz[j] = wwz[0] * qqz[0] + wwz[1] * qqz[1] + wwz[2] * qqz[2];
            i--;
        }
        ss1 = (vz[0] - vz[1]) / dx;

        i = ii;
        for (int j = 0; j < 2; j++)
        {
            b[0] = 13. / 12. * pow((*(ff + i + 1) - 2. * (*(ff + i + 2)) + *(ff + i + 3)), 2) + 0.25 * pow((3. * (*(ff + i + 1)) - 4. * (*(ff + i + 2)) + *(ff + i + 3)), 2);
            b[1] = 13. / 12. * pow((*(ff + i) - 2. * (*(ff + i + 1)) + *(ff + i + 2)), 2) + 0.25 * pow((*(ff + i) - *(ff + i + 2)), 2);
            b[2] = 13. / 12. * pow((*(ff + i - 1) - 2. * (*(ff + i)) + *(ff + i + 1)), 2) + 0.25 * pow((*(ff + i - 1) - 4. * (*(ff + i)) + 3. * (*(ff + i + 1))), 2);
            tau = abs(b[2] - b[0]);
            for (int jj = 0; jj < 3; jj++)
            {
                aaz[jj] = ccz[jj] * (1.0 + pow((tau / (b[jj] + epb)), 2));
            }
            sum_a = aaz[0] + aaz[1] + aaz[2];
            wwz[0] = aaz[0] / sum_a, wwz[1] = aaz[1] / sum_a, wwz[2] = aaz[2] / sum_a; //!正式权重
            qqz[0] = 11. / 6. * (*(ff + i + 1)) - 7. / 6. * (*(ff + i + 2)) + 2. / 6. * (*(ff + i + 3));
            qqz[1] = 2. / 6. * (*(ff + i)) + 5. / 6. * (*(ff + i + 1)) - 1. / 6. * (*(ff + i + 2));
            qqz[2] = -1. / 6. * (*(ff + i - 1)) + 5. / 6. * (*(ff + i)) + 2. / 6. * (*(ff + i + 1));

            vz[j] = wwz[0] * qqz[0] + wwz[1] * qqz[1] + wwz[2] * qqz[2];
            i--;
        }
        ss2 = (vz[0] - vz[1]) / dx;
         */

        /*WEMO_9

        double dd[4], da1[4], aa1[4][4], bb1[4][4], cc1[4][4], dd1[4][4], b[4], aerfa[4], q[4], w[4], vz[2];

        double ep = 1.0E-6, ss11 = 0.0, ss12 = 0.0, ss13 = 0.0, p = 2.0, all_aerfa;
        int j, k, mm;

        dd[0] = 1.0 / 35.0, dd[1] = 12.0 / 35.0, dd[2] = 18.0 / 35.0, dd[3] = 4.0 / 35.0;

        aa1[0][0] = -3.0 / 12.0, aa1[0][1] = 13.0 / 12.0, aa1[0][2] = -23.0 / 12.0, aa1[0][3] = 25.0 / 12.0;
        aa1[1][0] = 1.0 / 12.0, aa1[1][1] = -5.0 / 12.0, aa1[1][2] = 13.0 / 12.0, aa1[1][3] = 3.0 / 12.0;
        aa1[2][0] = -1.0 / 12.0, aa1[2][1] = 7.0 / 12.0, aa1[2][2] = 7.0 / 12.0, aa1[2][3] = -1.0 / 12.0;
        aa1[3][0] = 3.0 / 12.0, aa1[3][1] = 13.0 / 12.0, aa1[3][2] = -5.0 / 12.0, aa1[3][3] = 1.0 / 12.0;

        bb1[0][0] = -2.0 / 6.0, bb1[0][1] = 21.0 / 6.0, bb1[0][2] = -18.0 / 6.0, bb1[0][3] = 11.0 / 6.0;
        bb1[1][0] = 1.0 / 6.0, bb1[1][1] = -6.0 / 6.0, bb1[1][2] = 3.0 / 6.0, bb1[1][3] = 2.0 / 6.0;
        bb1[2][0] = -2.0 / 6.0, bb1[2][1] = -3.0 / 6.0, bb1[2][2] = 6.0 / 6.0, bb1[2][3] = -1.0 / 6.0;
        bb1[3][0] = -11.0 / 6.0, bb1[3][1] = 18.0 / 6.0, bb1[3][2] = -21.0 / 6.0, bb1[3][3] = 2.0 / 6.0;

        cc1[0][0] = -1.0, cc1[0][1] = 4.0, cc1[0][2] = -5.0, cc1[0][3] = 2.0;
        cc1[1][0] = 0.0, cc1[1][1] = 1.0, cc1[1][2] = -2.0, cc1[1][3] = 1.0;
        cc1[2][0] = 1.0, cc1[2][1] = -2.0, cc1[2][2] = 1.0, cc1[2][3] = 0.0;
        cc1[3][0] = 2.0, cc1[3][1] = -5.0, cc1[3][2] = 4.0, cc1[3][3] = -1.0;

        dd1[0][0] = -1.0, dd1[0][1] = 3.0, dd1[0][2] = -3.0, dd1[0][3] = 1.0;
        dd1[1][0] = -1.0, dd1[1][1] = 3.0, dd1[1][2] = -3.0, dd1[1][3] = 1.0;
        dd1[2][0] = -1.0, dd1[2][1] = 3.0, dd1[2][2] = -3.0, dd1[2][3] = 1.0;
        dd1[3][0] = -1.0, dd1[3][1] = 3.0, dd1[3][2] = -3.0, dd1[3][3] = 1.0;

        da1[0] = 1.0 / 35., da1[1] = 12. / 35., da1[2] = 18. / 35., da1[3] = 4. / 35.;
        for (mm = 1; mm > -0.5; mm--)
        {

            all_aerfa = 0.0;
            for (k = 0; k < 4; k++)
            {
                ss11 = 0.0, ss12 = 0.0, ss13 = 0.0, q[k] = 0;

                for (j = 0; j < 4; j++)
                {
                    ss11 = ss11 + bb1[k][j] * (*(fz + i + k + j - 3));
                    ss12 = ss12 + cc1[k][j] * (*(fz + i + k + j - 3));
                    ss13 = ss13 + dd1[k][j] * (*(fz + i + k + j - 3));
                    q[k] = q[k] + aa1[k][j] * (*(fz + i + k + j - 3)); //f_j+0.5
                }
                b[k] = ss11 * ss11 + 13.0 / 12.0 * ss12 * ss12 + 1043.0 / 2160.0 * ss13 * ss13 + 1.0 / 12.0 * ss11 * ss13;

                //b[k] = 0;
                aerfa[k] = da1[k] / (b[k] + ep) / (b[k] + ep);
                all_aerfa = all_aerfa + aerfa[k];
            }
            vz[mm] = 0;
            for (k = 0; k < 4; k++)
            {
                w[k] = aerfa[k] / all_aerfa;
                vz[mm] = vz[mm] + w[k] * q[k];
            }
            i--;
        }
        ss1 = (vz[1] - vz[0]) / dx;
         */

        // /*   weno7
        //     subroutine hh_weno7P(Ka,Kb,v,hh )           ! Ka=-3,  Kb=3
        //  Use OCFD_constants
        //  implicit none
        double S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33,
            ax0, ax1, ax2, ax3, am, q0, q1, q2, q3, vz[2];
        double CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
               ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
               ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
               ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
               ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
               b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
               b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
               d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
        double e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
               e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
               e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
               e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-4; //   !! WENO-JS
        for (int j = 0; j < 2; j++)
        {
            // ! 7th order WENO scheme
            // ! 1  阶导数      *(fz + i - 3)
            S10 = ax11 * *(fz + i - 3) + ax12 * *(fz + i - 2) + ax13 * *(fz + i - 1) + ax14 * *(fz + i);
            S11 = ax21 * *(fz + i - 2) - *(fz + i - 1) + ax23 * *(fz + i) + ax24 * *(fz + i + 1);
            S12 = ax31 * *(fz + i - 1) + ax32 * *(fz + i) + *(fz + i + 1) + ax34 * *(fz + i + 2);
            S13 = ax41 * *(fz + i) + ax42 * *(fz + i + 1) + ax43 * *(fz + i + 2) + ax44 * *(fz + i + 3);
            //  ! 2 阶导数
            S20 = -*(fz + i - 3) + b12 * *(fz + i - 2) + b13 * *(fz + i - 1) + b14 * *(fz + i);
            S21 = *(fz + i - 1) + b22 * *(fz + i) + *(fz + i + 1);
            S22 = *(fz + i) + b22 * *(fz + i + 1) + *(fz + i + 2);
            S23 = b41 * *(fz + i) + b42 * *(fz + i + 1) + b43 * *(fz + i + 2) - *(fz + i + 3);
            // ! 3 阶导数
            S30 = -*(fz + i - 3) + c12 * (*(fz + i - 2) - *(fz + i - 1)) + *(fz + i);
            S31 = -*(fz + i - 2) + c12 * (*(fz + i - 1) - *(fz + i)) + *(fz + i + 1);
            S32 = -*(fz + i - 1) + c12 * (*(fz + i) - *(fz + i + 1)) + *(fz + i + 2);
            S33 = -*(fz + i) + c12 * (*(fz + i + 1) - *(fz + i + 2)) + *(fz + i + 3);

            S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
            S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
            S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
            S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

            // !-------WENO J-S----------------------
            ax0 = CC0 / ((ep + S0) * (ep + S0));
            ax1 = CC1 / ((ep + S1) * (ep + S1));
            ax2 = CC2 / ((ep + S2) * (ep + S2));
            ax3 = CC3 / ((ep + S3) * (ep + S3));
            // !-----------------------------------------------
            am = ax0 + ax1 + ax2 + ax3;

            // !  4阶差分格式的通量
            q0 = e11 * *(fz + i - 3) + e12 * *(fz + i - 2) + e13 * *(fz + i - 1) + e14 * *(fz + i);
            q1 = e21 * *(fz + i - 2) + e22 * *(fz + i - 1) + e23 * *(fz + i) + e24 * *(fz + i + 1);
            q2 = e31 * *(fz + i - 1) + e32 * *(fz + i) + e33 * *(fz + i + 1) + e34 * *(fz + i + 2);
            q3 = e41 * *(fz + i) + e42 * *(fz + i + 1) + e43 * *(fz + i + 2) + e44 * *(fz + i + 3);

            // !  由4个4阶差分格式组合成1个7阶差分格式
            // !     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3;
            vz[j] = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;
            i--;
        }
        ss1 = (vz[0] - vz[1]) / dx;
        // -------------------------------
        i = ii; ////////////////////////
        // CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.;
        // ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.;
        // ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.;
        // ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.;
        // ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.;
        // b12 = 4., b13 = -5., b14 = 2., b22 = -2.;
        // b41 = 2., b42 = -5., b43 = 4., c12 = 3.;
        // d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
        // e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.;
        // e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.;
        // e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.;
        // e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.E-10; //

        for (int j = 0; j < 2; j++)
        {
            // !      7th order WENO scheme
            // ! 1  阶导数   *(ff + i - 2)
            S10 = ax11 * *(ff + i + 4) + ax12 * *(ff + i + 3) + ax13 * *(ff + i + 2) + ax14 * *(ff + i + 1);
            S11 = ax21 * *(ff + i + 3) - *(ff + i + 2) + ax23 * *(ff + i + 1) + ax24 * *(ff + i);
            S12 = ax31 * *(ff + i + 2) + ax32 * *(ff + i + 1) + *(ff + i) + ax34 * *(ff + i - 1);
            S13 = ax41 * *(ff + i + 1) + ax42 * *(ff + i) + ax43 * *(ff + i - 1) + ax44 * *(ff + i - 2);
            // ! 2 阶导数
            S20 = -*(ff + i + 4) + b12 * *(ff + i + 3) + b13 * *(ff + i + 2) + b14 * *(ff + i + 1);
            S21 = *(ff + i + 2) + b22 * *(ff + i + 1) + *(ff + i);
            S22 = *(ff + i + 1) + b22 * *(ff + i) + *(ff + i - 1);
            S23 = b41 * *(ff + i + 1) + b42 * *(ff + i) + b43 * *(ff + i - 1) - *(ff + i - 2);
            // ! 3 阶导数
            S30 = -*(ff + i + 4) + c12 * (*(ff + i + 3) - *(ff + i + 2)) + *(ff + i + 1);
            S31 = -*(ff + i + 3) + c12 * (*(ff + i + 2) - *(ff + i + 1)) + *(ff + i);
            S32 = -*(ff + i + 2) + c12 * (*(ff + i + 1) - *(ff + i)) + *(ff + i - 1);
            S33 = -*(ff + i + 1) + c12 * (*(ff + i) - *(ff + i - 1)) + *(ff + i - 2);

            S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
            S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
            S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
            S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

            ax0 = CC0 / ((ep + S0) * (ep + S0));
            ax1 = CC1 / ((ep + S1) * (ep + S1));
            ax2 = CC2 / ((ep + S2) * (ep + S2));
            ax3 = CC3 / ((ep + S3) * (ep + S3));

            // !-----------------------------------------------

            am = ax0 + ax1 + ax2 + ax3;

            // !  4阶差分格式的通量
            q0 = e11 * *(ff + i + 4) + e12 * *(ff + i + 3) + e13 * *(ff + i + 2) + e14 * *(ff + i + 1);
            q1 = e21 * *(ff + i + 3) + e22 * *(ff + i + 2) + e23 * *(ff + i + 1) + e24 * *(ff + i);
            q2 = e31 * *(ff + i + 2) + e32 * *(ff + i + 1) + e33 * *(ff + i) + e34 * *(ff + i - 1);
            q3 = e41 * *(ff + i + 1) + e42 * *(ff + i) + e43 * *(ff + i - 1) + e44 * *(ff + i - 2);

            // !  由4个4阶差分格式组合成1个7阶差分格式
            vz[j] = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

            i--;
        }
        ss2 = (vz[0] - vz[1]) / dx;
        // */

        /*   weno7-Z
        //     subroutine hh_weno7P(Ka,Kb,v,hh )           ! Ka=-3,  Kb=3
        //  Use OCFD_constants
        //  implicit none
        double S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33,
            ax0, ax1, ax2, ax3, am, q0, q1, q2, q3, vz[2], tau;
        double CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
               ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
               ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
               ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
               ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
               b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
               b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
               d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
        double e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
               e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
               e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
               e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-6; //   !! WENO-JS
        for (int j = 0; j < 2; j++)
        {
            // ! 7th order WENO scheme
            // ! 1  阶导数      *(fz + i - 3)
            S10 = ax11 * *(fz + i - 3) + ax12 * *(fz + i - 2) + ax13 * *(fz + i - 1) + ax14 * *(fz + i);
            S11 = ax21 * *(fz + i - 2) - *(fz + i - 1) + ax23 * *(fz + i) + ax24 * *(fz + i + 1);
            S12 = ax31 * *(fz + i - 1) + ax32 * *(fz + i) + *(fz + i + 1) + ax34 * *(fz + i + 2);
            S13 = ax41 * *(fz + i) + ax42 * *(fz + i + 1) + ax43 * *(fz + i + 2) + ax44 * *(fz + i + 3);
            //  ! 2 阶导数
            S20 = -*(fz + i - 3) + b12 * *(fz + i - 2) + b13 * *(fz + i - 1) + b14 * *(fz + i);
            S21 = *(fz + i - 1) + b22 * *(fz + i) + *(fz + i + 1);
            S22 = *(fz + i) + b22 * *(fz + i + 1) + *(fz + i + 2);
            S23 = b41 * *(fz + i) + b42 * *(fz + i + 1) + b43 * *(fz + i + 2) - *(fz + i + 3);
            // ! 3 阶导数
            S30 = -*(fz + i - 3) + c12 * (*(fz + i - 2) - *(fz + i - 1)) + *(fz + i);
            S31 = -*(fz + i - 2) + c12 * (*(fz + i - 1) - *(fz + i)) + *(fz + i + 1);
            S32 = -*(fz + i - 1) + c12 * (*(fz + i) - *(fz + i + 1)) + *(fz + i + 2);
            S33 = -*(fz + i) + c12 * (*(fz + i + 1) - *(fz + i + 2)) + *(fz + i + 3);

            S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
            S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
            S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
            S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

            // tau = abs(S0 - S1 - S2 + S3);
            tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
            // !-------WENO Z----------------------
            ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
            ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
            ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
            ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
            // !-------WENO J-S----------------------
            // ax0 = CC0 / ((ep + S0) * (ep + S0));
            // ax1 = CC1 / ((ep + S1) * (ep + S1));
            // ax2 = CC2 / ((ep + S2) * (ep + S2));
            // ax3 = CC3 / ((ep + S3) * (ep + S3));
            // !-----------------------------------------------
            am = ax0 + ax1 + ax2 + ax3;

            // !  4阶差分格式的通量
            q0 = e11 * *(fz + i - 3) + e12 * *(fz + i - 2) + e13 * *(fz + i - 1) + e14 * *(fz + i);
            q1 = e21 * *(fz + i - 2) + e22 * *(fz + i - 1) + e23 * *(fz + i) + e24 * *(fz + i + 1);
            q2 = e31 * *(fz + i - 1) + e32 * *(fz + i) + e33 * *(fz + i + 1) + e34 * *(fz + i + 2);
            q3 = e41 * *(fz + i) + e42 * *(fz + i + 1) + e43 * *(fz + i + 2) + e44 * *(fz + i + 3);

            // !  由4个4阶差分格式组合成1个7阶差分格式
            // !     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3;
            vz[j] = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;
            i--;
        }
        ss1 = (vz[0] - vz[1]) / dx;
        // -------------------------------
        i = ii; ////////////////////////

        for (int j = 0; j < 2; j++)
        {
            // !      7th order WENO scheme
            // ! 1  阶导数   *(ff + i - 2)
            S10 = ax11 * *(ff + i + 4) + ax12 * *(ff + i + 3) + ax13 * *(ff + i + 2) + ax14 * *(ff + i + 1);
            S11 = ax21 * *(ff + i + 3) - *(ff + i + 2) + ax23 * *(ff + i + 1) + ax24 * *(ff + i);
            S12 = ax31 * *(ff + i + 2) + ax32 * *(ff + i + 1) + *(ff + i) + ax34 * *(ff + i - 1);
            S13 = ax41 * *(ff + i + 1) + ax42 * *(ff + i) + ax43 * *(ff + i - 1) + ax44 * *(ff + i - 2);
            // ! 2 阶导数
            S20 = -*(ff + i + 4) + b12 * *(ff + i + 3) + b13 * *(ff + i + 2) + b14 * *(ff + i + 1);
            S21 = *(ff + i + 2) + b22 * *(ff + i + 1) + *(ff + i);
            S22 = *(ff + i + 1) + b22 * *(ff + i) + *(ff + i - 1);
            S23 = b41 * *(ff + i + 1) + b42 * *(ff + i) + b43 * *(ff + i - 1) - *(ff + i - 2);
            // ! 3 阶导数
            S30 = -*(ff + i + 4) + c12 * (*(ff + i + 3) - *(ff + i + 2)) + *(ff + i + 1);
            S31 = -*(ff + i + 3) + c12 * (*(ff + i + 2) - *(ff + i + 1)) + *(ff + i);
            S32 = -*(ff + i + 2) + c12 * (*(ff + i + 1) - *(ff + i)) + *(ff + i - 1);
            S33 = -*(ff + i + 1) + c12 * (*(ff + i) - *(ff + i - 1)) + *(ff + i - 2);

            S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
            S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
            S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
            S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

            tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
            // !-------WENO Z----------------------
            ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
            ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
            ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
            ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));

            // ax0 = CC0 / ((ep + S0) * (ep + S0));
            // ax1 = CC1 / ((ep + S1) * (ep + S1));
            // ax2 = CC2 / ((ep + S2) * (ep + S2));
            // ax3 = CC3 / ((ep + S3) * (ep + S3));

            // !-----------------------------------------------

            am = ax0 + ax1 + ax2 + ax3;

            // !  4阶差分格式的通量
            q0 = e11 * *(ff + i + 4) + e12 * *(ff + i + 3) + e13 * *(ff + i + 2) + e14 * *(ff + i + 1);
            q1 = e21 * *(ff + i + 3) + e22 * *(ff + i + 2) + e23 * *(ff + i + 1) + e24 * *(ff + i);
            q2 = e31 * *(ff + i + 2) + e32 * *(ff + i + 1) + e33 * *(ff + i) + e34 * *(ff + i - 1);
            q3 = e41 * *(ff + i + 1) + e42 * *(ff + i) + e43 * *(ff + i - 1) + e44 * *(ff + i - 2);

            // !  由4个4阶差分格式组合成1个7阶差分格式
            vz[j] = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

            i--;
        }
        ss2 = (vz[0] - vz[1]) / dx;
        */

        y = ss1 + ss2;
        return y;
    }

    double df_dx_char(double *ffz, double *fff, int ii)
    {
        double yy = 0.0, ss1, ss2;
        int i = ii;

        // /*  //weno_5-js  正通量//
        double ep = 1.E-4, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23;
        sn0 = 13. / 12. * pow((*(ffz + i - 2) - 2. * (*(ffz + i - 1)) + *(ffz + i)), 2) + 0.25 * pow((*(ffz + i - 2) - 4. * (*(ffz + i - 1)) + 3. * (*(ffz + i))), 2);
        sn1 = 13. / 12. * pow((*(ffz + i - 1) - 2. * (*(ffz + i)) + *(ffz + i + 1)), 2) + 0.25 * pow((*(ffz + i - 1) - *(ffz + i + 1)), 2);
        sn2 = 13. / 12. * pow((*(ffz + i) - 2. * (*(ffz + i + 1)) + *(ffz + i + 2)), 2) + 0.25 * pow((3. * (*(ffz + i)) - 4. * (*(ffz + i + 1)) + *(ffz + i + 2)), 2);
        az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
        W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
        q03 = 2. / 6. * (*(ffz + i - 2)) - 7. / 6. * (*(ffz + i - 1)) + 11. / 6. * (*(ffz + i));
        q13 = -1. / 6. * (*(ffz + i - 1)) + 5. / 6. * (*(ffz + i)) + 2. / 6. * (*(ffz + i + 1));
        q23 = 2. / 6. * (*(ffz + i)) + 5. / 6. * (*(ffz + i + 1)) - 1. / 6. * (*(ffz + i + 2));
        ss1 = W0 * q03 + W1 * q13 + W2 * q23;
        // weno_5-js 负通量
        //  double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        //  double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        i = ii, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        sn0 = 13. / 12. * pow((*(fff + i + 1) - 2. * (*(fff + i + 2)) + *(fff + i + 3)), 2) + 0.25 * pow((3. * (*(fff + i + 1)) - 4. * (*(fff + i + 2)) + *(fff + i + 3)), 2);
        sn1 = 13. / 12. * pow((*(fff + i) - 2. * (*(fff + i + 1)) + *(fff + i + 2)), 2) + 0.25 * pow((*(fff + i) - *(fff + i + 2)), 2);
        sn2 = 13. / 12. * pow((*(fff + i - 1) - 2. * (*(fff + i)) + *(fff + i + 1)), 2) + 0.25 * pow((*(fff + i - 1) - 4. * (*(fff + i)) + 3. * (*(fff + i + 1))), 2);
        az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
        W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
        q03 = 11. / 6. * (*(fff + i + 1)) - 7. / 6. * (*(fff + i + 2)) + 2. / 6. * (*(fff + i + 3));
        q13 = 2. / 6. * (*(fff + i)) + 5. / 6. * (*(fff + i + 1)) - 1. / 6. * (*(fff + i + 2));
        q23 = -1. / 6. * (*(fff + i - 1)) + 5. / 6. * (*(fff + i)) + 2. / 6. * (*(fff + i + 1));
        ss2 = W0 * q03 + W1 * q13 + W2 * q23;
        yy = ss1 + ss2;
        // */

        /*  //weno_5-Z  正通量//
        double ep = 1.E-5, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        double tau, sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23;
        sn0 = 13. / 12. * pow((*(ffz + i - 2) - 2. * (*(ffz + i - 1)) + *(ffz + i)), 2) + 0.25 * pow((*(ffz + i - 2) - 4. * (*(ffz + i - 1)) + 3. * (*(ffz + i))), 2);
        sn1 = 13. / 12. * pow((*(ffz + i - 1) - 2. * (*(ffz + i)) + *(ffz + i + 1)), 2) + 0.25 * pow((*(ffz + i - 1) - *(ffz + i + 1)), 2);
        sn2 = 13. / 12. * pow((*(ffz + i) - 2. * (*(ffz + i + 1)) + *(ffz + i + 2)), 2) + 0.25 * pow((3. * (*(ffz + i)) - 4. * (*(ffz + i + 1)) + *(ffz + i + 2)), 2);
        tau = abs(sn2 - sn0);
        az0 = C03 * (1.0 + pow(tau / (ep + sn0), 2));
        az1 = C13 * (1.0 + pow(tau / (ep + sn1), 2));
        az2 = C23 * (1.0 + pow(tau / (ep + sn2), 2));
        W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
        q03 = 2. / 6. * (*(ffz + i - 2)) - 7. / 6. * (*(ffz + i - 1)) + 11. / 6. * (*(ffz + i));
        q13 = -1. / 6. * (*(ffz + i - 1)) + 5. / 6. * (*(ffz + i)) + 2. / 6. * (*(ffz + i + 1));
        q23 = 2. / 6. * (*(ffz + i)) + 5. / 6. * (*(ffz + i + 1)) - 1. / 6. * (*(ffz + i + 2));
        ss1 = W0 * q03 + W1 * q13 + W2 * q23;
        // weno_5-Z 负通量
        //  double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        //  double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        i = ii, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        sn0 = 13. / 12. * pow((*(fff + i + 1) - 2. * (*(fff + i + 2)) + *(fff + i + 3)), 2) + 0.25 * pow((3. * (*(fff + i + 1)) - 4. * (*(fff + i + 2)) + *(fff + i + 3)), 2);
        sn1 = 13. / 12. * pow((*(fff + i) - 2. * (*(fff + i + 1)) + *(fff + i + 2)), 2) + 0.25 * pow((*(fff + i) - *(fff + i + 2)), 2);
        sn2 = 13. / 12. * pow((*(fff + i - 1) - 2. * (*(fff + i)) + *(fff + i + 1)), 2) + 0.25 * pow((*(fff + i - 1) - 4. * (*(fff + i)) + 3. * (*(fff + i + 1))), 2);
        tau = abs(sn2 - sn0);
        az0 = C03 * (1.0 + pow(tau / (ep + sn0), 2));
        az1 = C13 * (1.0 + pow(tau / (ep + sn1), 2));
        az2 = C23 * (1.0 + pow(tau / (ep + sn2), 2));
        W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
        q03 = 11. / 6. * (*(fff + i + 1)) - 7. / 6. * (*(fff + i + 2)) + 2. / 6. * (*(fff + i + 3));
        q13 = 2. / 6. * (*(fff + i)) + 5. / 6. * (*(fff + i + 1)) - 1. / 6. * (*(fff + i + 2));
        q23 = -1. / 6. * (*(fff + i - 1)) + 5. / 6. * (*(fff + i)) + 2. / 6. * (*(fff + i + 1));
        ss2 = W0 * q03 + W1 * q13 + W2 * q23;
        yy = ss1 + ss2;
        */

        /*   weno7================================
        double S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33,
            ax0, ax1, ax2, ax3, am, q0, q1, q2, q3;
        double CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
               ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
               ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
               ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
               ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
               b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
               b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
               d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
        double e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
               e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
               e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
               e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-4; //   !! WENO-JS

        // ! 7th order WENO scheme
        // ! 1  阶导数      *(ffz + i - 3)
        S10 = ax11 * *(ffz + i - 3) + ax12 * *(ffz + i - 2) + ax13 * *(ffz + i - 1) + ax14 * *(ffz + i);
        S11 = ax21 * *(ffz + i - 2) - *(ffz + i - 1) + ax23 * *(ffz + i) + ax24 * *(ffz + i + 1);
        S12 = ax31 * *(ffz + i - 1) + ax32 * *(ffz + i) + *(ffz + i + 1) + ax34 * *(ffz + i + 2);
        S13 = ax41 * *(ffz + i) + ax42 * *(ffz + i + 1) + ax43 * *(ffz + i + 2) + ax44 * *(ffz + i + 3);
        //  ! 2 阶导数
        S20 = -*(ffz + i - 3) + b12 * *(ffz + i - 2) + b13 * *(ffz + i - 1) + b14 * *(ffz + i);
        S21 = *(ffz + i - 1) + b22 * *(ffz + i) + *(ffz + i + 1);
        S22 = *(ffz + i) + b22 * *(ffz + i + 1) + *(ffz + i + 2);
        S23 = b41 * *(ffz + i) + b42 * *(ffz + i + 1) + b43 * *(ffz + i + 2) - *(ffz + i + 3);
        // ! 3 阶导数
        S30 = -*(ffz + i - 3) + c12 * (*(ffz + i - 2) - *(ffz + i - 1)) + *(ffz + i);
        S31 = -*(ffz + i - 2) + c12 * (*(ffz + i - 1) - *(ffz + i)) + *(ffz + i + 1);
        S32 = -*(ffz + i - 1) + c12 * (*(ffz + i) - *(ffz + i + 1)) + *(ffz + i + 2);
        S33 = -*(ffz + i) + c12 * (*(ffz + i + 1) - *(ffz + i + 2)) + *(ffz + i + 3);

        S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
        S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
        S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
        S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

        // !-------WENO J-S----------------------
        ax0 = CC0 / ((ep + S0) * (ep + S0));
        ax1 = CC1 / ((ep + S1) * (ep + S1));
        ax2 = CC2 / ((ep + S2) * (ep + S2));
        ax3 = CC3 / ((ep + S3) * (ep + S3));
        // !-----------------------------------------------
        am = ax0 + ax1 + ax2 + ax3;

        // !  4阶差分格式的通量
        q0 = e11 * *(ffz + i - 3) + e12 * *(ffz + i - 2) + e13 * *(ffz + i - 1) + e14 * *(ffz + i);
        q1 = e21 * *(ffz + i - 2) + e22 * *(ffz + i - 1) + e23 * *(ffz + i) + e24 * *(ffz + i + 1);
        q2 = e31 * *(ffz + i - 1) + e32 * *(ffz + i) + e33 * *(ffz + i + 1) + e34 * *(ffz + i + 2);
        q3 = e41 * *(ffz + i) + e42 * *(ffz + i + 1) + e43 * *(ffz + i + 2) + e44 * *(ffz + i + 3);

        // !  由4个4阶差分格式组合成1个7阶差分格式
        // !     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3;
        ss1 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

        // -------------------------------
        i = ii; ////////////////////////
        // !      7th order WENO scheme
        // ! 1  阶导数   *(fff + i - 2)
        S10 = ax11 * *(fff + i + 4) + ax12 * *(fff + i + 3) + ax13 * *(fff + i + 2) + ax14 * *(fff + i + 1);
        S11 = ax21 * *(fff + i + 3) - *(fff + i + 2) + ax23 * *(fff + i + 1) + ax24 * *(fff + i);
        S12 = ax31 * *(fff + i + 2) + ax32 * *(fff + i + 1) + *(fff + i) + ax34 * *(fff + i - 1);
        S13 = ax41 * *(fff + i + 1) + ax42 * *(fff + i) + ax43 * *(fff + i - 1) + ax44 * *(fff + i - 2);
        // ! 2 阶导数
        S20 = -*(fff + i + 4) + b12 * *(fff + i + 3) + b13 * *(fff + i + 2) + b14 * *(fff + i + 1);
        S21 = *(fff + i + 2) + b22 * *(fff + i + 1) + *(fff + i);
        S22 = *(fff + i + 1) + b22 * *(fff + i) + *(fff + i - 1);
        S23 = b41 * *(fff + i + 1) + b42 * *(fff + i) + b43 * *(fff + i - 1) - *(fff + i - 2);
        // ! 3 阶导数
        S30 = -*(fff + i + 4) + c12 * (*(fff + i + 3) - *(fff + i + 2)) + *(fff + i + 1);
        S31 = -*(fff + i + 3) + c12 * (*(fff + i + 2) - *(fff + i + 1)) + *(fff + i);
        S32 = -*(fff + i + 2) + c12 * (*(fff + i + 1) - *(fff + i)) + *(fff + i - 1);
        S33 = -*(fff + i + 1) + c12 * (*(fff + i) - *(fff + i - 1)) + *(fff + i - 2);

        S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
        S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
        S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
        S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;

        ax0 = CC0 / ((ep + S0) * (ep + S0));
        ax1 = CC1 / ((ep + S1) * (ep + S1));
        ax2 = CC2 / ((ep + S2) * (ep + S2));
        ax3 = CC3 / ((ep + S3) * (ep + S3));

        // !-----------------------------------------------

        am = ax0 + ax1 + ax2 + ax3;

        // !  4阶差分格式的通量
        q0 = e11 * *(fff + i + 4) + e12 * *(fff + i + 3) + e13 * *(fff + i + 2) + e14 * *(fff + i + 1);
        q1 = e21 * *(fff + i + 3) + e22 * *(fff + i + 2) + e23 * *(fff + i + 1) + e24 * *(fff + i);
        q2 = e31 * *(fff + i + 2) + e32 * *(fff + i + 1) + e33 * *(fff + i) + e34 * *(fff + i - 1);
        q3 = e41 * *(fff + i + 1) + e42 * *(fff + i) + e43 * *(fff + i - 1) + e44 * *(fff + i - 2);

        // !  由4个4阶差分格式组合成1个7阶差分格式
        ss2 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

        yy = ss1 + ss2;
        */

        /*   weno7-Z================================
        double S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33,
            ax0, ax1, ax2, ax3, am, q0, q1, q2, q3, tau;
        double CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
               ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
               ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
               ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
               ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
               b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
               b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
               d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
        double e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
               e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
               e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
               e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-6; //   !! WENO-JS

        // ! 7th order WENO scheme
        // ! 1  阶导数      *(ffz + i - 3)
        S10 = ax11 * *(ffz + i - 3) + ax12 * *(ffz + i - 2) + ax13 * *(ffz + i - 1) + ax14 * *(ffz + i);
        S11 = ax21 * *(ffz + i - 2) - *(ffz + i - 1) + ax23 * *(ffz + i) + ax24 * *(ffz + i + 1);
        S12 = ax31 * *(ffz + i - 1) + ax32 * *(ffz + i) + *(ffz + i + 1) + ax34 * *(ffz + i + 2);
        S13 = ax41 * *(ffz + i) + ax42 * *(ffz + i + 1) + ax43 * *(ffz + i + 2) + ax44 * *(ffz + i + 3);
        //  ! 2 阶导数
        S20 = -*(ffz + i - 3) + b12 * *(ffz + i - 2) + b13 * *(ffz + i - 1) + b14 * *(ffz + i);
        S21 = *(ffz + i - 1) + b22 * *(ffz + i) + *(ffz + i + 1);
        S22 = *(ffz + i) + b22 * *(ffz + i + 1) + *(ffz + i + 2);
        S23 = b41 * *(ffz + i) + b42 * *(ffz + i + 1) + b43 * *(ffz + i + 2) - *(ffz + i + 3);
        // ! 3 阶导数
        S30 = -*(ffz + i - 3) + c12 * (*(ffz + i - 2) - *(ffz + i - 1)) + *(ffz + i);
        S31 = -*(ffz + i - 2) + c12 * (*(ffz + i - 1) - *(ffz + i)) + *(ffz + i + 1);
        S32 = -*(ffz + i - 1) + c12 * (*(ffz + i) - *(ffz + i + 1)) + *(ffz + i + 2);
        S33 = -*(ffz + i) + c12 * (*(ffz + i + 1) - *(ffz + i + 2)) + *(ffz + i + 3);

        S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
        S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
        S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
        S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
        // tau = abs(S0 - S1 - S2 + S3);
        tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
        // !-------WENO Z----------------------
        ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
        ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
        ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
        ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
        // !-------WENO J-S----------------------
        // ax0 = CC0 / ((ep + S0) * (ep + S0));
        // ax1 = CC1 / ((ep + S1) * (ep + S1));
        // ax2 = CC2 / ((ep + S2) * (ep + S2));
        // ax3 = CC3 / ((ep + S3) * (ep + S3));
        // !-----------------------------------------------
        am = ax0 + ax1 + ax2 + ax3;
        // !  4阶差分格式的通量
        q0 = e11 * *(ffz + i - 3) + e12 * *(ffz + i - 2) + e13 * *(ffz + i - 1) + e14 * *(ffz + i);
        q1 = e21 * *(ffz + i - 2) + e22 * *(ffz + i - 1) + e23 * *(ffz + i) + e24 * *(ffz + i + 1);
        q2 = e31 * *(ffz + i - 1) + e32 * *(ffz + i) + e33 * *(ffz + i + 1) + e34 * *(ffz + i + 2);
        q3 = e41 * *(ffz + i) + e42 * *(ffz + i + 1) + e43 * *(ffz + i + 2) + e44 * *(ffz + i + 3);
        // !  由4个4阶差分格式组合成1个7阶差分格式
        // !     hj(0)=W0*q0+W1*q1+W2*q2+W3*q3;
        ss1 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;

        // -------------------------------
        i = ii; ////////////////////////
        // !      7th order WENO scheme
        // ! 1  阶导数   *(fff + i - 2)
        S10 = ax11 * *(fff + i + 4) + ax12 * *(fff + i + 3) + ax13 * *(fff + i + 2) + ax14 * *(fff + i + 1);
        S11 = ax21 * *(fff + i + 3) - *(fff + i + 2) + ax23 * *(fff + i + 1) + ax24 * *(fff + i);
        S12 = ax31 * *(fff + i + 2) + ax32 * *(fff + i + 1) + *(fff + i) + ax34 * *(fff + i - 1);
        S13 = ax41 * *(fff + i + 1) + ax42 * *(fff + i) + ax43 * *(fff + i - 1) + ax44 * *(fff + i - 2);
        // ! 2 阶导数
        S20 = -*(fff + i + 4) + b12 * *(fff + i + 3) + b13 * *(fff + i + 2) + b14 * *(fff + i + 1);
        S21 = *(fff + i + 2) + b22 * *(fff + i + 1) + *(fff + i);
        S22 = *(fff + i + 1) + b22 * *(fff + i) + *(fff + i - 1);
        S23 = b41 * *(fff + i + 1) + b42 * *(fff + i) + b43 * *(fff + i - 1) - *(fff + i - 2);
        // ! 3 阶导数
        S30 = -*(fff + i + 4) + c12 * (*(fff + i + 3) - *(fff + i + 2)) + *(fff + i + 1);
        S31 = -*(fff + i + 3) + c12 * (*(fff + i + 2) - *(fff + i + 1)) + *(fff + i);
        S32 = -*(fff + i + 2) + c12 * (*(fff + i + 1) - *(fff + i)) + *(fff + i - 1);
        S33 = -*(fff + i + 1) + c12 * (*(fff + i) - *(fff + i - 1)) + *(fff + i - 2);

        S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
        S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
        S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
        S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
        // tau = abs(S0 - S1 - S2 + S3);
        tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
        // !-------WENO Z----------------------
        ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
        ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
        ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
        ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
        // !-------WENO J-S----------------------
        // ax0 = CC0 / ((ep + S0) * (ep + S0));
        // ax1 = CC1 / ((ep + S1) * (ep + S1));
        // ax2 = CC2 / ((ep + S2) * (ep + S2));
        // ax3 = CC3 / ((ep + S3) * (ep + S3));
        // !-----------------------------------------------
        am = ax0 + ax1 + ax2 + ax3;
        // !  4阶差分格式的通量
        q0 = e11 * *(fff + i + 4) + e12 * *(fff + i + 3) + e13 * *(fff + i + 2) + e14 * *(fff + i + 1);
        q1 = e21 * *(fff + i + 3) + e22 * *(fff + i + 2) + e23 * *(fff + i + 1) + e24 * *(fff + i);
        q2 = e31 * *(fff + i + 2) + e32 * *(fff + i + 1) + e33 * *(fff + i) + e34 * *(fff + i - 1);
        q3 = e41 * *(fff + i + 1) + e42 * *(fff + i) + e43 * *(fff + i - 1) + e44 * *(fff + i - 2);
        // !  由4个4阶差分格式组合成1个7阶差分格式
        ss2 = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;
        yy = ss1 + ss2;
        */

        return yy;
    }

    double dfdxcent(double *g, int i)
    {
        double y;
        // 1up
        // y = (*(g + i) - *(g + i - 1)) / dx;
        // y = (  *(g + i - 2) *0.5 - *(g + i - 1)*2.0 + *(g + i) *1.5) / dx;
        // 4center
        y = (*(g + i - 2) * 1. - *(g + i - 1) * 8. + *(g + i + 1) * 8. - *(g + i + 2) * 1.) / 12. / dx;

        //   6 center
        // y = (*(g + i - 3) / (-60.) + *(g + i - 2) * (9. / 60.) + *(g + i - 1) * (-45. / 60.) + *(g + i + 1) * (45. / 60.) + *(g + i + 2) * (-9. / 60.) + *(g + i + 3) / (60.)) / dx;

        // 5up
        // y = (*(g + i - 3) / (-30.) + *(g + i - 2) / 4. - *(g + i - 1) + *(g + i) / 3. + *(g + i + 1) / 2. - *(g + i + 2) / 20.) / dx;
        // 7up
        // y = (*(g + i - 4) * 3. - *(g + i - 3) * 28. + *(g + i - 2) * 126. - *(g + i - 1) * 420. + *(g + i) * 105. + *(g + i + 1) * 252. - *(g + i + 2) * 42. + *(g + i + 3) * 4.) / 420. / dx;
        // 8up
        // y = (-*(g + i - 3) * 5. + *(g + i - 2) * 60. - *(g + i - 1) * 420. - *(g + i) * 378 + *(g + i + 1) * 1050. - *(g + i + 2) * 420. + *(g + i + 3) * 140. - *(g + i + 4) * 30. + *(g + i + 5) * 3.) / 840. / dx;
        // y = -(-*(g + i + 3) * 5. + *(g + i + 2) * 60. - *(g + i + 1) * 420. - *(g + i) * 378 + *(g + i - 1) * 1050. - *(g + i - 2) * 420. + *(g + i - 3) * 140. - *(g + i - 4) * 30. + *(g + i - 5) * 3.) / 840. / dx;

        // 8up_center
        //  y = (*(g + i - 4) * 3. - *(g + i - 3) * 32. + *(g + i - 2) * 168. - *(g + i - 1) * 672. - *(g + i) * 0.0 + *(g + i + 1) * 672. - *(g + i + 2) * 168. + *(g + i + 3) * 32. - *(g + i + 4) * 3.) / 840. / dx;
        // 9up
        //  y = (*(g + i - 4) - *(g + i - 3) * 12. + *(g + i - 2) * 72. - *(g + i - 1) * 336. - *(g + i) * 100.8 + *(g + i + 1) * 504. - *(g + i + 2) * 168. + *(g + i + 3) * 48. - *(g + i + 4) * 9. + *(g + i + 5) * 0.8) / 504. / dx;
        // y = (*(g + i - 5) * (-4.) + *(g + i - 4) * 35. - *(g + i - 3) * 140. + *(g + i - 2) * 350. - *(g + i - 1) * 700. + *(g + i) * 3221. + *(g + i + 1) * 140. - *(g + i + 2) * 10.) / 420. / dx;
        // y = (*(g + i + 1) - *(g + i - 1)) / 2.0 / dx;
        // y = (45. * (*(g + i + 1) - *(g + i - 1)) - 21. * (*(g + i + 2) - *(g + i - 2)) + (*(g + i + 3) - *(g + i - 3))) / 60. / dx;
        // y = -(*(g + i + 4) * 3. - *(g + i + 3) * 28. + *(g + i + 2) * 126. - *(g + i + 1) * 420. + *(g + i) * 105. + *(g + i - 1) * 252. - *(g + i - 2) * 42. + *(g + i - 3) * 4.) / 420. / dx;

        // 3 up
        //  y = (*(g + i - 1) / (-3.) + *(g + i) / (-2.) + *(g + i + 1) + *(g + i + 2) / (-6.)) / dx;
        // 4 center
        // y = (*(g + i - 2) * 1. - *(g + i - 1) * 8. + *(g + i + 1) * 8. - *(g + i + 2) * 1.) / 12. / dx;
        //   6 center
        //  y = (*(g + i - 3) / (-60.) + *(g + i - 2) * (9. / 60.) + *(g + i - 1) * (-45. / 60.) + *(g + i + 1) * (45. / 60.) + *(g + i + 2) * (-9. / 60.) + *(g + i + 3) / (60.)) / dx;

        /*  //weno_5-js 负通量  ture https://doi.org/10.1007/s42967-019-0001-3
        double ep = 1.E-6, C03 = 1. / 10., C13 = 3. / 5., C23 = 3. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(g + i + 1) - 2. * (*(g + i + 2)) + *(g + i + 3)), 2) + 0.25 * pow((3. * (*(g + i + 1)) - 4. * (*(g + i + 2)) + *(g + i + 3)), 2);
            sn1 = 13. / 12. * pow((*(g + i) - 2. * (*(g + i + 1)) + *(g + i + 2)), 2) + 0.25 * pow((*(g + i) - *(g + i + 2)), 2);
            sn2 = 13. / 12. * pow((*(g + i - 1) - 2. * (*(g + i)) + *(g + i + 1)), 2) + 0.25 * pow((*(g + i - 1) - 4. * (*(g + i)) + 3. * (*(g + i + 1))), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 11. / 6. * (*(g + i + 1)) - 7. / 6. * (*(g + i + 2)) + 2. / 6. * (*(g + i + 3));
            q13 = 2. / 6. * (*(g + i)) + 5. / 6. * (*(g + i + 1)) - 1. / 6. * (*(g + i + 2));
            q23 = -1. / 6. * (*(g + i - 1)) + 5. / 6. * (*(g + i)) + 2. / 6. * (*(g + i + 1));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
        */

        /*  //weno_5-js
        double ep = 1.E-6, C03 = 3. / 10., C13 = 3. / 5., C23 = 1. / 10.;
        double sn0, sn1, sn2, az0, az1, az2, W0, W1, W2, q03, q13, q23, vz[2];
        for (int j = 0; j < 2; j++)
        {
            sn0 = 13. / 12. * pow((*(g + i) - 2. * (*(g + i + 1)) + *(g + i + 2)), 2) + 0.25 * pow((3. * (*(g + i)) - 4. * (*(g + i + 1)) + *(g + i + 2)), 2);
            sn1 = 13. / 12. * pow((*(g + i - 1) - 2. * (*(g + i)) + *(g + i + 1)), 2) + 0.25 * pow((*(g + i - 1) - *(g + i + 1)), 2);
            sn2 = 13. / 12. * pow((*(g + i - 2) - 2. * (*(g + i - 1)) + *(g + i)), 2) + 0.25 * pow((*(g + i - 2) - 4. * (*(g + i - 1)) + 3. * (*(g + i))), 2);
            az0 = C03 / (pow((ep + sn0), 2)), az1 = C13 / (pow((ep + sn1), 2)), az2 = C23 / (pow((ep + sn2), 2));
            W0 = az0 / (az0 + az1 + az2), W1 = az1 / (az0 + az1 + az2), W2 = az2 / (az0 + az1 + az2);
            q03 = 1. / 3. * (*(g + i)) + 5. / 6. * (*(g + i + 1)) - 1. / 6. * (*(g + i + 2));
            q13 = -1. / 6. * (*(g + i - 1)) + 5. / 6. * (*(g + i)) + 1. / 3. * (*(g + i + 1));
            q23 = 1. / 3. * (*(g + i - 2)) - 7. / 6. * (*(g + i - 1)) + 11. / 6. * (*(g + i));
            vz[j] = W0 * q03 + W1 * q13 + W2 * q23;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
        */

        /*   weno7-Z---负通量------------------------------
      double S0, S1, S2, S3, S10, S11, S12, S13, S20, S21, S22, S23, S30, S31, S32, S33,
          ax0, ax1, ax2, ax3, am, q0, q1, q2, q3, vz[2], tau;
      double CC0 = 1. / 35., CC1 = 12. / 35., CC2 = 18. / 35., CC3 = 4. / 35.,
             ax11 = -2. / 6., ax12 = 9. / 6., ax13 = -18. / 6., ax14 = 11. / 6.,
             ax21 = 1. / 6., ax23 = 3. / 6., ax24 = 2. / 6.,
             ax31 = -2. / 6., ax32 = -3. / 6., ax34 = -1. / 6.,
             ax41 = -11. / 6., ax42 = 18. / 6., ax43 = -9. / 6., ax44 = 2. / 6.,
             b12 = 4., b13 = -5., b14 = 2., b22 = -2.,
             b41 = 2., b42 = -5., b43 = 4., c12 = 3.,
             d12 = 13. / 12., d13 = 1043. / 960., d14 = 1. / 12.;
      double e11 = -3. / 12., e12 = 13. / 12., e13 = -23. / 12., e14 = 25. / 12.,
             e21 = 1. / 12., e22 = -5. / 12., e23 = 13. / 12., e24 = 3. / 12.,
             e31 = -1. / 12., e32 = 7. / 12., e33 = 7. / 12., e34 = -1. / 12.,
             e41 = 3. / 12., e42 = 13. / 12., e43 = -5. / 12., e44 = 1. / 12., ep = 1.e-6; //   !! WENO-JS

      for (int j = 0; j < 2; j++)
      {
          // !      7th order WENO scheme
          // ! 1  阶导数   *(g + i - 2)
          S10 = ax11 * *(g + i + 4) + ax12 * *(g + i + 3) + ax13 * *(g + i + 2) + ax14 * *(g + i + 1);
          S11 = ax21 * *(g + i + 3) - *(g + i + 2) + ax23 * *(g + i + 1) + ax24 * *(g + i);
          S12 = ax31 * *(g + i + 2) + ax32 * *(g + i + 1) + *(g + i) + ax34 * *(g + i - 1);
          S13 = ax41 * *(g + i + 1) + ax42 * *(g + i) + ax43 * *(g + i - 1) + ax44 * *(g + i - 2);
          // ! 2 阶导数
          S20 = -*(g + i + 4) + b12 * *(g + i + 3) + b13 * *(g + i + 2) + b14 * *(g + i + 1);
          S21 = *(g + i + 2) + b22 * *(g + i + 1) + *(g + i);
          S22 = *(g + i + 1) + b22 * *(g + i) + *(g + i - 1);
          S23 = b41 * *(g + i + 1) + b42 * *(g + i) + b43 * *(g + i - 1) - *(g + i - 2);
          // ! 3 阶导数
          S30 = -*(g + i + 4) + c12 * (*(g + i + 3) - *(g + i + 2)) + *(g + i + 1);
          S31 = -*(g + i + 3) + c12 * (*(g + i + 2) - *(g + i + 1)) + *(g + i);
          S32 = -*(g + i + 2) + c12 * (*(g + i + 1) - *(g + i)) + *(g + i - 1);
          S33 = -*(g + i + 1) + c12 * (*(g + i) - *(g + i - 1)) + *(g + i - 2);

          S0 = S10 * S10 + d12 * S20 * S20 + d13 * S30 * S30 + d14 * S10 * S30;
          S1 = S11 * S11 + d12 * S21 * S21 + d13 * S31 * S31 + d14 * S11 * S31;
          S2 = S12 * S12 + d12 * S22 * S22 + d13 * S32 * S32 + d14 * S12 * S32;
          S3 = S13 * S13 + d12 * S23 * S23 + d13 * S33 * S33 + d14 * S13 * S33;
          tau = abs(S0 + 3. * S1 - 3. * S2 - S3);
          // !-------WENO Z----------------------
          ax0 = CC0 * (1.0 + pow((tau / (S0 + ep)), 2));
          ax1 = CC1 * (1.0 + pow((tau / (S1 + ep)), 2));
          ax2 = CC2 * (1.0 + pow((tau / (S2 + ep)), 2));
          ax3 = CC3 * (1.0 + pow((tau / (S3 + ep)), 2));
          // ax0 = CC0 / ((ep + S0) * (ep + S0));
          // ax1 = CC1 / ((ep + S1) * (ep + S1));
          // ax2 = CC2 / ((ep + S2) * (ep + S2));
          // ax3 = CC3 / ((ep + S3) * (ep + S3));
          // !-----------------------------------------------
          am = ax0 + ax1 + ax2 + ax3;
          // !  4阶差分格式的通量
          q0 = e11 * *(g + i + 4) + e12 * *(g + i + 3) + e13 * *(g + i + 2) + e14 * *(g + i + 1);
          q1 = e21 * *(g + i + 3) + e22 * *(g + i + 2) + e23 * *(g + i + 1) + e24 * *(g + i);
          q2 = e31 * *(g + i + 2) + e32 * *(g + i + 1) + e33 * *(g + i) + e34 * *(g + i - 1);
          q3 = e41 * *(g + i + 1) + e42 * *(g + i) + e43 * *(g + i - 1) + e44 * *(g + i - 2);
          // !  由4个4阶差分格式组合成1个7阶差分格式
          vz[j] = (ax0 * q0 + ax1 * q1 + ax2 * q2 + ax3 * q3) / am;
          i--;
      }
      y = (vz[0] - vz[1]) / dx;
      */

        /*  //weno_9
        double uim4, uim3, uim2, uim1, ui, uip1, uip2, uip3, uip4, vz[2];
        double uz1, uz2, uz3, uz4, uz5;

        double beta0, beta1, beta2, beta3, beta4, omt0, omt1, omt2, omt3, omt4, omts;

        const double g0 = 1.0 / 126.0, g1 = 10. / 63., g2 = 10. / 21., g3 = 20. / 63., g4 = 5. / 126., eps = 1e-6;

        for (int j = 0; j < 2; j++)
        {
            uim4 = *(g + i - 4), uim3 = *(g + i - 3), uim2 = *(g + i - 2), uim1 = *(g + i - 1);
            ui = *(g + i), uip1 = *(g + i + 1), uip2 = *(g + i + 2), uip3 = *(g + i + 3), uip4 = *(g + i + 4);

            beta0 = uim4 * (22658. * uim4 - 208501. * uim3 + 364863. * uim2 - 288007. * uim1 + 86329. * ui) +
                    uim3 * (482963. * uim3 - 1704396. * uim2 + 1358458. * uim1 - 411487. * ui) +
                    uim2 * (1521393. * uim2 - 2462076. * uim1 + 758823. * ui) +
                    uim1 * (1020563. * uim1 - 649501. * ui) + 107918. * pow(ui, 2);
            beta1 = uim3 * (6908. * uim3 - 60871. * uim2 + 99213. * uim1 - 70237. * ui + 18079. * uip1) +
                    uim2 * (138563. * uim2 - 464976. * uim1 + 337018. * ui - 88297. * uip1) +
                    uim1 * (406293. * uim1 - 611976 * ui + 165153. * uip1) +
                    ui * (242723. * ui - 140251. * uip1) + 22658. * pow(uip1, 2);
            beta2 = uim2 * (6908. * uim2 - 51001. * uim1 + 67923. * ui - 38947. * uip1 + 8209. * uip2) +
                    uim1 * (104963. * uim1 - 299076. * ui + 179098. * uip1 - 38947. * uip2) +
                    ui * (231153. * ui - 299076. * uip1 + 67923. * uip2) +
                    uip1 * (104963. * uip1 - 51001. * uip2) + 6908. * pow(uip2, 2);
            beta3 = uim1 * (22658. * uim1 - 140251. * ui + 165153. * uip1 - 88297. * uip2 + 18079. * uip3) +
                    ui * (242723. * ui - 611976. * uip1 + 337018. * uip2 - 70237. * uip3) +
                    uip1 * (406293. * uip1 - 464976. * uip2 + 99213. * uip3) +
                    uip2 * (138563. * uip2 - 60871. * uip3) + 6908. * pow(uip3, 2);
            beta4 = ui * (107918. * ui - 649501. * uip1 + 758823. * uip2 - 411487. * uip3 + 86329. * uip4) +
                    uip1 * (1020563. * uip1 - 2462076. * uip2 + 1358458. * uip3 - 288007 * uip4) +
                    uip2 * (1521393. * uip2 - 1704396. * uip3 + 364863. * uip4) +
                    uip3 * (482963. * uip3 - 208501. * uip4) + 22658. * pow(uip4, 2);

            uz1 = (1. / 5.) * uim4 + (-21. / 20.) * uim3 + (137. / 60.) * uim2 + (-163. / 60.) * uim1 + (137. / 60.) * ui;
            uz2 = (-1. / 20.) * uim3 + (17. / 60.) * uim2 + (-43. / 60.) * uim1 + (77. / 60.) * ui + (1. / 5.) * uip1;
            uz3 = (1. / 30.) * uim2 + (-13. / 60.) * uim1 + (47. / 60.) * ui + (9. / 20.) * uip1 + (-1. / 20.) * uip2;
            uz4 = (-1. / 20.) * uim1 + (9. / 20.) * ui + (47. / 60.) * uip1 + (-13. / 60.) * uip2 + (1. / 30.) * uip3;
            uz5 = (1. / 5.) * ui + (77. / 60.) * uip1 + (-43. / 60.) * uip2 + (17. / 60.) * uip3 + (-1. / 20.) * uip4;

            omt0 = g0 * pow(eps + beta0, -2), omt1 = g1 * pow(eps + beta1, -2), omt2 = g2 * pow(eps + beta2, -2), omt3 = g3 * pow(eps + beta3, -2);
            omt4 = g4 * pow(eps + beta4, -2), omts = omt0 + omt1 + omt2 + omt3 + omt4;
            vz[j] = (omt0 * uz1 + omt1 * uz2 + omt2 * uz3 + omt3 * uz4 + omt4 * uz5) / omts;
            i--;
        }
        y = (vz[0] - vz[1]) / dx;
          */

        return y;
    }

    double dfdx(double *g, int i)
    {
        double y;
        // 1up
        y = (*(g + i + 1) - *(g + i)) / dx;

        return y;
    }

    double dfdxx(double *g, int i)
    {
        double y;
        // 2ord
        y = (*(g + i + 1) - *(g + i) * 2.0 + *(g + i - 1)) / dx / dx;
        return y;
    }

    //---------------------------------------------------------
    void compt_CFL()
    {
        int i, j;
        double *umax;
        i = n1 - nb * 2;
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            uc_max[j] = abs(u[j]) + c[j];
        }
        umax = max_element(uc_max + nb, uc_max + n1 - nb);

        dt = CFL * dx / (*umax);
        // dt = CFL * dx;
        // dt = CFL * dx / 2.75; // vacuum
        // dt = CFL * dx / 5.0; // Shu-Osher
        // dt = CFL * dx / 2.; // Sod
        // cout << "; \n   umax: " << *umax << ";    dt: " << dt;
        // dt=0.02/1.0;        dt = CFL;
    }

    void maxTV()
    {
        // double TV = 0.0;
        // for (int i = 0; i < n1 - 1; i++) //间断
        // {
        //     // TV = abs(a[i + 1] - a[i]) + TV;
        //     TV = abs(p[i + 1] - p[i]) + TV;

        // }
        // // TV = abs(TV - 0.875);
        // TV = abs(TV - 0.9);
        // // cout << " \n TV: " << TV;
        // if (ii1 > -10) // tt
        // {
        //     if (TV > TVmax)
        //     {
        //         TVmax = TV;
        //     }
        // }
        TVmax = 1.0;
        for (int i = 0; i < n1 - 1; i++) // 间断
        {
            if (p[i + 1] > 0.0 & a[i + 1] > 0.0)
            {
                if (p[i + 1] < TVmax)
                    TVmax = p[i + 1];
                if (a[i + 1] < TVmax)
                    TVmax = a[i + 1];
            }
            else
            {
                TVmax = 11111.0;
                break;
            }
        }
    }

    void TVabs(int m)
    {
        double TV = 0.0;
        TV = TVmax;

        // if (TV < 2.0E-16)
        // {
        //     TV = 2.0E-16;
        // }
        cout << "\n  CFL: " << CFL << ",  TV: " << TV;
        string Title = "TV_CFL.plt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=CFL,TV " << endl;
            ofs.precision(5), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
            ofs.precision(5), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
    }
};
////-----------------------------------------
void comput_tt(cfda &);
int comput_main(int, double, double, int);
int main()
{
    int n, oder_n = 2, begin_cell = 200; // 200
    double t = 0.2, CFL = 1.2;           // 计算时间总长，CFL数 0.6  1.1  0.6
    n = begin_cell + 1;
    for (int m = 1; m < oder_n; m++)
    {
        comput_main(n, t, CFL, m);
        // CFL = CFL + 0.0001;
    }
    cout << " \n All comput time: ";
    printf("%d ms", clock()); // 输出运行所费时间，单位毫秒ms
    // system("pause");
    return 2;
}

int comput_main(int n, double t, double CFL, int m)
{
    //
    int nb = 8; // 边界网格结点数（有限差分）
    int n1 = n + 8 * 2;
    double main_df[3][n + 8 * 2], u_excat[3][n + 8 * 2];
    double Ltt_n0[3][n + 8 * 2], x[n + 8 * 2], Ltt_n1[3][n + 8 * 2], Lttx_n1[3][n + 8 * 2]; // 声明变长数组

    double L_n0[3][n + 8 * 2], L_n1[3][n + 8 * 2], L_n2[3][n + 8 * 2];       // f=L
    double Lt_n0[3][n + 8 * 2], Lt_n1[3][n + 8 * 2], Lt_n2[3][n + 8 * 2];    // L_t
    double Ltx_n0[3][n + 8 * 2], Ltx_n1[3][n + 8 * 2], Ltx_n2[3][n + 8 * 2]; // L_t

    double u_n0[3][n + 8 * 2], u_n1[3][n + 8 * 2], u_n2[3][n + 8 * 2], u_n3[3][n + 8 * 2]; // 推进步中间时刻
    double u_n4[3][n + 8 * 2], u_n5[3][n + 8 * 2], u_n6[3][n + 8 * 2], u_n7[3][n + 8 * 2]; // 推进步中间时刻
    double f_n0[3][n + 8 * 2], f_n1[3][n + 8 * 2], f_n2[3][n + 8 * 2], f_n3[3][n + 8 * 2]; // 推进步中间时刻
    double f_n4[3][n + 8 * 2], f_n5[3][n + 8 * 2], f_n6[3][n + 8 * 2], f_n7[3][n + 8 * 2]; // 推进步中间时刻

    double u_nn[3][n + 8 * 2], f_nn[3][n + 8 * 2], Lttx_n0[3][n + 8 * 2]; // 下一时刻（n+1)

    double TL_n0[3][n + 8 * 2], TL_n1[3][n + 8 * 2], TL_n2[3][n + 8 * 2];    // f=L
    double TLt_n0[3][n + 8 * 2], TLt_n1[3][n + 8 * 2], TLt_n2[3][n + 8 * 2]; // TL_t

    double Tu1[3][n + 8 * 2], Tu2[3][n + 8 * 2], uc_max[n + 8 * 2];
    double fz[3][n + 8 * 2], ff[3][n + 8 * 2];                     // 守恒通量
    double u[n + 8 * 2], c[n + 8 * 2], p[n + 8 * 2], a[n + 8 * 2]; // 原始变量  u速度，c声速，p压强，a密度
    double fluxz[3][n + 8 * 2], fluxf[3][n + 8 * 2], flux[3][n + 8 * 2];
    class cfda cff;
    cff.kt = 1;
    cff.A1 = 1. / 3.0;
    cff.A2 = 2. / 3.0;
    cff.t = t;          // 1.0 / cff.PI; //计算总时间
    cff.nx_begin = 0.0; // 网格左端点
    cff.nx_end = 1.0;   // 网格右端点
    cff.CFL = CFL, cff.TVmax = 0.0;

    cff.n = n, cff.nb = nb, cff.n1 = n1, cff.m = m;
    cff.Ltt_n0 = Ltt_n0[0], cff.x = x, cff.Ltt_n1 = Ltt_n1[0], cff.Lttx_n1 = Lttx_n1[0]; // 声明变长数组
    cff.L_n0 = L_n0[0], cff.L_n1 = L_n1[0], cff.L_n2 = L_n2[0];                          // f=L
    cff.Lt_n0 = Lt_n0[0], cff.Lt_n1 = Lt_n1[0], cff.Lt_n2 = Lt_n2[0];                    // L_t
    cff.Ltx_n0 = Ltx_n0[0], cff.Ltx_n1 = Ltx_n1[0], cff.Ltx_n2 = Ltx_n2[0];              // L_t

    cff.u_n0 = u_n0[0], cff.u_n1 = u_n1[0], cff.u_n2 = u_n2[0], cff.u_n3 = u_n3[0]; // 推进步中间时刻
    cff.u_n4 = u_n4[0], cff.u_n5 = u_n5[0], cff.u_n6 = u_n6[0], cff.u_n7 = u_n7[0]; // 推进步中间时刻
    cff.f_n0 = f_n0[0], cff.f_n1 = f_n1[0], cff.f_n2 = f_n2[0], cff.f_n3 = f_n3[0]; // 推进步中间时刻
    cff.f_n4 = f_n4[0], cff.f_n5 = f_n5[0], cff.f_n6 = f_n6[0], cff.f_n7 = f_n7[0]; // 推进步中间时刻
    cff.u_nn = u_nn[0], cff.f_nn = f_nn[0], cff.Lttx_n0 = Lttx_n0[0];               // 下一时刻（n+1)

    cff.TL_n0 = TL_n0[0], cff.TL_n1 = TL_n1[0], cff.TL_n2 = TL_n2[0];
    cff.TLt_n0 = TLt_n0[0], cff.TLt_n1 = TLt_n1[0], cff.TLt_n2 = TLt_n2[0];

    cff.Tu1 = Tu1[0], cff.Tu2 = Tu2[0], cff.uc_max = uc_max;
    cff.df = main_df[0], cff.u_excat = u_excat[0],
    cff.u = u, cff.c = c, cff.p = p, cff.a = a;                     // 原始变量  u速度，c声速，p压强，a密度
    cff.fz = fz[0], cff.ff = ff[0];                                 // 守恒通量
    cff.flux = flux[0], cff.fluxz = fluxz[0], cff.fluxf = fluxf[0]; // 守恒通量

    cff.intc();
    comput_tt(cff);
    cff.Store_obj(cff.put_obj, cff.u, cff.a, cff.p); // 存储数据
    cff.Write_obj(-1);
    cff.TVabs(m);
    return 2;
}

////-----------------------------------------

void comput_tt(cfda &cff1)
{
    int n_all, i3 = 0, tt_flag = 0;
    cff1.tt = 0, n_all = cff1.n1;
    for (int i1 = 0; i1 < 4000000; i1++) // 200 50   15
    {
        // cff1.Store_obj(cff1.put_obj, cff1.u, cff1.a, cff1.p); //存储数据
        // cff1.Write_obj(-1);
        if (tt_flag > 1)
        {
            cout.precision(18); // 精度为18，正常为6
            cout << " \n comput time: " << cff1.tt << " comput time n:   " << i1;
            break;
        }
        cff1.compt_CFL();
        if (cff1.t < cff1.tt + cff1.dt + 1E-10)
        {
            tt_flag = 3;
            cff1.dt = cff1.t - cff1.tt; //
        }
        if (i1 < 50000000) // 4  5
        {
 
            // cff1.compt_ThD13();
            // cff1.compt_ThD24();
            // cff1.compt_ThD25();

            // cff1.RK33_compt();
            // cff1.SSPRK54_compt_Shu();
            // cff1.LSRK65_compt2();

            // cff1.compt23_LLttnn();
            // cff1.compt24_LLttnn(); // cff1.compt24_LLttn12();
            // cff1.compt35_LLttnn();

            // cff1.compt_ThD13_Char();
            // cff1.compt_ThD23_Char();
            cff1.compt_ThD24_Char();
            // cff1.compt_ThD25_Char();
            // cff1.compt_ThD36_Char();

            // cff1.RK33_compt_Char();
            // cff1.SSPRK54_compt_Shu_Char();
            // cff1.LSRK65_compt2_Char();

            // cff1.compt23_LLttnn_Char();
            // cff1.compt24_LLttnn_Char();
 
        }
 
        cff1.tt = cff1.tt + cff1.dt;
        //___________________________________________________________________
        cff1.carr(cff1.u_n0, cff1.u_nn, n_all * 3);
        cff1.uu_to_pc(cff1.u_nn);
        // cout << " \n comput time: " << cff1.tt << " comput time n:   " << i1;
        // cff1.border();
        cff1.maxTV();
        // if (cff1.TVmax > 11000.0)
        //     break;
        // if (i3 * cff1.t / cff1.onet_out < cff1.tt + 0.00000001)
        // {
        //     // cff1.Store_obj(cff1.put_one_n[i3], cff1.u_nn); //存储数据
        //     i3++;
        // }
    }
}
