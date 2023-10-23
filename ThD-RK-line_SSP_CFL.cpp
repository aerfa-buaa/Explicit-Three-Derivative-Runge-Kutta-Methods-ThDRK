#include <math.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <string>
#include <fstream>
#include <iterator>
#include <valarray>
using namespace std;
// example:line 3-4-L-W   增加CFL数输出误差   Line Excample-1  精简版本只有1up,2cneter  U_t+U_x=0
// g++ -O2 -std=c++11 -Wall "-Wl,--stack=268435456" Taylor3-RK-line_SSP_CFL.cpp -o a.exe
class cfda
{

private:
    int w;

public:
    double PI = 3.141592653589793;
    // double PI = 3.1415926535897932385;

    int n;                                         // 网格结点数（有限差分）
    int nb;                                        // 边界网格结点数（有限差分）
    double t, tt, CFL;                             // 计算总时 间t,实际计算时间tt;
    int nt, kt, ij, ijj, onet_out = 2;             // 计算时间步长;
    double nx_begin, nx_end;                       // 网格长度;
    int n1;                                        // 总的网格结点数（有限差分）
    double *Ltt_n0, *x;                            // 声明变长数组
    double *L_n0, *L_n1, *L_n2, *L_n3, *L_n4;      // f=L
    double *Lt_n0, *Lt_n1, *Lt_n2, *Lt_n3, *Lt_n4; // L_t

    double *TL_n0, *TL_n1, *TL_n2, *TL_n3, *TL_n4;      // Tf=L
    double *TLt_n0, *TLt_n1, *TLt_n2, *TLt_n3, *TLt_n4; // TL_t

    double *u_n0, *u_n1, *u_n2, *u_n3, *u_n4, *u_n5, *u_n6, *u_n7, *u_n8, *u_n9; // 推进步中间时刻
    double *f_n0, *f_n1, *f_n2, *f_n3, *f_n4, *f_n5, *f_n6, *f_n7, *f_n8, *f_n9; // 推进步中间时刻
    double *u_nn, *f_nn, *Lttx_n0;                                               // 下一时刻（n+1)
    double *Tu1, *Tu2;
    // double  u_excat[513 + 6 * 2];
    double *df, *u_excat; // du[513 + 6 * 2],df[513 + 6 * 2],
    ////-----------------------------------------
    double dx;
    double dt;
    double *dd; // ss
    double A1, A2, A3, A4, A5, B0, B1, B2, B3, bb1, bb2, bb3, bbb1, bbb2, bbb3, C0, C1, C2, C3, D0, D1, D2, D3, maxuu, TVmax, k;
    double a21, a32, a31, a41, a42, a43, aa12, aa21, aa31, aa32, aa41, aa42, aa43, aaa31, aaa32, aaa21, w1, w2, v1, v2, v3, v4, vv1, vv2, vv3, vv4, ww1, ww2;
    double put_obj[2][1300];   //
    double put_one_n[2][1300]; // onet_out = 10;
    ////-----------------------------------------
    void intc()
    {
        dx = (nx_end - nx_begin) / (1.0 * (n - 1));

        // for (int i = 0; i < n1; i++)
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     u_n0[i] = intc_fun(x[i]);
        //     df[i] = 1;
        // }

        for (int i = 0; i < n1; i++) // 间断
        {
            x[i] = nx_begin + (i - nb) * dx;
            u_n0[i] = 0.0;
            if (x[i] >= -0.5 && x[i] <= 0.5)
            {
                u_n0[i] = 1.0;
            }
            // df[i] = 1;
        }
        // Write_obj(0);
    }

    double intc_fun(double xxf)
    {
        double yyf;
        // yyf = (sin(2 * PI * xxf) / 2.0 + 0.5) * 1.0; //  / PIbergure
        // yyf = sin(PI * xxf) / 2.0 + 0.25;//bergure
        yyf = sin(PI * xxf); // line
        return yyf;
    }

    // //-----------------------------------------
    void border()
    {
        borderfun(u_n0), borderfun(u_n1), borderfun(u_n2), borderfun(u_n3), borderfun(u_n4), borderfun(u_n5);
        borderfun(u_n6), borderfun(u_n7), borderfun(u_n8), borderfun(u_n9), borderfun(u_nn);
        borderfun(f_n0), borderfun(f_n1), borderfun(f_n2), borderfun(f_n3), borderfun(f_n4), borderfun(f_n5);
        borderfun(f_n6), borderfun(f_n7), borderfun(f_n8), borderfun(f_n9), borderfun(f_nn);
        borderfun(L_n0), borderfun(L_n1), borderfun(L_n2), borderfun(L_n3), borderfun(L_n4);
        borderfun(Lt_n0), borderfun(Lt_n1), borderfun(Lt_n2), borderfun(Lt_n3), borderfun(Lt_n4);
        borderfun(Ltt_n0), borderfun(Lttx_n0);
    }

    void borderfun(double *ffbc)
    {

        for (int i = 0; i < nb; i++)
        {

            *(ffbc + i) = *(ffbc + n1 - 2 * nb + i - 1);  // 周期边条
            *(ffbc + n1 - nb + i) = *(ffbc + nb + i + 1); // 周期边条

            // *(ffbc + i) = 0.;
            // *(ffbc + n1 - nb + i) = 0.;

            // f_n0[i] = f_n0[nb];                   //恒定边条
            // f_n0[n1 - 1 - i] = f_n0[n1 - 1 - nb]; //恒定边条

            // *(ffbc + i) = *(ffbc + nb);                   //恒定边条
            // *(ffbc + n1 - 1 - i) = *(ffbc + n1 - 1 - nb); //恒定边条

            // *(ffbc + i) = 0.0;                   //恒定边条
            // *(ffbc + n1 - 1 - i) = 0.0; //恒定边条
            // ffbc++;
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
    void compt_ThD_13_line_SSP()
    {
        int i, j;
        a21 = 0.594223212099088, v2 = 0.306027487008159;
        ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.0;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            u_nn[j] = u_n0[j] + dt * L_n0[j] + dt * dt * Lt_n0[j] / 2.0 + dt * dt * dt * TLt_n0[j] / 6.0;
        }
    }

    void compt_ThD_23_line_SSP()
    {
        int i, j;
        a21 = 0.500;
        B2 = 0.5;
        bb2 = 0.125;
        bbb1 = 0.020833333333333;
        bbb2 = 1. / 6. * (1. - 3. * a21 * (a21 * B2 + 2. * bb2) - 6. * bbb1);
        B1 = 1. - B2;
        bb1 = 1. / 2. - a21 * B2 - bb2;

        // a21 = 0.614823;
        // B2 = 0.402390;
        // bb2 = 0.082793;
        // bbb1 = 0.0227428;
        // bbb2 = 0.0169676; // 1 / 6 * (1 - 3 * a21 * (a21 * B2 + 2 * bb2) - 6 * bbb1);
        // B1 = 1. - B2, bb1 = 1. / 2. - a21 * B2 - bb2;

        A1 = a21, A2 = a21 * a21 / 2.0, A3 = a21 * a21 * a21 / 6.0;

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            f_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A2 * dt * dt * Lt_n0[j] + A3 * dt * dt * dt * TLt_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            TLt_n1[j] = dfdxxx(f_n1, j);
            u_nn[j] = u_n0[j] + dt * (B1 * L_n0[j] + B2 * L_n1[j]) +
                      dt * dt * (bb1 * Lt_n0[j] + bb2 * Lt_n1[j]) +
                      dt * dt * dt * (bbb1 * TLt_n0[j] + bbb2 * TLt_n1[j]);
        }
    }

    void compt_ThD_24_line_SSP()
    {
        int i, j;
        a21 = 0.614823, A1 = a21, A2 = a21 * a21 / 2.0, A3 = a21 * a21 * a21 / 6.0;
        // B2 = 0.40031, bb2 = 0.090365;
        B2 = 0.402390, bb2 = 0.082793;
        B1 = 1. - B2;
        bb1 = 1. / 2. - a21 * B2 - bb2;
        bbb1 = 1. / 24. * (4. - 1. / a21 - 8. * a21 * a21 * B2 - 12. * a21 * bb2);
        bbb2 = 1. / (24 * a21) - 1. / 6. * a21 * (a21 * B2 + 3. * bb2);

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            f_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A2 * dt * dt * Lt_n0[j] + A3 * dt * dt * dt * TLt_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            TLt_n1[j] = dfdxxx(f_n1, j);
            u_nn[j] = u_n0[j] + dt * (B1 * L_n0[j] + B2 * L_n1[j]) +
                      dt * dt * (bb1 * Lt_n0[j] + bb2 * Lt_n1[j]) +
                      dt * dt * dt * (bbb1 * TLt_n0[j] + bbb2 * TLt_n1[j]);
        }
    }

    void compt_ThD_25_line_SSP()
    {
        int i, j;
        a21 = 0.60, A1 = a21, A2 = a21 * a21 / 2.0, A3 = a21 * a21 * a21 / 6.0;

        B1 = 1., B2 = 0.;
        bb1 = (2. - 5. * a21 + 10. * a21 * a21 * a21) / (20. * a21 * a21 * a21), bb2 = (-2. + 5. * a21) / (20. * a21 * a21 * a21);
        bbb1 = (3. + 10. * (-1. + a21) * a21) / (60. * a21 * a21), bbb2 = (3. - 5. * a21) / (60. * a21 * a21);

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            f_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A2 * dt * dt * Lt_n0[j] + A3 * dt * dt * dt * TLt_n0[j];
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            TLt_n1[j] = dfdxxx(f_n1, j);
            u_nn[j] = u_n0[j] + dt * (B1 * L_n0[j] + B2 * L_n1[j]) +
                      dt * dt * (bb1 * Lt_n0[j] + bb2 * Lt_n1[j]) +
                      dt * dt * dt * (bbb1 * TLt_n0[j] + bbb2 * TLt_n1[j]);
        }
    }

    void compt_ThD_36_line_SSP()
    {
        int i, j;

        // double xc[14] = {111., 0.869270422229783, 0.451450450187708, 0, 0.0879532579603533,
        //                  0.0139504965269889, 0.00320807794004533, 0, 0.197279956148903, 0.0392094809402299, 0.263510562910867,
        //                  0.0136210623590151, 0.0, 0.};
        // a21 = xc[1], a31 = xc[2], a32 = xc[3], aa31 = xc[4], aa32 = xc[5], aaa31 = xc[6], aaa32 = xc[7];
        // B1 = 1.0, B2 = 0., B3 = 0., bb1 = xc[8], bb2 = xc[9], bb3 = xc[10], bbb1 = xc[11], bbb2 = xc[12], bbb3 = xc[13];

        a21 = 0.869270422229783, A1 = a21, A2 = a21 * a21 / 2.0, A3 = a21 * a21 * a21 / 6.0;
        a21 = 1.01640357753173;
        B1 = 1., B2 = 0., B3 = 0.;
        aaa32 = 0., a32 = 0., bbb2 = 0., bbb3 = 0;
        bb1 = (-2. + 5. * pow(a21, 2.) * (16. + a21 * (-42. + 29. * a21))) / (60. * pow(2. - 3. * a21, 2) * pow(a21, 2.));
        bb2 = 1. / (60. * pow(a21, 2.) * (2. + a21 * (-6. + 5. * a21)));
        bb3 = pow(3. - 5. * a21, 4.) / (60. * pow(2. - 3. * a21, 2.) * (2. + a21 * (-6. + 5. * a21)));
        bbb1 = (1. + 5. * (-1. + a21) * a21) / (60. * a21 * (-2. + 3. * a21));
        a31 = (2. - 3. * a21) / (3. - 5. * a21);
        aa31 = (pow(2. - 3. * a21, 2.) * (-2. + a21 * (6. + a21 * (22. + 15. * a21 * (-6. + 5. * a21))))) / (6. * pow(3. - 5. * a21, 4.) * pow(a21, 2.));
        aa32 = (pow(2. - 3. * a21, 2.) * (2. + a21 * (-6. + 5. * a21))) / (6. * pow(3. - 5. * a21, 4.) * pow(a21, 2.));
        aaa31 = (pow(2. - 3. * a21, 2.) * (-2. + 3. * a21 * (4 + a21 * (-8 + 5. * a21)))) / (6. * pow(3. - 5. * a21, 4.) * a21);

        // a31 = 0.451450450187708;
        // bb1 = (-2. + 5. * a31 * a31 * (16. + a31 * (-42. + 29. * a31))) / (60. * pow(2. - 3. * a31, 2.) * a31 * a31);
        // bb2 = pow(3. - 5. * a31, 4) / (60. * pow(2. - 3. * a31, 2) * (2 + a31 * (-6. + 5. * a31)));
        // bb3 = 1. / (60. * a31 * a31 * (2. + a31 * (-6. + 5. * a31)));
        // bbb1 = (1. + 5. * (-1 + a31) * a31) / (60. * a31 * (-2. + 3. * a31));
        // a21 = (2. - 3. * a31) / (3. - 5. * a31);
        // aa31 = (a31 * a31 * (-6. + a31 * (78. + a31 * (-248. + 25. * (12. - 5. * a31) * a31)))) / (6. * pow(2. - 3. * a31, 2));
        // aa32 = (pow(3. - 5. * a31, 2.) * a31 * a31 * (2. + a31 * (-6. + 5. * a31))) / (6. * pow(2. - 3. * a31, 2));
        // aaa31 = (a31 * a31 * (6. + a31 * (-30. + (48. - 25. * a31) * a31))) / (6. * (-2. + 3. * a31));
        // B1 = 1., B2 = 0., B3 = 0.;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            f_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A2 * dt * dt * Lt_n0[j] + A3 * dt * dt * dt * TLt_n0[j];
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            TLt_n1[j] = dfdxxx(f_n1, j);
            f_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) +
                      dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]) +
                      dt * dt * dt * (aaa31 * TLt_n0[j] + aaa32 * TLt_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            TLt_n2[j] = dfdxxx(f_n2, j);
            u_nn[j] = u_n0[j] + dt * (B1 * L_n0[j] + B2 * L_n1[j] + B3 * L_n2[j]) +
                      dt * dt * (bb1 * Lt_n0[j] + bb2 * Lt_n1[j] + bb3 * Lt_n2[j]) +
                      dt * dt * dt * (bbb1 * TLt_n0[j] + bbb2 * TLt_n1[j] + bbb3 * TLt_n2[j]);
        }
    }

    // //-----------------------------------------
    void compt_2_stage_int()
    {
        int i, j;
        a21 = 0.765; // 0.532  0.468; //0.765;   // 1. / 210. * (93. + sqrt(2139.));
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);

            TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt23_LLttnn_line_SSP()
    {

        int i, j;

        a21 = 0.594223212099088, v2 = 0.306027487008159;
        ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.0;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            // Lt_n0[j] = -1. * dfdxcent(L_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt24_LLttnn_line_SSP()
    {
        int i, j;

        A1 = 0.5, B0 = 1., B1 = 0., C0 = 1. / 3., C1 = 2. / 3.;
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A1 * A1 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_nn[j] = u_n0[j] + dt * (B0 * L_n0[j] + B1 * L_n1[j]) + dt * dt / 2.0 * (C0 * Lt_n0[j] + C1 * Lt_n1[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt34_LLttnn_line_SSP()
    {
        int i, j;
        // A1 = 1. / 4., A2 = 1. / 2., B0 = 1. / 3.; //稳定性足够  B0 = 1. / 5.;  B0 = 1. / 3.;
        // A1 = 1. / 4., A2 = 2. / 3., B0 = 1. / 3.;
        A1 = 0.5, A2 = 0.8, B0 = 0.6;
        B1 = -((pow(A2, 3) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3)));
        B2 = -((pow(A1, 3) * (-1. + B0)) / (pow(A1, 3) - pow(A2, 3)));
        C0 = -((-1. + 2. * A1 + 2. * A2 - 6. * A1 * A2) / (6. * A1 * A2)) + (A2 * (-pow(A1, 3) + A1 * A2 * A2) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));
        C1 = -((-1. + 2. * A2) / (6. * A1 * (A1 - A2))) + (A1 * pow(A2, 3) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));
        C2 = -((1. - 2. * A1) / (6. * (A1 - A2) * A2)) - (pow(A1, 3) * A2 * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A1 * A1 / 2.0 * dt * dt * Lt_n0[j];
            u_n2[j] = u_n0[j] + A2 * dt * L_n0[j] + A2 * A2 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]), f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            L_n2[j] = -dfdx(f_n2, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            Lt_n2[j] = dfdxx(f_n2, j);

            u_nn[j] = u_n0[j] + dt * (B0 * L_n0[j] + B1 * L_n1[j] + B2 * L_n2[j]) + dt * dt / 2.0 * (C0 * Lt_n0[j] + C1 * Lt_n1[j] + C2 * Lt_n2[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
            // Tu2[j] = Tu1[j];
            // Tu1[j] = u_nn[j];
        }
    }

    void compt35_LLttnn_line_SSP()
    {
        int i, j;

        a21 = 0.7504;
        aa12 = -1. + 2. * a21;
        vv1 = (1. + 2. * a21 * (-4. + 5. * a21)) / (12. * a21 * (-3. + 5. * a21));
        vv2 = 1. / (12. * a21 * (3. + 10. * (-1. + a21) * a21));
        vv3 = (25. * pow(aa12, 3)) / (12. * (-3. + 5. * a21) * (3. + 10. * (-1 + a21) * a21));
        v1 = 1., v2 = 0., v3 = 0.;
        a31 = (3. - 5. * a21) / (5. - 10. * a21), a32 = 0.;
        aa31 = ((-3. + 5. * a21) * (-3. + 10. * a21) * (1. + 5. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));
        aa32 = ((-3. + 5. * a21) * (3. + 10. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + v3 * L_n2[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j]);

            // TL_n1a[j] = L_n1[j], TL_n1[j] = L_n0[j];
            // TLt_n1a[j] = Lt_n1[j], TLt_n1[j] = Lt_n0[j];
        }
    }

    void compt46_LLttnn_line_SSP()
    {
        int i, j;

        a21 = 0.188882174845445, aa41 = 0.106470264959241;
        a31 = 0.478952613047175, aa42 = 0.069604114078399;
        a41 = 0.793700870016536, aa43 = 0.138906156494863;
        v1 = 1., vv1 = 0.059536043418564;
        aa21 = 0.017838237987173, vv2 = 0.220386636642689;
        aa31 = 0.001782796208054, vv3 = 0.157700637816082;
        aa32 = 0.112915006564305, vv4 = 0.062376682122665;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + aa21 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            u_n3[j] = u_n0[j] + dt * (a41 * L_n0[j] + a42 * L_n1[j] + a43 * L_n2[j]) +
                      dt * dt * (aa41 * Lt_n0[j] + aa42 * Lt_n1[j] + aa43 * Lt_n2[j]);
            f_n3[j] = f_u(u_n3[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n3[j] = dfdx(f_n3, j);
            Lt_n3[j] = dfdxx(f_n3, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j] + vv4 * Lt_n3[j]);

            // TL_n1a[j] = L_n1[j], TL_n1[j] = L_n0[j];
            // TLt_n1a[j] = Lt_n1[j], TLt_n1[j] = Lt_n0[j];
        }
    }


    //------------------------------------------------
    void f_eq_u(double *ff, double *uu)
    {
        for (int i = 0; i < n1; i++)
        {

            *ff++ = f_u(*uu++);
        }
    }
    double f_u(double xx)
    {
        double y;
        // y = xx * xx / 2.0;//bergure
        y = xx; // line
        return y;
    }
    void Store_obj(double *ss, double *ff)
    {
        for (int ik = 0; ik < n1; ik++) //
        {
            *ss++ = *ff++;
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
            ofs << " VARIABLES=x,u " << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (int j = 0; j < kt; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_obj[j][i];
                }

                ofs << endl;
            }
            ofs.close();
        }

        if (j < 0.5 && j > -0.8)
        {
            // string title = "result_t_burgers_excat"; //+ to_string(t);
            string title = to_string(n);
            string Title = "excat.txt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << "Title=" << title << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                ofs << "  ";
                ofs.precision(18);
                ofs << *(u_excat + i);
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
                    ofs << put_one_n[j][i];
                }

                ofs << endl;
            }
            ofs.close();
        }
    }
    ////-----------------------------------------
    void RK_compt_t1()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j = 0, i;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n0, nb + j); //-uu_t0[nb + j]; -dfdx
            df[nb + j] = dfdx(f_n0, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_nn[nb + j] = f_u(u_nn[nb + j]);
            j++;
        }
    }

    void RK_compt_Second()
    {
        int i, j;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n0[j] = dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            TLt_n0[j] = dfdxxx(f_n0, j);
            // u_nn[j] = u_n0[j] + dt * L_n0[j] + dt * dt * Lt_n0[j] / 2.0 + dt * dt * dt * TLt_n0[j] / 6.0;
            u_nn[j] = u_n0[j] + dt * dt * Lt_n0[j];
        }
    }

    void RK3_compt_t1()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j = 0, i;
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n0, nb + j); //-uu_t0[nb + j]; -dfdx
            df[nb + j] = dfdx(f_n0, nb + j);
            u_n1[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        j = 0, i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n1, nb + j); //-uu_t1[nb + j];
            df[nb + j] = dfdx(f_n1, nb + j);
            u_n2[nb + j] = u_n0[nb + j] * 0.75 + u_n1[nb + j] * 0.25 + 0.25 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        j = 0, i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            // df[nb + j] = -df_dx(f_n2, nb + j);
            df[nb + j] = dfdx(f_n2, nb + j);
            u_nn[nb + j] = u_n0[nb + j] * 1.0 / 3.0 + u_n2[nb + j] * 2.0 / 3.0 + 2.0 / 3.0 * dt * df[nb + j];
            j++;
        }
    }

    void SSPRK54_compt()
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 0.39175222700392 * dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = 0.44437049406734 * u_n0[nb + j] + 0.55562950593266 * u_n1[nb + j] + 0.36841059262959 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = 0.62010185138540 * u_n0[nb + j] + 0.37989814861460 * u_n2[nb + j] + 0.25189177424738 * dt * df[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            Ltt_n0[nb + j] = df[nb + j];
            u_n4[nb + j] = 0.17807995410773 * u_n0[nb + j] + 0.82192004589227 * u_n3[nb + j] + 0.54497475021237 * dt * df[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }

        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n4, nb + j);
            u_nn[nb + j] = 0.00683325884039 * u_n0[nb + j] + 0.51723167208978 * u_n2[nb + j] + 0.12759831133288 * u_n3[nb + j] +
                           0.34833675773694 * u_n4[nb + j] + 0.08460416338212 * dt * Ltt_n0[nb + j] + 0.22600748319395 * dt * df[nb + j];
            j++;
        }
    }

 
    void RK65LS_compt2() // 正确 Numerical Methods for Ordinary Differential Equations
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 1. / 4. * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] + 1. / 8. * dt * L_n0[nb + j] + 1. / 8. * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = u_n0[nb + j] + 0.0 * dt * L_n0[nb + j] + 0.0 * dt * L_n1[nb + j] +
                           1. / 2. * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = u_n0[nb + j] + 3. / 16. * dt * L_n0[nb + j] - 3. / 8. * dt * L_n1[nb + j] +
                           3. / 8. * dt * L_n2[nb + j] + 9. / 16. * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n4[nb + j] = -dfdx(f_n4, nb + j); //-uu_t1[nb + j];
            u_n5[nb + j] = u_n0[nb + j] - 3. / 7. * dt * L_n0[nb + j] + 8. / 7. * dt * L_n1[nb + j] + 6. / 7. * dt * L_n2[nb + j] -
                           12. / 7. * dt * L_n3[nb + j] + 8. / 7. * dt * L_n4[nb + j];
            f_n5[nb + j] = f_u(u_n5[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n5, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + 7. / 90. * dt * L_n0[nb + j] + 0. * dt * L_n1[nb + j] +
                           16. / 45. * dt * L_n2[nb + j] + 2. / 15. * dt * L_n3[nb + j] +
                           16. / 45. * dt * L_n4[nb + j] + 7. / 90. * dt * df[nb + j];
            j++;
        }
    }

    void RK76LS_compt() // 208 Numerical Methods for Ordinary Differential Equations
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 1. / 3. * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] + 0. * dt * L_n0[nb + j] + 2. / 3. * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = u_n0[nb + j] + 1. / 12. * dt * L_n0[nb + j] + 1. / 3. * dt * L_n1[nb + j] -
                           1. / 12. * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = u_n0[nb + j] + 25. / 48. * dt * L_n0[nb + j] - 55. / 24. * dt * L_n1[nb + j] +
                           35. / 48. * dt * L_n2[nb + j] + 15. / 8. * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n4[nb + j] = -dfdx(f_n4, nb + j); //-uu_t1[nb + j];
            u_n5[nb + j] = u_n0[nb + j] + 3. / 20. * dt * L_n0[nb + j] - 11. / 24. * dt * L_n1[nb + j] - 1. / 8. * dt * L_n2[nb + j] +
                           1. / 2. * dt * L_n3[nb + j] + 1. / 10. * dt * L_n4[nb + j];
            f_n5[nb + j] = f_u(u_n5[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            f_n7[nb + j] = -dfdx(f_n5, nb + j);
            u_n6[nb + j] = u_n0[nb + j] - 261. / 260. * dt * L_n0[nb + j] + 33. / 13. * dt * L_n1[nb + j] +
                           43. / 156. * dt * L_n2[nb + j] - 118. / 39. * dt * L_n3[nb + j] +
                           32. / 195. * dt * L_n4[nb + j] + 80. / 39. * dt * f_n7[nb + j];
            f_n6[nb + j] = f_u(u_n6[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n6, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + 13. / 200. * dt * L_n0[nb + j] + 0. * dt * L_n1[nb + j] +
                           11. / 40. * dt * L_n2[nb + j] + 11. / 40. * dt * L_n3[nb + j] +
                           4. / 25. * dt * L_n4[nb + j] + 4. / 25. * dt * f_n7[nb + j] + 13. / 200. * dt * df[nb + j];
            j++;
        }
    }

    //---------------------------------------------------------

    int Perimeter()
    {
        cout << " yes ";
        return 2 * t;
    }

    double dfdx(double *g, int i)
    {
        double y;
        // 1up
        // y = (*(g + i) - *(g + i - 1)) / dx; //u_t+u_x=0
        y = (*(g + i + 1) - *(g + i)) / dx; // u_t-u_x=0
        return y;
    }

    double dfdxx(double *g, int i)
    {
        double y;
        // 2ord
        y = (*(g + i + 1) - *(g + i) * 2.0 + *(g + i - 1)) / dx / dx;
        // y = (*(g + i ) - *(g + i - 1) * 2.0 + *(g + i-2)) / dx / dx;
        return y;
    }

    double dfdxxx(double *g, int i)
    {
        double y;
        // 2ord
        y = (*(g + i + 2) - *(g + i + 1) * 3.0 + *(g + i) * 3.0 - *(g + i - 1)) / dx / dx / dx;
        // y = (*(g + i ) - *(g + i - 1) * 2.0 + *(g + i-2)) / dx / dx;
        return y;
    }

    void compt_CFL()
    {
        double *umax;
        // umax = max_element(u_n0 + 1, u_n0 + n1 - 1);
        // cout << *umax;
        // dt = CFL * dx / (*umax);
        // dt=0.02/1.0;
        dt = CFL * dx;
    }

    void maxTV()
    {
        double TV = 0.0;
        for (int i = 0; i < n1 - 1; i++) // 间断
        {
            TV = abs(u_nn[i + 1] - u_nn[i]) + TV;
        }
        TV = abs(TV - 2.0);
        // cout << " \n TV: " << TV;
        if (tt > 0.0)
        {
            if (TV > TVmax)
            {
                TVmax = TV;
            }
        }
    }

    void TVabs(int m)
    {
        double TV = 0.0;
        TV = TVmax;

        if (TV < 2.0E-16)
        {
            TV = 2.0E-16;
        }
        cout << "\n  CFL: " << CFL << ",  TV: " << TV << "   \n";
        string Title = "TV_CFL.plt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=CFL,TV " << endl;
            ofs.precision(4), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
            ofs.precision(4), ofs << CFL, ofs << "  ";
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
    // int n, oder_n = 2, begin_cell = 1600;
    int n, oder_n = 200, begin_cell = 4800;
    double t = 100.0, CFL = 0.701; // 计算时间总长，CFL数  CFL = 0.66 1.45
    n = begin_cell + 1;
    for (int m = 1; m < oder_n; m++)
    {
        comput_main(n, t, CFL, m);
        CFL = CFL + 0.01;
    }
    // system("pause");
    return 2;
}
int comput_main(int n, double t, double CFL, int m)
{
    int nb = 6; // 边界网格结点数（有限差分）
    int n1 = n + 6 * 2;
    double main_df[n + 6 * 2], u_excat[n + 6 * 2];
    double Ltt_n0[n + 6 * 2], x[n + 6 * 2];                                                          // 声明变长数组
    double L_n0[n + 6 * 2], L_n1[n + 6 * 2], L_n2[n + 6 * 2], L_n3[n + 6 * 2], L_n4[n + 6 * 2];      // f=L
    double Lt_n0[n + 6 * 2], Lt_n1[n + 6 * 2], Lt_n2[n + 6 * 2], Lt_n3[n + 6 * 2], Lt_n4[n + 6 * 2]; // L_t
    double u_n0[n + 6 * 2], u_n1[n + 6 * 2], u_n2[n + 6 * 2], u_n3[n + 6 * 2], u_n4[n + 6 * 2];      // 推进步中间时刻
    double u_n5[n + 6 * 2], u_n6[n + 6 * 2], u_n7[n + 6 * 2], u_n8[n + 6 * 2], u_n9[n + 6 * 2];      // 推进步中间时刻
    double f_n0[n + 6 * 2], f_n1[n + 6 * 2], f_n2[n + 6 * 2], f_n3[n + 6 * 2], f_n4[n + 6 * 2];      // 推进步中间时刻
    double f_n5[n + 6 * 2], f_n6[n + 6 * 2], f_n7[n + 6 * 2], f_n8[n + 6 * 2], f_n9[n + 6 * 2];      // 推进步中间时刻
    double u_nn[n + 6 * 2], f_nn[n + 6 * 2], Lttx_n0[n + 6 * 2];                                     // 下一时刻（n+1)

    double TL_n0[n + 6 * 2], TL_n1[n + 6 * 2], TL_n2[n + 6 * 2], TL_n3[n + 6 * 2], TL_n4[n + 6 * 2];
    double TLt_n0[n + 6 * 2], TLt_n1[n + 6 * 2], TLt_n2[n + 6 * 2], TLt_n3[n + 6 * 2], TLt_n4[n + 6 * 2];
    double Tu1[n + 6 * 2], Tu2[n + 6 * 2];
    class cfda cff;
    cff.kt = 1;
    cff.A1 = 1 / 3.0;
    cff.A2 = 2 / 3.0;
    cff.t = t;           // 1.0 / cff.PI; //计算总时间
    cff.nx_begin = -2.0; // 网格右端点
    cff.nx_end = 1.0;    // 网格左端点
    cff.CFL = CFL, cff.TVmax = 0.0;

    cff.n = n, cff.nb = nb, cff.n1 = n1;
    cff.Ltt_n0 = Ltt_n0, cff.x = x;
    cff.u_n0 = u_n0, cff.u_n1 = u_n1, cff.u_n2 = u_n2, cff.u_n3 = u_n3, cff.u_n4 = u_n4;                     // 推进步中间时刻
    cff.u_n5 = u_n5, cff.u_n6 = u_n6, cff.u_n7 = u_n7, cff.u_n8 = u_n8, cff.u_n9 = u_n9;                     // 推进步中间时刻
    cff.f_n0 = f_n0, cff.f_n1 = f_n1, cff.f_n2 = f_n2, cff.f_n3 = f_n3, cff.f_n4 = f_n4;                     // 推进步中间时刻
    cff.f_n5 = f_n5, cff.f_n6 = f_n6, cff.f_n7 = f_n7, cff.f_n8 = f_n8, cff.f_n9 = f_n9;                     // 推进步中间时刻
                                                                                                             // 声明变长数组
    cff.u_nn = u_nn, cff.f_nn = f_nn, cff.Lttx_n0 = Lttx_n0;                                                 // 下一时刻（n+1)
    cff.L_n0 = L_n0, cff.L_n1 = L_n1, cff.L_n2 = L_n2, cff.L_n3 = L_n3, cff.L_n4 = L_n4;                     // f=L
    cff.Lt_n0 = Lt_n0, cff.Lt_n1 = Lt_n1, cff.Lt_n2 = Lt_n2, cff.Lt_n3 = Lt_n3, cff.Lt_n4 = Lt_n4;           // L_t
                                                                                                             // 存储一时刻（n-1)
    cff.TL_n0 = TL_n0, cff.TL_n1 = TL_n1, cff.TL_n2 = TL_n2, cff.TL_n3 = TL_n3, cff.TL_n4 = TL_n4;           // f=TL
    cff.TLt_n0 = TLt_n0, cff.TLt_n1 = TLt_n1, cff.TLt_n2 = TLt_n2, cff.TLt_n3 = TLt_n3, cff.TLt_n4 = TLt_n4; // TL_t

    cff.Tu1 = Tu1, cff.Tu2 = Tu2;
    cff.df = main_df, cff.u_excat = u_excat, cff.intc();
    for (cff.ij = 0; cff.ij < cff.kt; cff.ij++)
    {
        // cff.nt = 1000 * pow(2, cff.ij); //计算时间步数

        comput_tt(cff);

        cff.Store_obj(cff.put_obj[cff.ij], cff.u_nn); // 存储数据
        // cff.Write_obj(cff.ij + 1);                    //存储数据
    }
    cff.Write_obj(-1);
    cff.TVabs(m);
    return 2;
}
////-----------------------------------------
void comput_tt(cfda &cff1)
{
    int n_all, i3 = 0, tt_flag = 0;
    cff1.tt = 0, n_all = cff1.n1, cff1.border();
    cff1.f_eq_u(cff1.f_n0, cff1.u_n0), cff1.TVmax = 0.0;
    for (int i1 = 0; i1 < 500; i1++) // 500 200 50
    {
        cff1.compt_CFL();
        cff1.tt = cff1.tt + cff1.dt;
        if (i1 < 20000000)
        {
 
            // cff1.compt_ThD_13_line_SSP();
            cff1.compt_ThD_23_line_SSP();
            // cff1.compt_ThD_24_line_SSP();
            // cff1.compt_ThD_25_line_SSP();
            // cff1.compt_ThD_36_line_SSP();

            // cff1.compt23_LLttnn_line_SSP();
            // cff1.compt24_LLttnn_line_SSP(); // cff1.compt34_LLttnn_line_SSP();
            // cff1.compt35_LLttnn_line_SSP();
            // cff1.compt46_LLttnn_line_SSP();

            // cff1.RK_compt_t1(); // 3-RK
            // cff1.RK_compt_Second(); // 3-RK
            // cff1.RK3_compt_t1(); // 3-RK

            // cff1.SSPRK54_compt(); //
            // cff1.RK65LS_compt2();// cff1.RK65LS_compt();
            // cff1.RK76LS_compt(); //
        }
 

        //___________________________________________________________________
        if (cff1.tt > cff1.t)
        {
            cff1.f_eq_u(cff1.f_nn, cff1.u_nn);
            cout.precision(18); // 精度为18，正常为6
            cout << " \n comput time: " << cff1.tt - cff1.dt << " comput time n:   " << i1;
            break;
        }
        cff1.carr(cff1.u_n0, cff1.u_nn, n_all);
        cff1.f_eq_u(cff1.f_n0, cff1.u_nn);

        // cout << " \n comput time: " << cff1.tt << " comput time n:   " << i1;

        cff1.border();
        cff1.maxTV();

        if (i3 * cff1.t / cff1.onet_out < cff1.tt + 0.00000001)
        {
            // cff1.Store_obj(cff1.put_one_n[i3], cff1.u_nn); //存储数据
            i3++;
        }
    }
}
