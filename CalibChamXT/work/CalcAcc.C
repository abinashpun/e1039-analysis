/// @file CalcAcc.C: Simple ROOT macro to calculate the tracklet acceptance.
using namespace std;

void HodoPosition(const int idx, double& x, double& y, double& z)
{
  switch (idx) {
  case 1: x = 162/2; y =  70; z =  683; break;
  case 2: x = 203/2; y = 152; z = 1420; break;
  case 3: x = 228/2; y = 168; z = 1958; break;
  case 4: x = 304/2; y = 183; z = 2240; break;
  }
}

void CalcOneChamber(const string name, const int hodo1, const int hodo2, const double x_ch, const double y_ch, const double z_ch)
{
  double x1, y1, z1;
  HodoPosition(hodo1, x1, y1, z1);
  double x2, y2, z2;
  HodoPosition(hodo2, x2, y2, z2);

  double  x_max = x1 + (z_ch - z1) * (x2 - x1) / (z2 - z1);
  double  y_max = y1 + (z_ch - z1) * (y2 - y1) / (z2 - z1);
  double tx_max = (x1 + x2) / (z2 - z1);
  double ty_max = (y1 + y2) / (z2 - z1);

  if (x_max > x_ch) x_max = x_ch;
  if (y_max > y_ch) y_max = y_ch;

  cout << setw(6) << name << " "
       << setw(8) <<  x_max << " " << setw(8) <<  y_max << " "
       << setw(8) << tx_max << " " << setw(8) << ty_max
       << endl;
}

void CalcAcc()
{
  double x_d0 = 102/2;  double y_d0 = 122/2;  double z_d0 =  620;
  double x_d2 = 233/2;  double y_d2 = 264/2;  double z_d2 = 1345;
  double x_d3 = 320/2;  double y_d3 = 166  ;  double z_d3 = 1900;

  CalcOneChamber("D0", 1, 2, x_d0, y_d0, z_d0);
  CalcOneChamber("D2", 2, 4, x_d2, y_d2, z_d2);
  CalcOneChamber("D3", 2, 4, x_d3, y_d3, z_d3);

  exit(0);
}

/*

use user_e1039_geom_plane

select det_name, x0, y0, z0, width, height 
  from param_2022022101 
 where det_name like 'H%T' or det_name like 'H%B' order by z0;
+----------+-----------+----------+---------+---------+--------+
| det_name | x0        | y0       | z0      | width   | height |
+----------+-----------+----------+---------+---------+--------+
| H1B      |  -0.76518 |  -35.062 |     683 | 162.008 |  69.85 |
| H1T      |  -0.83482 |   34.788 |     683 | 162.008 |  69.85 |
| H2B      | -0.937406 | -76.0406 | 1420.95 | 203.238 |    152 |
| H2T      |  -1.38415 |  75.9594 | 1421.27 | 203.238 |    152 |
| H3B      | 0.0165354 | -84.1908 | 1958.34 | 227.518 | 167.64 |
| H3T      |  0.105385 |  83.4492 |  1958.9 | 227.518 | 167.64 |
| H4T      | -0.294042 |  90.7328 | 2234.29 | 304.518 | 182.88 |
| H4B      | -0.274918 | -92.0383 | 2251.71 | 304.518 | 182.88 |
+----------+-----------+----------+---------+---------+--------+

select det_name, x0, y0, z0, width, height
  from param_2022022101
 where det_name like 'D%X' order by z0;
+----------+-----------+------------+---------+---------+--------+
| det_name | x0        | y0         | z0      | width   | height |
+----------+-----------+------------+---------+---------+--------+
| D0X      |    -0.552 |      2.743 | 623.432 |   101.6 | 121.92 |
| D1X      |     0.349 |     -0.173 | 689.814 |   152.4 | 137.16 |
| D2X      | -0.816646 | -0.0619792 | 1347.34 | 233.274 | 264.16 |
| D3mX     |   -2.6796 |   -79.5878 | 1894.71 |     320 |    166 |
| D3pX     |  -1.02271 |    78.6947 |  1931.3 |     320 |    166 |
+----------+-----------+------------+---------+---------+--------+

*/
