#ifndef _UTIL_ASYM__H_
#define _UTIL_ASYM__H_
class SQDimuon;
class SQDimuonVector;
class TVector3;
class TLorentzVector;

namespace UtilAsym {
  void CalcAngle(const SQDimuon* dim, const int pol_sign, double& xb, double& xt, double& phi_s_tf);
  void CalcAngle(const TVector3& mu0, const TVector3& mu1, const int pol_sign, double& xb, double& xt, double& phi_s_tf);
  void CalcAngle(const TLorentzVector& mu0, const TLorentzVector& mu1, const int pol_sign, double& xb, double& xt, double& phi_s_tf);
}

#endif // _UTIL_ASYM__H_
