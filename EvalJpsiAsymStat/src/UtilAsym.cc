#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <interface_main/SQDimuon.h>
#include "UtilAsym.h"
using namespace std;

void UtilAsym::CalcAngle(const SQDimuon* dim, const int pol_sign, double& xb, double& xt, double& phi_s_tf)
{
  CalcAngle(dim->get_mom_pos(), dim->get_mom_neg(), pol_sign, xb, xt, phi_s_tf);
}

void UtilAsym::CalcAngle(const TVector3& mu0, const TVector3& mu1, const int pol_sign, double& xb, double& xt, double& phi_s_tf)
{
  const double M_MU = 0.10565836668;
  CalcAngle(TLorentzVector(mu0, M_MU), TLorentzVector(mu1, M_MU), pol_sign, xb, xt, phi_s_tf);
}

/// Convert the momenta of two muons in the lab frame to "phi_s_tf", "xb" and "xt".
/** 
 * Imported from "asymmetry.C" written by Abinash.
 */
void UtilAsym::CalcAngle(const TLorentzVector& mu0, const TLorentzVector& mu1, const int pol_sign, double& xb, double& xt, double& phi_s_tf)
{
  //Initial 4-vector definition
  const float pro_mass = 0.93827;
  const float beamE = 120.0;
  TLorentzVector hadron_b, hadron_t;
  hadron_b.SetPxPyPzE(0.0,0.0,TMath::Sqrt(beamE*beamE-pro_mass*pro_mass),beamE);
  hadron_t.SetPxPyPzE(0.0,0.0,0.0,pro_mass);
  TLorentzVector initial = hadron_b + hadron_t;

  TLorentzVector spin;
  spin.SetPxPyPzE(0.0, (pol_sign > 0 ? +1.0 : -1.0), 0.0, 0.0);

  TLorentzVector dimu = mu0 + mu1;
  float dimu_gmass = dimu.Mag();
    
  //boost to centre of mass
  TLorentzRotation cm_r;
  TVector3 cm_v = initial.BoostVector();
  cm_r.Boost(-cm_v);
          
  TLorentzVector hadron_b_cm = cm_r * hadron_b;
  TLorentzVector hadron_t_cm = cm_r * hadron_t;
  TLorentzVector mu0_cm = cm_r * mu0;
  TLorentzVector mu1_cm = cm_r * mu1;
  TLorentzVector dimu_cm = cm_r * dimu;
  TLorentzVector spin_cm = cm_r * spin;
          
  float s = 2.0*beamE*pro_mass + 2.0*pro_mass*pro_mass;
  float sqrt_s = TMath::Sqrt(s);
  float dimu_E = dimu_cm.E();
  float dimu_pl = dimu_cm.Pz();
  float dimu_pl_max = sqrt_s/2.0*(1-(dimu_gmass*dimu_gmass)/s);
  float xf = dimu_pl/dimu_pl_max;
  float y = 0.5*TMath::Log((dimu_E+dimu_pl)/(dimu_E-dimu_pl));
  float tau = (dimu_gmass*dimu_gmass)/s;
  xb = 0.5*(TMath::Sqrt(xf*xf+4*tau)+xf);
  xt = 0.5*(TMath::Sqrt(xf*xf+4*tau)-xf);
          
  //boost to dilepton restframe
  TLorentzRotation cm_dimu;
  TVector3 cm_v2 = dimu_cm.BoostVector();
  cm_dimu.Boost(-cm_v2);
          
  TLorentzVector hadron_b_cm_dimu = cm_dimu * hadron_b_cm;
  TLorentzVector hadron_t_cm_dimu = cm_dimu * hadron_t_cm;
  TLorentzVector mu0_cm_dimu = cm_dimu * mu0_cm;
  TLorentzVector mu1_cm_dimu = cm_dimu * mu1_cm;
  TLorentzVector dimu_cm_dimu = cm_dimu * dimu_cm;
  TLorentzVector spin_cm_dimu = cm_dimu * spin_cm;
          
  TVector3 zaxis = (hadron_b_cm_dimu.Vect()).Unit()-(hadron_t_cm_dimu.Vect()).Unit();
  TVector3 zunit = zaxis.Unit();
  TVector3 hadron_plane_norm = zunit.Cross((hadron_b_cm_dimu.Vect()).Unit());
  TVector3 lepton_plane_norm = zunit.Cross((mu0_cm_dimu.Vect()).Unit());
  TVector3 yunit = lepton_plane_norm.Unit();
  TVector3 xunit = (yunit.Cross(zunit)).Unit();

  //below are phi_s_tf from target rest frame. The frame for the Asymmetry calculation
  TVector3 spin_vect_tf = spin.Vect();
  TVector3 z_prime_unit = (hadron_b.Vect()).Unit();
  TVector3 y_prime_unit = (z_prime_unit.Cross(dimu.Vect())).Unit(); 
  TVector3 x_prime_unit = (y_prime_unit.Cross(z_prime_unit)).Unit(); 
  float cosphi_s_tf = (x_prime_unit.Dot(spin_vect_tf))/((x_prime_unit.Mag())*(spin_vect_tf.Mag()));
  //float phi_s_tf;
  if ( acos( spin_vect_tf.Dot(y_prime_unit) / ( (spin_vect_tf).Mag() * y_prime_unit.Mag() ) ) > TMath::Pi()/2 )
  {
    phi_s_tf = -acos(cosphi_s_tf);
  }
  else
  {
    phi_s_tf = acos(cosphi_s_tf);
  } 
}
