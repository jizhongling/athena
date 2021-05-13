#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
//////////////////////////////////
// Support structure for ALl-silicon
//////////////////////////////////

using namespace std;
using namespace dd4hep;

// Info from
// https://github.com/reynier0611/g4lblvtx/blob/master/source/AllSi_vtx_serv_2lyr_Detector.cc
// TODO: this is quite incomplete, should probably wait for official word
//       from he tracking WG

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();
  bool      reflect = x_det.reflect();
  const int sign    = reflect ? -1 : 1;

  // second vertexing layer
  std::vector<double> z_det   = {15 * cm, 20 * cm};
  std::vector<double> rin_l2  = {5.48 * cm, 14.8 * cm};
  std::vector<double> rout_l2 = {0, 0};

  // first vertexing layer
  std::vector<double> rin_l1  = {3.30 * cm, 14.36 * cm};
  std::vector<double> rout_l1 = {0, 0};

  const int nzplanes_l2 = z_det.size();
  const int nzplanes_l1 = z_det.size();

  for (int i = 0; i < nzplanes_l2; i++) {
    rout_l2[i] = rin_l2[i] + 0.44;
    z_det[i] *= sign / abs(sign);
  }
  for (int i = 0; i < nzplanes_l1; i++) {
    rout_l1[i] = rin_l1[i] + 0.44;
    z_det[i] *= sign / abs(sign);
  }
  // mother volume
  std::vector<double> rin_mo  = rin_l1;
  std::vector<double> rout_mo = rout_l2;

  DetElement det(detName, detID);
  Material   Vacuum = desc.material("Vacuum");
  Polycone   empty_cone("empty_cone", 0.0, 360 * degree, z_det, rin_mo, rout_mo);
  Volume     detVol("empty_cone", empty_cone, Vacuum);
  detVol.setVisAttributes(desc.invisible());

  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, 0.));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  detPV.addPhysVolID("barrel", 1);
  det.setPlacement(detPV);

  Material Al       = desc.material("Al");
  Material Graphite = desc.material("Graphite");

  // cb_DIRC_bars_Logic.setVisAttributes(desc.visAttributes(x_det.visStr()));

  Polycone polycone_l2("polycone_l2", 0, 360 * degree, z_det, rin_l2, rout_l2);
  Volume   logical_l2("polycone_l2_logic", polycone_l2, Al);
  logical_l2.setVisAttributes(desc.visAttributes(x_det.visStr()));
  detVol.placeVolume(logical_l2, tr);
  Polycone polycone_l1("polycone_l1", 0, 360 * degree, z_det, rin_l1, rout_l1);
  Volume   logical_l1("polycone_l1_logic", polycone_l1, Al);
  logical_l1.setVisAttributes(desc.visAttributes(x_det.visStr()));
  detVol.placeVolume(logical_l1, tr);

  return det;
}

DECLARE_DETELEMENT(allsilicon_support, createDetector)
