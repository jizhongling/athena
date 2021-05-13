#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
//////////////////////////////////
// Central Barrel DIRC
//////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det   = e;
  string    detName = x_det.nameStr();
  int       detID   = x_det.id();

  xml_dim_t dim   = x_det.dimensions();
  xml_dim_t pos   = x_det.position();
  double    RIn   = dim.rmin();
  double    ROut  = dim.rmax();
  double    SizeZ = dim.length();

  Material Vacuum = desc.material("Vacuum");
  Tube     cb_DIRC_Barrel_GVol_Solid(RIn, ROut, SizeZ / 2.0, 0., 360.0 * deg);
  Volume   detVol("cb_DIRC_GVol_Solid_Logic", cb_DIRC_Barrel_GVol_Solid, Vacuum);
  detVol.setVisAttributes(desc.invisible());

  DetElement   det(detName, detID);
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, pos.z()));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  detPV.addPhysVolID("barrel", 1);
  det.setPlacement(detPV);

  //////////////////
  // DIRC Bars
  //////////////////

  double dR = 83.65 * cm;

  double cb_DIRC_bars_DZ = SizeZ;
  double cb_DIRC_bars_DY = 42. * cm;
  double cb_DIRC_bars_DX = 1.7 * cm;
  double myL             = 2 * M_PI * dR;
  int    NUM             = myL / cb_DIRC_bars_DY;

  double cb_DIRC_bars_deltaphi = 2 * 3.1415926 / NUM;

  Material cb_DIRC_bars_Material = desc.material("Quartz");

  Box    cb_DIRC_bars_Solid(cb_DIRC_bars_DX / 2., cb_DIRC_bars_DY / 2., cb_DIRC_bars_DZ / 2.);
  Volume cb_DIRC_bars_Logic("cb_DIRC_bars_Logix", cb_DIRC_bars_Solid, cb_DIRC_bars_Material);
  cb_DIRC_bars_Logic.setVisAttributes(desc.visAttributes(x_det.visStr()));
  sens.setType("photoncounter");
  cb_DIRC_bars_Logic.setSensitiveDetector(sens);

  for (int ia = 0; ia < NUM; ia++) {
    double phi = (ia * (cb_DIRC_bars_deltaphi));
    double x   = -dR * cos(phi);
    double y   = -dR * sin(phi);

    Transform3D  tr(RotationZ(cb_DIRC_bars_deltaphi * ia), Position(x, y, 0));
    PlacedVolume barPV = detVol.placeVolume(cb_DIRC_bars_Logic, tr);
    barPV.addPhysVolID("module", ia);
  }
  return det;
}

DECLARE_DETELEMENT(cb_DIRC, createDetector)
