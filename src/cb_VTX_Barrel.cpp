#include <XML/Helper.h>
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
//////////////////////////////////
// Central Barrel Vertex Detector
//////////////////////////////////

using namespace std;
using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h e, SensitiveDetector sens)
{
  xml_det_t  x_det      = e;
  string     detName    = x_det.nameStr();
  int        detID      = x_det.id();

  xml_dim_t  dim        = x_det.dimensions();
  double     RIn        = dim.rmin();
  double     ROut       = dim.rmax();
  double     SizeZ      = dim.length();

  xml_dim_t  pos        = x_det.position();

  Material   Vacuum     = desc.material("Vacuum");

  // Create Global Volume 
  Tube cb_VTX_Barrel_GVol_Solid(RIn, ROut, SizeZ / 2.0, 0., 360.0 * deg);
  Volume detVol("cb_VTX_Barrel_GVol_Logic", cb_VTX_Barrel_GVol_Solid, Vacuum);
  detVol.setVisAttributes(desc.visAttributes(x_det.visStr()));

  //////////////////
  // Barrel Ladder
  //////////////////
  xml_comp_t x_layer = x_det.child(_U(layer));
  const int repeat   = x_layer.repeat();
  
  xml_comp_t x_slice = x_layer.child(_U(slice));
  Material slice_mat = desc.material(x_slice.materialStr());

  double x = 0.0 * cm;
  double y = 0.0 * cm;
  double z = 0.0 * cm;
  int FDIV = 0;
  double dR = 0.0;
  double length = 0.0;
  double phi = 0.0;
  // Ladder Layer Parameters
  double lay_Dx[6];
  double lay_Dy[6];
  double lay_Dz[6];
  double lay_Rin[6];
  lay_Dx[0] = 0.050 * mm; lay_Dy[0] = 2.0 * cm; lay_Dz[0] = 10.0 * cm; lay_Rin[0] = 3.5 * cm;
  lay_Dx[1] = 0.050 * mm; lay_Dy[1] = 2.0 * cm; lay_Dz[1] = 11.0 * cm; lay_Rin[1] = 4.5 * cm;
  lay_Dx[2] = 0.150 * mm; lay_Dy[2] = 4.0 * cm; lay_Dz[2] = 18.0 * cm; lay_Rin[2] = 6.5 * cm;
  lay_Dx[3] = 0.150 * mm; lay_Dy[3] = 4.0 * cm; lay_Dz[3] = 24.0 * cm; lay_Rin[3] = 10.5 * cm;
  lay_Dx[4] = 0.150 * mm; lay_Dy[4] = 4.0 * cm; lay_Dz[4] = 36.0 * cm; lay_Rin[4] = 13.5 * cm;
  lay_Dx[5] = 0.150 * mm; lay_Dy[5] = 4.0 * cm; lay_Dz[5] = 48.0 * cm; lay_Rin[5] = 15.5 * cm;

  int i_layer = 0;
  int i_module = 0;
  // Loop over layers
  for(int i = 0; i < repeat; i++) {
    double cb_VTX_Barrel_ladder_DZ = lay_Dz[i];
    double cb_VTX_Barrel_ladder_DY = lay_Dy[i];
    double cb_VTX_Barrel_ladder_Thickness = lay_Dx[i];
    dR = lay_Rin[i]; 
    length = 2.0 * 3.1415 * dR;
    int laddersCount = length / cb_VTX_Barrel_ladder_DY;

    for (int i = 0; i < 2; i++) {
      double LN = cb_VTX_Barrel_ladder_DY * laddersCount;
      double LN1 = cb_VTX_Barrel_ladder_DY * (laddersCount + 1.0 + i);

      if (LN/LN1 > 0.8)
        laddersCount = laddersCount + 1;
    }

    double cb_VTX_Barrel_ladder_deltaphi = 2.0 * 3.1415926 / laddersCount;

    string ladderBoxName = detName + _toString(i, "_ladder_Solid_%d");
    string ladderName = detName + _toString(i, "_ladder_Logic_%d");
    Volume ladderVol(ladderName, Box(cb_VTX_Barrel_ladder_Thickness * 0.5, cb_VTX_Barrel_ladder_DY * 0.5, cb_VTX_Barrel_ladder_DZ * 0.5), slice_mat);
    ladderVol.setVisAttributes(desc,x_layer.visStr());
    sens.setType("tracker");
    ladderVol.setSensitiveDetector(sens);
    i_layer++;

    for (int ia = 0; ia < laddersCount; ia++) {
      phi = (ia * (cb_VTX_Barrel_ladder_deltaphi));
      x = - dR * cos(phi);
      y = - dR * sin(phi);

      RotationZYX ladder_rot = RotationZYX(cb_VTX_Barrel_ladder_deltaphi * ia, 0.0, 0.0);
      Position ladder_pos = Position(x, y, z);
      string ladderName = detName + _toString(i, "_ladder_Phys_%d") + _toString(ia, "_%d"); 
      PlacedVolume ladderPV = detVol.placeVolume(ladderVol, Transform3D(ladder_rot, ladder_pos));
      i_module++;
      ladderPV.addPhysVolID("layer", i_layer).addPhysVolID("module", i_module);
    }
  }

  // TODO: Pixels

  DetElement   det(detName, detID);
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, pos.z()));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  detPV.addPhysVolID("barrel", 1);
  det.setPlacement(detPV);
  return det;
}

DECLARE_DETELEMENT(cb_VTX_Barrel, createDetector)
