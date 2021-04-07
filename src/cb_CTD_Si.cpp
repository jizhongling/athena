#include <XML/Helper.h>
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
//////////////////////////////////
// Central Barrel Tracker Silicon
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
  double     SizeZCut   = dim.zmax();
  double     SiLayerGap = dim.gap();

  Material   Vacuum     = desc.material("Vacuum");

  // Create Global Volume 
  Tube cb_CTD_GVol_Solid(RIn, ROut, SizeZ / 2.0, 0., 360.0 * deg);
  Volume detVol("cb_CTD_GVol_Logic", cb_CTD_GVol_Solid, Vacuum);
  detVol.setVisAttributes(desc.visAttributes(x_det.visStr()));

  // Construct Silicon Layers
  xml_comp_t x_layer = x_det.child(_U(layer));
  const int repeat   = x_layer.repeat();
  
  xml_comp_t x_slice = x_layer.child(_U(slice));
  Material slice_mat = desc.material(x_slice.materialStr());
  
  double layerRIn[100];
  double layerROut[100];
  
  // Loop over layers
  for(int i = 0; i < repeat; i++) {
    layerRIn[i]  = RIn + (SiLayerGap * i);
    layerROut[i] = RIn + (0.01 + SiLayerGap * i);

    if (layerROut[i] > ROut)
      continue; 

    string logic_layer_name = detName + _toString(i, "_Logic_lay_%d");
    Volume layerVol(logic_layer_name,Tube(layerRIn[i], layerROut[i], SizeZ / 2.0, 0.0, 360.0 * deg), slice_mat);
    layerVol.setVisAttributes(desc,x_layer.visStr());
    sens.setType("tracker");
    layerVol.setSensitiveDetector(sens);

    Position     layer_pos = Position(0.0, 0.0, 0.0);
    PlacedVolume layerPV = detVol.placeVolume(layerVol, layer_pos);
    layerPV.addPhysVolID("layer", i+1);
  }

  DetElement   det(detName, detID);
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(0.0, 0.0, 0.0), Position(0.0, 0.0, 0.0));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  detPV.addPhysVolID("barrel", 1);
  det.setPlacement(detPV);
  return det;
}

DECLARE_DETELEMENT(cb_CTD_Si, createDetector)
