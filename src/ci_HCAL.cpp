#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include <XML/Helper.h>
///////////////////////////
// Central Ion GEM
///////////////////////////

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{

  xml::DetElement detElem  = handle;
  std::string     detName  = detElem.nameStr();
  int             detID    = detElem.id();
  xml::Component  dims     = detElem.dimensions();
  double          RIn      = dims.rmin();
  double          ROut     = dims.rmax();
  double          SizeZ    = dims.z_length();
  double          ShiftZ   = dims.z_offset();
  double          PosZ     = dims.z();
  double          lay_RIn  = dims.rmin1();
  double          lay_ROut = dims.rmax1();
  double          lay_dz   = dims.dz();
  Material        mat_iron = desc.material("Iron");
  Material        mat_vac  = desc.material("Vacuum");

  // Outer volume
  Tube   ci_Hcal_Solid(RIn, ROut, SizeZ / 2., 0., 360 * dd4hep::deg);
  Volume envelopeVol("ci_Hcal_Logic", ci_Hcal_Solid, mat_vac);

  // Iron tube for the layers
  Tube ci_Hcal_detSolid(lay_RIn, lay_ROut, lay_dz / 2., 0., 360 * dd4hep::deg);

  // Adding layers to placed detector volume
  for (xml_coll_t li(detElem, _U(layer)); li; ++li) {
    xml_comp_t  x_layer    = li;
    std::string layer_name = detName + _toString(x_layer.id(), "_layer%d");
    Volume      layer_vol(layer_name, ci_Hcal_detSolid, mat_iron);
    layer_vol.setVisAttributes(detElem.visStr());
    sens.setType("calorimeter");
    layer_vol.setSensitiveDetector(sens);
    Position     layer_pos(0, 0, x_layer.z());
    PlacedVolume layer_phv = envelopeVol.placeVolume(layer_vol, layer_pos);
    layer_phv.addPhysVolID("layer", x_layer.id());
  }

  DetElement   det(detName, detID);
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(0, 0, 0), Position(0, 0, ShiftZ));
  PlacedVolume detPV = motherVol.placeVolume(envelopeVol, tr);
  detPV.addPhysVolID("system", detID);
  det.setPlacement(detPV);

  return det;
}
// clang-format off
DECLARE_DETELEMENT(ci_HCAL, createDetector)
