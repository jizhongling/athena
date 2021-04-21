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

  xml::DetElement detElem   = handle;
  std::string     detName   = detElem.nameStr();
  int             detID     = detElem.id();
  xml::Component  dims      = detElem.dimensions();
  xml::Component  pos       = detElem.position();
  double          SizeZ     = dims.z_length(); // Size in Z direction
  double          ROut      = dims.rmax();     // Outer radius
  double          RIn       = dims.rmin();     // Inner radius
  double          ShiftZ    = dims.z_offset();
  double          X         = dims.x();
  double          Z         = dims.z();
  int             Nlayers   = dims.number();
  double          HCAL_rmin = dims.rmax1(); // Maximum radius that the layer can be
  Material        mat       = desc.material(detElem.materialStr());
  Material        vac       = desc.material("Vacuum");

  // Outer Volume
  Tube   ci_GEM_GVol_Solid(RIn, ROut, SizeZ / 2., 0., 360 * deg);
  Volume detVol("ci_GEM_GVol_Logic", ci_GEM_GVol_Solid, vac);

  // Adding layers to placed volume
  for (xml_coll_t li(detElem, _U(layer)); li; ++li) {
    xml_comp_t  x_layer    = li;
    std::string layer_name = detName + _toString(x_layer.id(), "_layer%d");
    double      outer_r    = x_layer.outer_r();
    if (outer_r > HCAL_rmin) {
      outer_r = HCAL_rmin;
    }

    Volume layer_vol(layer_name, Tube(x_layer.inner_r(), outer_r, x_layer.dz()), mat);
    layer_vol.setVisAttributes(desc.visAttributes(detElem.visStr()));
    sens.setType("tracker");
    layer_vol.setSensitiveDetector(sens);
    Position     layer_pos(0, 0, x_layer.z());
    PlacedVolume layer_phv = detVol.placeVolume(layer_vol, layer_pos);
    layer_phv.addPhysVolID("layer", x_layer.id());
  }

  DetElement   det(detName, detID);
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(0.0, 0.0, 0.0), Position(pos.x(), pos.x(), pos.z() + SizeZ / 2.0));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  detPV.addPhysVolID("system", detID);
  det.setPlacement(detPV);

  return det;
}
// clang-format off
DECLARE_DETELEMENT(ci_GEM, createDetector)
