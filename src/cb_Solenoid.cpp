#include <XML/Helper.h>
///////////////////////////
// Central Barrel Solenoid
///////////////////////////

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string     detName = detElem.nameStr();
  int             detID   = detElem.id();

  xml::Component dims = detElem.dimensions();
  double SizeZ  = dims.z();    // Size in Z direction
  double ROut   = dims.rmax(); // Outer diameter
  double RIn    = dims.rmin(); // Inner diameter
  double ShiftZ = dims.delta();
  // TODO: Material
  //Material mat  = desc.material("Steel235");
  Material mat = desc.material(detElem.materialStr());

  Tube   cb_Solenoid_GVol_Solid(RIn, ROut, SizeZ / 2., 0., 360 * deg);
  Volume detVol("cb_Solenoid_GVol_Logic", cb_Solenoid_GVol_Solid, mat);

  detVol.setVisAttributes(desc.visAttributes(detElem.visStr()));

  DetElement   det(detName, detID);
  Volume       motherVol = desc.pickMotherVolume(det);
  Transform3D  tr(RotationZYX(0., 0., 0.), Position(0., 0., ShiftZ));
  PlacedVolume detPV = motherVol.placeVolume(detVol, tr);
  det.setPlacement(detPV);
  return det;
}

// clang-format off
DECLARE_DETELEMENT(cb_Solenoid, createDetector)
