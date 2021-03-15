#include <XML/Helper.h>

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens) {
  xml::DetElement detElem = handle;
  std::string     detName = detElem.nameStr();
  int             detID   = detElem.id();

  xml::Component dims      = detElem.dimensions();
  double         rInner    = dims.inner_radius();
  double         rMin      = dims.rmin();
  double         thickness = dims.thickness();
  double         innerZ    = dims.inner_z();
  double         angle     = dims.angle();

  Material mat = desc.material(detElem.materialStr());

  Tube             outerTubeShape(rMin, rInner + thickness, innerZ + thickness);
  Tube             innerTubeShape(0, rInner, innerZ);
  SubtractionSolid unchamferedShape(outerTubeShape, innerTubeShape);
  Cone             chamferShape(thickness, 0, rMin, 0, rMin + 2 * tan(angle) * thickness);
  SubtractionSolid detShape(unchamferedShape, chamferShape, Position(0, 0, innerZ + thickness));
  Volume           detVol(detName, detShape, mat);

  detVol.setVisAttributes(desc.visAttributes(detElem.visStr()));

  DetElement   det(detName, detID);
  Volume       motherVol = desc.pickMotherVolume(det);
  PlacedVolume detPV     = motherVol.placeVolume(detVol);
  det.setPlacement(detPV);
  return det;
}

// clang-format off
DECLARE_DETELEMENT(TestDetector, createDetector)
