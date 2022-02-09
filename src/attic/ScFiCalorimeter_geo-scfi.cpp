#include <iostream>
#include <limits>
#include <cmath>
#include <tuple>
#include <algorithm>

#include <XML/Helper.h>
#include "DD4hep/DetFactoryHelper.h"
#include "GeometryHelpers.h"

using namespace dd4hep;
using Point =  ROOT::Math::XYPoint;
const double eps = std::numeric_limits<float>::epsilon();

std::tuple<Volume, Position> build_module(const Detector &desc, const xml::Component &mod_xml, SensitiveDetector &sens);

// helper function to get x, y, z if defined in a xml component
template<class XmlComp>
Position get_xml_xyz(const XmlComp &comp, dd4hep::xml::Strng_t name)
{
  Position pos(0., 0., 0.);
  if (comp.hasChild(name)) {
    auto child = comp.child(name);
    pos.SetX(dd4hep::getAttrOrDefault<double>(child, _Unicode(x), 0.));
    pos.SetY(dd4hep::getAttrOrDefault<double>(child, _Unicode(y), 0.));
    pos.SetZ(dd4hep::getAttrOrDefault<double>(child, _Unicode(z), 0.));
  }
  return pos;
}

// main
static Ref_t create_detector(Detector &desc, xml::Handle_t handle, SensitiveDetector sens)
{
  xml::DetElement detElem = handle;
  std::string detName = detElem.nameStr();
  int detID = detElem.id();
  DetElement det(detName, detID);
  sens.setType("calorimeter");
  auto dim = detElem.dimensions();
  auto length = dim.length();

  // build module
  auto [modVol, modSize] = build_module(desc, detElem.child(_Unicode(module)), sens);
  const xml::Component &glue_xml = detElem.child(_Unicode(glue));
  const auto glue_width = dd4hep::getAttrOrDefault<double>(glue_xml, _Unicode(width), 0.1*mm);
  const double glue_x = modSize.x()*2. + glue_width;
  const double glue_y = modSize.y()*2. + glue_width;
  auto glueMat = desc.material(glue_xml.attr<std::string>(_Unicode(material)));
  Box GlueHorizShape(glue_x/2., glue_width/2., modSize.z()/2.);
  Box GlueVertShape(glue_width/2., glue_y/2., modSize.z()/2.);
  UnionSolid GlueShape(GlueHorizShape, GlueVertShape);
  Volume modGlue("modGlue_vol", GlueShape, glueMat);

  Assembly blockVol("block_vol");
  int modid = 1;
  for (int ix = 0; ix < 2; ix++)
    for (int iy = 0; iy < 2; iy++) {
      auto modPV = blockVol.placeVolume(modVol,
          Position{(modSize.x()+glue_width)*(ix-0.5), (modSize.y()+glue_width)*(iy-0.5), 0});
      modPV.addPhysVolID("module", modid++);
    }
  blockVol.placeVolume(modGlue, Position{0, 0, 0});
  blockVol->Voxelize("");

  const xml::Component &block_xml = detElem.child(_Unicode(block));
  const auto block_x = block_xml.attr<double>(_Unicode(sizex));
  const auto block_y = block_xml.attr<double>(_Unicode(sizey));
  // calorimeter block z-offsets (as blocks are shorter than the volume length)
  const double block_offset = -length/2. + modSize.z()/2.;
  Box envShape(block_x*8./2., block_y*8./2., length/2.);
  Volume envVol(detName + "_envelope", envShape, desc.material("Vacuum"));
  envVol.setVisAttributes(desc.visAttributes(detElem.visStr()));
  int blockid = 1;
  for (int ix = 0; ix < 8; ix++)
    for (int iy = 0; iy < 8; iy++) {
      auto blockPV = envVol.placeVolume(blockVol,
          Position{block_x*(ix-3.5), block_y*(iy-3.5), block_offset});
      blockPV.addPhysVolID("block", blockid++);
    }

  desc.add(Constant(detName + "_NModules", std::to_string((blockid-1)*(modid-1)*4)));

  // detector position and rotation
  auto pos = get_xml_xyz(detElem, _Unicode(position));
  auto rot = get_xml_xyz(detElem, _Unicode(rotation));
  Volume motherVol = desc.pickMotherVolume(det);
  Transform3D tr = Translation3D(pos.x(), pos.y(), pos.z()) * RotationZYX(rot.z(), rot.y(), rot.x());
  PlacedVolume envPV = motherVol.placeVolume(envVol, tr);
  envPV.addPhysVolID("system", detID);
  det.setPlacement(envPV);
  return det;
}

// helper function to build module with scintillating fibers
std::tuple<Volume, Position> build_module(const Detector &desc, const xml::Component &mod_xml, SensitiveDetector &sens)
{
  auto sx = mod_xml.attr<double>(_Unicode(sizex));
  auto sy = mod_xml.attr<double>(_Unicode(sizey));
  auto sz = mod_xml.attr<double>(_Unicode(sizez));
  auto modMat = desc.material(mod_xml.attr<std::string>(_Unicode(material)));

  auto fiber_xml = mod_xml.child(_Unicode(fiber));
  auto fr = fiber_xml.attr<double>(_Unicode(radius));
  auto fdistx = fiber_xml.attr<double>(_Unicode(spacex));
  auto fdisty = fiber_xml.attr<double>(_Unicode(spacey));
  auto foffx = dd4hep::getAttrOrDefault<double>(fiber_xml, _Unicode(offsetx), 0.5*mm);
  auto foffy = dd4hep::getAttrOrDefault<double>(fiber_xml, _Unicode(offsety), 0.5*mm);
  auto fiberMat = desc.material(fiber_xml.attr<std::string>(_Unicode(material)));

  Box modShape(sx/2., sy/2., sz/2.);
  Volume modVol("module_vol", modShape, modMat);
  if (mod_xml.hasAttr(_Unicode(vis)))
    modVol.setVisAttributes(desc.visAttributes(mod_xml.attr<std::string>(_Unicode(vis))));

  Tube fiberShape(0., 0.94*fr, sz/2.);
  Volume fiberVol("~fiber_vol", fiberShape, fiberMat);
  Tube fiberCladShape(0.94*fr, fr, sz/2.);
  Volume fiberCladVol("fiberClad_vol", fiberCladShape, fiberMat);
  fiberVol.setSensitiveDetector(sens);

  // place the fibers
  int nx = int(sx / fdistx) + 1;
  int ny = int(sy / fdisty) + 1;
  for (int iy = 0; iy < ny; iy++) {
    double y = foffy + fdisty * iy;
    // about to touch the boundary
    if (sy - y < foffy - eps) break;
    for (int ix = 0; ix < nx; ix++) {
      double x = foffx + fdistx/2. * (iy%2) + fdistx * ix;
      // about to touch the boundary
      if (sx - x < foffx - eps) break;
      auto fiberPV = modVol.placeVolume(fiberVol, Position{-sx/2.+x, sy/2.-y, 0});
      modVol.placeVolume(fiberCladVol, Position{-sx/2.+x, sy/2.-y, 0});
      fiberPV.addPhysVolID("~fiber_x", ix + 1).addPhysVolID("~fiber_y", iy + 1);
    }
  }

  return std::make_tuple(modVol, Position{sx, sy, sz});
}

DECLARE_DETELEMENT(ScFiCalorimeter, create_detector)
