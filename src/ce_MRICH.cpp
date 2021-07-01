
#include <XML/Helper.h>
#include "TMath.h"
#include "TString.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "GeometryHelpers.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;


void addModules(Volume &mother, xml::DetElement &detElem, Detector &desc, SensitiveDetector &sens);


// create the detector
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
    xml::DetElement detElem = handle;
    std::string detName = detElem.nameStr();
    int detID = detElem.id();

    DetElement det(detName, detID);
    xml::Component dims = detElem.dimensions();
    // xml::Component rads = detElem.child(_Unicode(radiator));

    auto rmin = dims.rmin();
    auto rmax = dims.rmax();
    auto length = dims.length();
    auto z0 = dims.z();

    auto gasMat = desc.material("AirOptical");

    // detector envelope
    Tube envShape(rmin, rmax, length/2., 0., 2*M_PI);
    Volume envVol("ce_MRICH_GVol", envShape, gasMat);
    envVol.setVisAttributes(desc.visAttributes(detElem.visStr()));

    // modules
    addModules(envVol, detElem, desc, sens);

    // place envelope
    Volume motherVol = desc.pickMotherVolume(det);
    PlacedVolume envPV = motherVol.placeVolume(envVol, Position(0, 0, z0));
    envPV.addPhysVolID("system", detID);
    det.setPlacement(envPV);
    return det;
}


void addModules(Volume &mother, xml::DetElement &detElem, Detector &desc, SensitiveDetector &sens)
{
    xml::Component dims = detElem.dimensions();
    xml::Component mods = detElem.child(_Unicode(modules));

    auto rmin = dims.rmin();
    auto rmax = dims.rmax();

    auto mThick = mods.attr<double>(_Unicode(thickness));
    auto mWidth = mods.attr<double>(_Unicode(width));
    auto mGap = mods.attr<double>(_Unicode(gap));

    auto modMat = desc.material(mods.materialStr());
    auto gasMat = desc.material("AirOptical");

    // single module
    Box mShape(mWidth/2., mWidth/2., mThick/2. - 0.1*mm);
    Volume mVol("ce_MRICH_mod_Solid", mShape, modMat);

    // a thin gas layer to detect optical photons
    Box modShape(mWidth/2., mWidth/2., mThick/2.);
    Volume modVol("ce_MRICH_mod_Solid_v", modShape, gasMat);
    // thin gas layer is on top (+z) of the material
    modVol.placeVolume(mVol, Position(0., 0., -0.1*mm));

    modVol.setVisAttributes(desc.visAttributes(mods.visStr()));
    sens.setType("photoncounter");
    modVol.setSensitiveDetector(sens);

    // place modules in the sectors (disk)
    auto points = athena::geo::fillSquares({0., 0.}, mWidth + mGap, rmin - mGap, rmax + mGap);

    // determine module direction, always facing z = 0
    double roty = dims.z() > 0. ? M_PI/2. : -M_PI/2.;
    int imod = 1;
    for (auto &p : points) {
        // operations are inversely ordered
        Transform3D tr = Translation3D(p.x(), p.y(), 0.)        // move to position
                       * RotationY(roty);                       // facing z = 0.
        auto modPV = mother.placeVolume(modVol, tr);
        modPV.addPhysVolID("sector", 1).addPhysVolID("module", imod ++);
    }
}

// clang-format off
DECLARE_DETELEMENT(refdet_ce_MRICH, createDetector)

