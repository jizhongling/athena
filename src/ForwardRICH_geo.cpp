//==========================================================================
//  Forward Ring Imaging Cherenkov Detector
//--------------------------------------------------------------------------
//
// Author: C. Peng (ANL)
// Date: 09/30/2020
//
//==========================================================================

#include <XML/Helper.h>
#include "TMath.h"
#include "TString.h"
#include "Math/Point2D.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

typedef ROOT::Math::XYPoint Point;


// check if a square in a ring
inline bool in_ring(const Point &pt, double side, double rmin, double rmax, double phmin, double phmax)
{
    if (pt.r() > rmax || pt.r() < rmin) {
        return false;
    }

    // check four corners
    std::vector<Point> pts {
        Point(pt.x() - side/2., pt.y() - side/2.),
        Point(pt.x() - side/2., pt.y() + side/2.),
        Point(pt.x() + side/2., pt.y() - side/2.),
        Point(pt.x() + side/2., pt.y() + side/2.),
    };
    for (auto &p : pts) {
        if (p.r() > rmax || p.r() < rmin || p.phi() > phmax || p.phi() < phmin) {
            return false;
        }
    }
    return true;
}

// check if a square is overlapped with the others
inline bool overlap(const Point &pt, double side, const std::vector<Point> &pts)
{
    for (auto &p : pts) {
        auto pn = (p - pt)/side;
        if ((std::abs(pn.x()) < 1. - 1e-6) && (std::abs(pn.y()) < 1. - 1e-6)) {
            return true;
        }
    }
    return false;
}

// a helper function to recursively fill square in a ring
void add_square(Point p, std::vector<Point> &res, double lside, double rmin, double rmax,
                double phmin, double phmax)
{
    // outside of the ring or overlapping
    if (!in_ring(p, lside, rmin, rmax, phmin, phmax) || overlap(p, lside, res)) {
        return;
    }

    res.emplace_back(p);

    // check adjacent squares
    add_square(Point(p.x() + lside, p.y()), res, lside, rmin, rmax, phmin, phmax);
    add_square(Point(p.x() - lside, p.y()), res, lside, rmin, rmax, phmin, phmax);
    add_square(Point(p.x(), p.y() + lside), res, lside, rmin, rmax, phmin, phmax);
    add_square(Point(p.x(), p.y() - lside), res, lside, rmin, rmax, phmin, phmax);
}

// fill squares
std::vector<Point> fill_squares(Point ref, double lside, double rmin, double rmax,
                                double phmin = 0., double phmax = 2.*M_PI)
{
    // start with a seed square and find one in the ring
    // move to center
    ref = ref - Point(int(ref.x()/lside)*lside, int(ref.y()/lside)*lside);

    auto find_seed = [] (const Point &ref, int n, double side, double rmin, double rmax, double phmin, double phmax) {
        for (int ix = -n; ix < n; ++ix) {
            for (int iy = -n; iy < n; ++iy) {
                Point pt(ref.x() + ix*side, ref.y() + iy*side);
                if (in_ring(pt, side, rmin, rmax, phmin, phmax)) {
                    return pt;
                }
            }
        }
        return ref;
    };

    std::vector<Point> res;
    ref = find_seed(ref, int(rmax/lside) + 2, lside, rmin, rmax, phmin, phmax);
    add_square(ref, res, lside, rmin, rmax, phmin, phmax);
    return res;
}

// create the detector
static Ref_t createDetector(Detector& desc, xml::Handle_t handle, SensitiveDetector sens)
{
    xml::DetElement detElem = handle;

    std::string detName = detElem.nameStr();
    int detID = detElem.id();

    DetElement det(detName, detID);
    xml::Component dims = detElem.dimensions();
    xml::Component rads = detElem.child(_Unicode(radiator));
    xml::Component mir = detElem.child(_Unicode(mirror));
    xml::Component mcp = detElem.child(_Unicode(mcppmt));

    // dimensions
    double z0 = dims.z0();
    double length = dims.length();
    double rmin = dims.rmin();
    double rmax1 = dims.attr<double>(_Unicode(rmax1));
    double rmax2 = dims.attr<double>(_Unicode(rmax2));

    // mirror setting
    auto mThick = mir.thickness();
    auto mirZ = mir.attr<double>(_Unicode(zdiff));

    // mcppmt setting
    auto pRmin = mcp.rmin();
    auto pRmax = mcp.rmax();
    auto pThick = mcp.thickness();
    auto pSize = mcp.attr<double>(_Unicode(module_size));
    auto pGap = mcp.attr<double>(_Unicode(module_gap));
    auto pTol = mcp.attr<double>(_Unicode(rtol));
    auto pZ = mcp.attr<double>(_Unicode(zdiff));

    // materials
    auto mirMat = desc.material(mir.materialStr());
    auto gasMat = desc.material(rads.materialStr());
    auto mcpMat = desc.material(mcp.materialStr());

    // constants
    auto richCenterAngle = std::atan((rmin + (rmax2 - rmin)/2.)/mirZ);
    //std::cout << richCenterAngle*180./M_PI << std::endl;

    // an envelope for the detector
    // use a complicated shape to avoid conflict with the other parts
    // cone for radiator and the first set of mirrors
    double halfLength = length/2.;
    Cone env1(halfLength, rmin, rmax1, rmin, rmax2);
    // envelope for detection plane
    // Cone env2(halfLength - pZ/2., rmin, pRmax, rmin, rmax2);
    Tube env2(rmin, pRmax + pTol + pGap + 1.0*cm, (length - pZ)/2., 0., 2*M_PI);

    UnionSolid envShape(env1, env2, Position(0., 0., pZ));

    Volume envVol(detName + "_envelope", envShape, gasMat);
    envVol.setVisAttributes(desc.visAttributes(detElem.visStr()));

    // ---------------
    // spherical mirrors inside it
    int ilayer = 1;

    // optical surface
    OpticalSurfaceManager surfMgr = desc.surfaceManager();
    OpticalSurface mirSurf  = surfMgr.opticalSurface("MirrorOpticalSurface");
    // mirror slices
    int imod = 1;
    for (xml::Collection_t sl(mir, _Unicode(slice)); sl; ++sl, ++imod) {
        auto focus = sl.attr<double>(_Unicode(focus));
        auto wphi = sl.attr<double>(_Unicode(phiw));
        auto rotZ = sl.attr<double>(_Unicode(rotz));
        auto mRmin = sl.attr<double>(_Unicode(rmin));
        auto mRmax = sl.attr<double>(_Unicode(rmax));
        double curve = 0.;
        if (sl.hasAttr(_Unicode(curve))) {
            curve = sl.attr<double>(_Unicode(curve));
        }
        // geometry of mirror slice
        PlacedVolume mirPV;
        Volume mirVol(Form("mirror_v_dummy%d", imod));
        mirVol.setMaterial(mirMat);
        mirVol.setVisAttributes(desc.visAttributes(mir.visStr()));
        // spherical mirror
        if (curve > 0.) {
            // somehow geant4 does not support -wphi/2. to wphi/2., so additonal rotation in Z
            double mTheta1 = std::asin(mRmin/curve);
            double mTheta2 = std::asin(mRmax/curve);
            double rotY = -std::asin(focus/curve);
            mirVol.setSolid(Sphere(curve, curve + mThick, mTheta1, mTheta2, 0., wphi));
            // action is in a reverse order
            Transform3D tr = Translation3D(0., 0., mirZ - halfLength)   // move for z position
                           * RotationZ(rotZ)                            // rotate phi angle
                           * RotationY(rotY)                            // rotate for focus point
                           * RotationX(180*degree)
                           * Translation3D(0., 0., -curve)              // move spherical shell to origin
                           * RotationZ(-wphi/2.);                       // center phi angle to 0. (-wphi/2., wphi/2.)
            mirPV = envVol.placeVolume(mirVol, tr);
        // plane mirror
        } else {
            mirVol.setSolid(Tube(mRmin, mRmax, mThick/2.0, 0., wphi));
            Transform3D tr = Translation3D(0., 0., mirZ - halfLength)   // move for z position
                           * RotationZ(rotZ)                            // rotate phi angle
                           * RotationZ(-wphi/2.);                       // center phi angle to 0. (-wphi/2., wphi/2.)
            mirPV = envVol.placeVolume(mirVol, tr);
        }
        mirPV.addPhysVolID("layer", ilayer).addPhysVolID("module", imod);
        DetElement mirDE(det, Form("Mirror_DE%d", imod), imod);
        mirDE.setPlacement(mirPV);
        SkinSurface mirSurfBorder(desc, mirDE, Form("RICHmirror%d", imod), mirSurf, mirVol);
        mirSurfBorder.isValid();
    }
    ilayer++;

    // ---------------
    // photo-detector unit
    // Fill the photo-detection plane with square shape MCP-PMTs
    Box mcpShape1(pSize/2.0, pSize/2.0, pThick/2.0);
    Volume mcpVol1("mcppmt_v_material", mcpShape1, mcpMat);

    // a thin layer of cherenkov gas for accepting optical photons
    Box mcpShape(pSize/2.0, pSize/2.0, pThick/2.0 + 0.1*mm);
    Volume mcpVol("mcppmt_v", mcpShape, gasMat);
    mcpVol.placeVolume(mcpVol1, Position(0., 0., -0.1*mm));

    mcpVol.setVisAttributes(desc.visAttributes(mcp.visStr()));
    sens.setType("photoncounter");
    mcpVol.setSensitiveDetector(sens);

    // photo-detector plane envelope
    for (size_t ipd = 0; ipd < 6; ++ipd) {
        double phmin = -M_PI/6.;
        double phmax = M_PI/6.;
        Tube pdEnvShape(pRmin - pTol - pGap, pRmax + pTol + pGap, pThick/2.0 + 0.1*cm, phmin, phmax);
        Volume pdVol("pd_envelope", pdEnvShape, desc.material("AirOptical"));
        auto points = fill_squares(Point(0., 0.), pSize + pGap, pRmin - pTol - pGap, pRmax + pTol + pGap, phmin, phmax);
        for (size_t i = 0; i < points.size(); ++i) {
            auto pt = points[i];
            auto mcpPV = pdVol.placeVolume(mcpVol, Position(pt.x(), pt.y(), 0.));
            mcpPV.addPhysVolID("layer", ilayer).addPhysVolID("module", i + 1);
            DetElement mcpDE(det, Form("MCPPMT_DE%d_%d", ipd + 1, i + 1), i + 1);
            mcpDE.setPlacement(mcpPV);
        }
        Transform3D tr = Translation3D(0., 0., -halfLength + pZ + pThick/2.0)   // move for z position
                        * RotationZ(ipd*M_PI/3.)        // rotate phi angle
                        * RotationY(-richCenterAngle);  // rotate to perpendicular position
        auto pdPV = envVol.placeVolume(pdVol, tr);
        pdPV.addPhysVolID("layer", ilayer).addPhysVolID("piece", ipd + 1);
    }
    Volume motherVol = desc.pickMotherVolume(det);
    PlacedVolume envPV = motherVol.placeVolume(envVol, Position(0, 0, z0 + halfLength));
    envPV.addPhysVolID("system", detID);
    det.setPlacement(envPV);

    return det;
}
//@}

// clang-format off
DECLARE_DETELEMENT(refdet_ForwardRICH, createDetector)

