/** \addtogroup VertexTracker Vertex Trackers
 * \brief Type: **SiVertexBarrel**.
 * \author W. Armstrong
 * \ingroup trackers
 *
 *
 * \code
 * \endcode
 *
 * @{
 */
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/Shapes.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "XML/Layering.h"
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Definitions/Units.hpp"


using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector& description, xml_h e, SensitiveDetector sens) {
  typedef vector<PlacedVolume> Placements;
  xml_det_t                    x_det    = e;
  Material                     air      = description.air();
  int                          det_id   = x_det.id();
  string                       det_name = x_det.nameStr();
  DetElement                   sdet(det_name, det_id);
  //Assembly                     assembly(det_name);
  map<string, Volume>          volumes;
  map<string, Placements>      sensitives;
  map<string, std::vector<VolPlane>>        volplane_surfaces;
  map<string, xml_h>      xmleles;
  PlacedVolume                 pv;
  dd4hep::xml::Dimension dimensions(x_det.dimensions());

  map<string, std::array<double, 2>> module_thicknesses;

  Acts::ActsExtension* detWorldExt = new Acts::ActsExtension();
  detWorldExt->addType("barrel", "detector");
  sdet.addExtension<Acts::ActsExtension>(detWorldExt);
  Tube topVolumeShape(dimensions.rmin(), dimensions.rmax(), dimensions.length() * 0.5);
  Volume assembly(det_name,topVolumeShape,air);


//  The Cold Plate is approximately 30 mm wide and is based on the same carbon-ply layup
//as for the IB Stave. Two pipes with an inner diameter of 2.67 mm and a wall thickness
//of 64 μm have been used. The two pipes are interconnected at one end of the Cold Plate
//providing a loop, whose inlet and outlet are on the same side and correspond to the
//
//requirements have led to an equilateral section of the frame with a 42 mm wide side, that
//provides almost the same rigidity for all the possible Stave positions.
//
//                     module mat                um
//
//Aluminium                                      50
//Polyimide                                      100
//Carbon fibre                                   120
//Silicon                                        50
//Eccobond-45                                    100
//
//Metal layers                 Aluminium        200
//Insulating layers            Polyimide        200
//Glue Cooling tube wall       Eccobond-45      100
//
//Carbon fleece            40
//Carbon paper             30
//Polyimide                64
//Water                    
//Carbon fibre             120
//Eccobond-45              100



  sens.setType("tracker");

  // loop over the modules
  for (xml_coll_t mi(x_det, _U(module)); mi; ++mi) {
    xml_comp_t x_mod = mi;
    xml_comp_t m_env = x_mod.child(_U(frame));
    string     m_nam = x_mod.nameStr();
    xmleles[m_nam]  = x_mod;

    // triangular volume envelope
    double frame_thickness = m_env.thickness();
    double frame_width     = m_env.width();
    double frame_height    = getAttrOrDefault<double>(m_env, _U(height), 5.0 * mm);
    double tanth           = frame_height/(frame_width/2.0);
    double frame_height2   = frame_height-frame_thickness-frame_thickness/tanth;
    double frame_width2    = 2.0*frame_height2/tanth;

    Trd1 moduleframe_part1(frame_width / 2, 0.001 * mm, m_env.length() / 2,
                           frame_height / 2);
    Trd1 moduleframe_part2(frame_width2/2, 0.001 * mm,
                           m_env.length() / 2 + 0.01 * mm, frame_height2/2);
    SubtractionSolid moduleframe(moduleframe_part1, moduleframe_part2,Position(0.0,frame_thickness,0.0));
    Volume v_module(m_nam+"_vol", moduleframe, description.material(m_env.materialStr()));
    v_module.setVisAttributes(description, m_env.visStr());

    int ncomponents   = 0;
    int sensor_number = 1;

    if (volumes.find(m_nam) != volumes.end()) {
      printout(ERROR, "SiTrackerBarrel", "Logics error in building modules.");
      throw runtime_error("Logics error in building modules.");
    }
    double total_thickness = 0;
    xml_coll_t ci(x_mod, _U(module_component));
    for (ci.reset(), total_thickness = 0.0; ci; ++ci) {
      total_thickness += xml_comp_t(ci).thickness();
    }
    // module assembly
    Assembly m_vol( m_nam );
    m_vol.placeVolume(v_module, Position(0.0,0.0,frame_height/2+total_thickness/2.0));
    volumes[m_nam] = m_vol;
    m_vol.setVisAttributes(description.visAttributes(x_mod.visStr()));

    double thickness_so_far = 0.0;
    double thickness_sum = -total_thickness/2.0;
    for (xml_coll_t ci(x_mod, _U(module_component)); ci; ++ci, ++ncomponents) {
      xml_comp_t x_comp = ci;
      xml_comp_t x_pos  = x_comp.position(false);
      xml_comp_t x_rot  = x_comp.rotation(false);
      string     c_nam  = _toString(ncomponents, "component%d");
      Box        c_box(x_comp.width() / 2, x_comp.length() / 2, x_comp.thickness() / 2);
      Volume     c_vol(c_nam, c_box, description.material(x_comp.materialStr()));

      if (x_pos && x_rot) {
        Position c_pos(x_pos.x(0), x_pos.y(0), x_pos.z(0));
        RotationZYX c_rot(x_rot.z(0), x_rot.y(0), x_rot.x(0));
        pv = m_vol.placeVolume(c_vol, Transform3D(c_rot, c_pos));
      } else if (x_rot) {
        Position c_pos(0, 0, thickness_sum + x_comp.thickness() / 2.0);
        pv = m_vol.placeVolume(c_vol, Transform3D(RotationZYX(x_rot.z(0), x_rot.y(0), x_rot.x(0)),c_pos));
      } else if (x_pos) {
        pv = m_vol.placeVolume(c_vol, Position(x_pos.x(0), x_pos.y(0), x_pos.z(0)));
      } else {
        pv = m_vol.placeVolume(c_vol, Position(0,0,thickness_sum+x_comp.thickness()/2.0));
      }
      c_vol.setRegion(description, x_comp.regionStr());
      c_vol.setLimitSet(description, x_comp.limitsStr());
      c_vol.setVisAttributes(description, x_comp.visStr());
      if (x_comp.isSensitive()) {
        pv.addPhysVolID("sensor", sensor_number++);
        c_vol.setSensitiveDetector(sens);
        sensitives[m_nam].push_back(pv);
        module_thicknesses[m_nam] = {thickness_so_far + x_comp.thickness()/2.0, total_thickness-thickness_so_far - x_comp.thickness()/2.0};


          // -------- create a measurement plane for the tracking surface attched to the sensitive volume -----
          Vector3D u( 0. , 1. , 0. ) ;
          Vector3D v( 0. , 0. , 1. ) ;
          Vector3D n( 1. , 0. , 0. ) ;
          //    Vector3D o( 0. , 0. , 0. ) ;

          // compute the inner and outer thicknesses that need to be assigned to the tracking surface
          // depending on wether the support is above or below the sensor
          double inner_thickness = module_thicknesses[m_nam][0];
          double outer_thickness = module_thicknesses[m_nam][1];

          SurfaceType type( SurfaceType::Sensitive ) ;

          //if( isStripDetector )
          //  type.setProperty( SurfaceType::Measurement1D , true ) ;

          VolPlane surf( c_vol , type , inner_thickness , outer_thickness , u,v,n ) ; //,o ) ;
          volplane_surfaces[m_nam].push_back(surf);

    //--------------------------------------------


      }
      thickness_sum += x_comp.thickness();
      thickness_so_far += x_comp.thickness();
    }
  }

  // now build the layers
  for (xml_coll_t li(x_det, _U(layer)); li; ++li) {
    xml_comp_t x_layer  = li;
    xml_comp_t x_barrel = x_layer.child(_U(barrel_envelope));
    xml_comp_t x_layout = x_layer.child(_U(rphi_layout));
    xml_comp_t z_layout = x_layer.child(_U(z_layout)); // Get the <z_layout> element.
    int        lay_id   = x_layer.id();
    string     m_nam    = x_layer.moduleStr();
    string     lay_nam  = _toString(x_layer.id(), "layer%d");
    Tube       lay_tub(x_barrel.inner_r(), x_barrel.outer_r(), x_barrel.z_length() / 2.0);
    Volume     lay_vol(lay_nam, lay_tub, air); // Create the layer envelope volume.
    lay_vol.setVisAttributes(description.visAttributes(x_layer.visStr()));

    double      phi0       = x_layout.phi0();     // Starting phi of first module.
    double      phi_tilt   = x_layout.phi_tilt(); // Phi tilt of a module.
    double      rc         = x_layout.rc();       // Radius of the module center.
    int         nphi       = x_layout.nphi();     // Number of modules in phi.
    double      rphi_dr    = x_layout.dr();       // The delta radius of every other module.
    double      phi_incr   = (M_PI * 2) / nphi;   // Phi increment for one module.
    double      phic       = phi0;                // Phi of the module center.
    double      z0         = z_layout.z0();       // Z position of first module in phi.
    double      nz         = z_layout.nz();       // Number of modules to place in z.
    double      z_dr       = z_layout.dr();       // Radial displacement parameter, of every other module.

    Volume      module_env = volumes[m_nam];
    DetElement  lay_elt(sdet, _toString(x_layer.id(), "layer%d"), lay_id);
    Placements& sensVols = sensitives[m_nam];

    // the local coordinate systems of modules in dd4hep and acts differ
    // see http://acts.web.cern.ch/ACTS/latest/doc/group__DD4hepPlugins.html
    Acts::ActsExtension* layerExtension = new Acts::ActsExtension();
    layerExtension->addType("sensitive cylinder", "layer");
    //layerExtension->addValue(0, "r_min", "envelope");
    //layerExtension->addValue(0, "r_max", "envelope");
    //layerExtension->addValue(0, "z_min", "envelope");
    //layerExtension->addValue(0, "z_max", "envelope");
    //layerExtension->addType("axes", "definitions", "XzY");
    lay_elt.addExtension<Acts::ActsExtension>(layerExtension);

    // Z increment for module placement along Z axis.
    // Adjust for z0 at center of module rather than
    // the end of cylindrical envelope.
    double z_incr = nz > 1 ? (2.0 * z0) / (nz - 1) : 0.0;
    // Starting z for module placement along Z axis.
    double module_z = -z0;
    int module = 1;

    // Loop over the number of modules in phi.
    for (int ii = 0; ii < nphi; ii++) {
      double dx = z_dr * std::cos(phic + phi_tilt);  // Delta x of module position.
      double dy = z_dr * std::sin(phic + phi_tilt);  // Delta y of module position.
      double x = rc * std::cos(phic);                // Basic x module position.
      double y = rc * std::sin(phic);                // Basic y module position.

      // Loop over the number of modules in z.
      for (int j = 0; j < nz; j++) {
        string module_name = _toString(module, "module%d");
        DetElement mod_elt(lay_elt, module_name, module);

        Transform3D tr(RotationZYX(0, ((M_PI / 2) - phic - phi_tilt), -M_PI / 2),
                       Position(x, y, module_z));

        pv = lay_vol.placeVolume(module_env, tr);
        pv.addPhysVolID("module", module);
        mod_elt.setPlacement(pv);
        for (size_t ic = 0; ic < sensVols.size(); ++ic) {
          PlacedVolume sens_pv = sensVols[ic];
          DetElement comp_de(mod_elt, std::string("de_") + sens_pv.volume().name(), module);
          comp_de.setPlacement(sens_pv);
          Acts::ActsExtension* sensorExtension = new Acts::ActsExtension();
          //sensorExtension->addType("sensor", "detector");
          comp_de.addExtension<Acts::ActsExtension>(sensorExtension);
          //comp_de.setAttributes(description, sens_pv.volume(), x_layer.regionStr(), x_layer.limitsStr(),
          //                      xml_det_t(xmleles[m_nam]).visStr());
          //

      volSurfaceList( comp_de )->push_back( volplane_surfaces[m_nam][ic] ) ;


        }

        /// Increase counters etc.
        module++;
        // Adjust the x and y coordinates of the module.
        x += dx;
        y += dy;
        // Flip sign of x and y adjustments.
        dx *= -1;
        dy *= -1;
        // Add z increment to get next z placement pos.
        module_z += z_incr;
      }
      phic += phi_incr;  // Increment the phi placement of module.
      rc += rphi_dr;     // Increment the center radius according to dr parameter.
      rphi_dr *= -1;     // Flip sign of dr parameter.
      module_z = -z0;    // Reset the Z placement parameter for module.
    }
    // Create the PhysicalVolume for the layer.
    pv = assembly.placeVolume(lay_vol);  // Place layer in mother
    pv.addPhysVolID("layer", lay_id);    // Set the layer ID.
    lay_elt.setAttributes(description, lay_vol, x_layer.regionStr(), x_layer.limitsStr(),
                          x_layer.visStr());
    lay_elt.setPlacement(pv);
  }
  sdet.setAttributes(description, assembly, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  assembly.setVisAttributes(description.invisible());
  pv = description.pickMotherVolume(sdet).placeVolume(assembly);
  pv.addPhysVolID("system", det_id);  // Set the subdetector system ID.
  pv.addPhysVolID("barrel", 1);       // Flag this as a barrel subdetector.
  sdet.setPlacement(pv);
  return sdet;
}

//@}
// clang-format off
DECLARE_DETELEMENT(BarrelTrackerWithFrame, create_detector)
DECLARE_DETELEMENT(athena_TrackerBarrel, create_detector)
DECLARE_DETELEMENT(athena_VertexBarrel,  create_detector)
