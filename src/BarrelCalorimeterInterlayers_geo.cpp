// Detector plugin to support a hybrid central barrel calorimeter
// The detector consists of interlayers of Pb/ScFi (segmentation in global r, phi) and W/Si (segmentation in local x, y)
// Assembly is used as the envelope so two different detectors can be interlayered with each other
//
//
// Implementation of the Sci Fiber geometry: M. Żurek 06/19/2021
// Support interlayers between multiple detectors: C. Peng 07/09/2021


#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "Math/Point2D.h"
#include "TGeoPolygon.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

typedef ROOT::Math::XYPoint Point;
// headers for helper functions, defined in BarrelCalorimeterHybrid_geo
vector<vector<Point>> fiberPositions(double radius, double x_spacing, double z_spacing,
                                     double x, double z, double phi, double spacing_tol = 1e-2);
std::pair<int, int> getNdivisions(double x, double z, double dx, double dz);
vector<tuple<int, Point, Point, Point, Point>> gridPoints(int div_x, int div_z, double x, double z, double phi);
void buildFibers(Detector& desc, SensitiveDetector sens, Volume &s_vol, xml_comp_t x_fiber,
                 std::tuple<double, double, double, double> dimensions);


// barrel ecal layers contained in an assembly
static Ref_t create_detector(Detector& desc, xml_h e, SensitiveDetector sens)  {
  Layering      layering (e);
  xml_det_t     x_det     = e;
  Material      air       = desc.air();
  int           det_id    = x_det.id();
  string        det_name  = x_det.nameStr();
  xml_comp_t    x_staves  = x_det.staves();
  double        offset    = x_det.attr<double>(_Unicode(offset));
  xml_comp_t    x_dim     = x_det.dimensions();
  int           nsides    = x_dim.numsides();
  double        inner_r   = x_dim.rmin();
  double        dphi      = (2*M_PI/nsides);
  double        hphi      = dphi/2;

  DetElement    sdet      (det_name, det_id);
  Volume        motherVol = desc.pickMotherVolume(sdet);

  Assembly      envelope  (det_name);
  Transform3D   tr        = Translation3D(0, 0, offset) * RotationZ(hphi);
  PlacedVolume  env_phv   = motherVol.placeVolume(envelope, tr);
  sens.setType("calorimeter");

  env_phv.addPhysVolID("system",det_id);
  sdet.setPlacement(env_phv);

  // build a single stave
  DetElement    stave_det("stave0", det_id);
  Assembly      mod_vol("stave");

  // keep tracking of the total thickness
  double l_pos_z = inner_r;
  { // =====  buildBarrelStave(desc, sens, module_volume) =====
    // Parameters for computing the layer X dimension:
    double tan_hphi = std::tan(hphi);
    double l_dim_y  = x_dim.z()/2.;

    // Loop over the sets of layer elements in the detector.
    int l_num = 1;
    for(xml_coll_t li(x_det, _U(layer)); li; ++li)  {
      xml_comp_t x_layer = li;
      int repeat = x_layer.repeat();
      double l_space_between = dd4hep::getAttrOrDefault(x_layer, _Unicode(space_between), 0.);
      double l_space_before = dd4hep::getAttrOrDefault(x_layer, _Unicode(space_before), 0.);
      l_pos_z += l_space_before;
      // Loop over number of repeats for this layer.
      for (int j = 0; j < repeat; j++)    {
        string l_name = Form("layer%d", l_num);
        double l_thickness = layering.layer(l_num - 1)->thickness();  // Layer's thickness.
        double l_dim_x = tan_hphi* l_pos_z;
        l_pos_z += l_thickness;

        Position   l_pos(0, 0, l_pos_z - l_thickness/2.);      // Position of the layer.
	    double l_trd_x1 = l_dim_x;
	    double l_trd_x2 = l_dim_x + l_thickness*tan_hphi;
	    double l_trd_y1 = l_dim_y;
	    double l_trd_y2 = l_trd_y1;
	    double l_trd_z  = l_thickness/2;
        Trapezoid  l_shape(l_trd_x1, l_trd_x2, l_trd_y1, l_trd_y2, l_trd_z);
        Volume     l_vol(l_name, l_shape, air);
        DetElement layer(stave_det, l_name, det_id);

        // Loop over the sublayers or slices for this layer.
        int s_num = 1;
        double s_pos_z = -(l_thickness / 2.);
        for(xml_coll_t si(x_layer,_U(slice)); si; ++si)  {
          xml_comp_t x_slice = si;
          string     s_name  = Form("slice%d", s_num);
          double     s_thick = x_slice.thickness();
	      double s_trd_x1 = l_dim_x + (s_pos_z + l_thickness/2)*tan_hphi;
	      double s_trd_x2 = l_dim_x + (s_pos_z + l_thickness/2 + s_thick)*tan_hphi;
	      double s_trd_y1 = l_trd_y1;
	      double s_trd_y2 = s_trd_y1;
	      double s_trd_z  = s_thick/2.;
          Trapezoid  s_shape(s_trd_x1, s_trd_x2, s_trd_y1, s_trd_y2, s_trd_z);
          Volume     s_vol(s_name, s_shape, desc.material(x_slice.materialStr()));
          DetElement slice(layer, s_name, det_id);

          // build fibers
          if (x_slice.hasChild(_Unicode(fiber))) {
              buildFibers(desc, sens, s_vol, x_slice.child(_Unicode(fiber)), {s_trd_x1, s_thick, l_dim_y, hphi});
          }

          if ( x_slice.isSensitive() ) {
            s_vol.setSensitiveDetector(sens);
          }
          s_vol.setAttributes(desc, x_slice.regionStr(), x_slice.limitsStr(), x_slice.visStr());

          // Slice placement.
          PlacedVolume slice_phv = l_vol.placeVolume(s_vol, Position(0, 0, s_pos_z + s_thick/2));
          slice_phv.addPhysVolID("slice", s_num);
          slice.setPlacement(slice_phv);
          // Increment Z position of slice.
          s_pos_z += s_thick;
          ++s_num;
        }

        // Set region, limitset, and vis of layer.
        l_vol.setAttributes(desc, x_layer.regionStr(), x_layer.limitsStr(), x_layer.visStr());

        PlacedVolume layer_phv = mod_vol.placeVolume(l_vol, l_pos);
        layer_phv.addPhysVolID("layer", l_num);
        layer.setPlacement(layer_phv);
        // Increment to next layer Z position. Do not add space_between for the last layer
        if (j < repeat - 1) {
          l_pos_z += l_space_between;
        }
        ++l_num;
      }
    }
  }
  // Phi start for a stave.
  double phi = M_PI / nsides;
  // Create nsides staves.
  for (int i = 0; i < nsides; i++, phi -= dphi)      { // i is module number
    // Compute the stave position
    Transform3D tr(RotationZYX(0, phi, M_PI*0.5), Translation3D(0, 0, 0));
    PlacedVolume pv = envelope.placeVolume(mod_vol, tr);
    pv.addPhysVolID("module", i + 1);
    DetElement sd = (i == 0) ? stave_det : stave_det.clone(Form("stave%d", i));
    sd.setPlacement(pv);
    sdet.add(sd);
  }

  Solid  support_frame_s;
  // optional stave support
  if (x_staves.hasChild("support")) {
    xml_comp_t  x_support           = x_staves.child(_U(support));
    double      support_thickness   = getAttrOrDefault(x_support, _U(thickness), 5.0 * cm);
    double      trd_x1_support      = (2 * std::tan(hphi) * l_pos_z + support_thickness)/2;
    // is the support on the inside surface?
    bool        is_inside_support   = getAttrOrDefault<bool>(x_support, _Unicode(inside), true);
    double      trd_x1              = std::tan(hphi) * inner_r;
    double      trd_x2              = std::tan(hphi) * (l_pos_z + support_thickness);
    double      trd_y1              = x_dim.z()/2.;

    // number of "beams" running the length of the stave.
    int    n_beams          = getAttrOrDefault<int>(x_support, _Unicode(n_beams), 3);
    double beam_thickness   = support_thickness / 4.0; // maybe a parameter later...
    trd_x1_support          = (2 * std::tan(hphi) * (l_pos_z + beam_thickness)) / 2.;
    double grid_size        = getAttrOrDefault(x_support, _Unicode(grid_size), 25.0 * cm);
    double beam_width       = 2.0 * trd_x1_support / (n_beams + 1); // quick hack to make some gap between T beams

    double cross_beam_thickness = support_thickness/4.0;
    //double trd_x1_support     = (2 * std::tan(hphi) * (inner_r + beam_thickness)) / 2.;
    double trd_x2_support       = trd_x2;

    int n_cross_supports = std::floor((trd_y1-cross_beam_thickness)/grid_size);

    Box        beam_vert_s(beam_thickness / 2.0 , trd_y1, support_thickness / 2.0 );
    Box        beam_hori_s(beam_width / 2.0, trd_y1, beam_thickness / 2.0);
    UnionSolid T_beam_s(beam_vert_s, beam_hori_s, Position(0, 0, -support_thickness / 2.0 + beam_thickness / 2.0));

    // cross supports
    Trapezoid  trd_support(trd_x1_support,trd_x2_support,
                           beam_thickness / 2.0, beam_thickness / 2.0,
                           support_thickness / 2.0 - cross_beam_thickness/2.0);
    UnionSolid support_array_start_s(T_beam_s,trd_support,Position(0,0,cross_beam_thickness/2.0));
    for (int isup = 0; isup < n_cross_supports; isup++) {
      support_array_start_s = UnionSolid(support_array_start_s, trd_support,
                                         Position(0, -1.0 * isup * grid_size, cross_beam_thickness/2.0));
      support_array_start_s = UnionSolid(support_array_start_s, trd_support,
                                         Position(0, 1.0 * isup * grid_size, cross_beam_thickness/2.0));
    }
    support_array_start_s =
        UnionSolid(support_array_start_s, beam_hori_s,
                   Position(-1.8 * 0.5*(trd_x1+trd_x2_support) / n_beams, 0, -support_thickness / 2.0 + beam_thickness / 2.0));
    support_array_start_s =
        UnionSolid(support_array_start_s, beam_hori_s,
                   Position(1.8 * 0.5*(trd_x1+trd_x2_support) / n_beams, 0, -support_thickness / 2.0 + beam_thickness / 2.0));
    support_array_start_s =
        UnionSolid(support_array_start_s, beam_vert_s, Position(-1.8 * 0.5*(trd_x1+trd_x2_support) / n_beams, 0, 0));
    support_array_start_s =
        UnionSolid(support_array_start_s, beam_vert_s, Position(1.8 * 0.5*(trd_x1+trd_x2_support) / n_beams, 0, 0));

    support_frame_s = support_array_start_s;

    Material support_mat = desc.material(x_support.materialStr());
    Volume   support_vol("support_frame_v", support_frame_s, support_mat);
    support_vol.setVisAttributes(desc,x_support.visStr());

    // figure out how to best place
    auto pv = mod_vol.placeVolume(support_vol, Position(0.0, 0.0, l_pos_z + support_thickness / 2.0));
  }

  //l_pos_z += support_thickness;
  // Set envelope volume attributes.
  envelope.setAttributes(desc, x_det.regionStr(), x_det.limitsStr(), x_det.visStr());
  return sdet;
}

void buildFibers(Detector& desc, SensitiveDetector sens, Volume &s_vol, xml_comp_t x_fiber,
                 std::tuple<double, double, double, double> dimensions)
{
  auto [s_trd_x1, s_thick, s_length, hphi] = dimensions;
  double f_radius = getAttrOrDefault(x_fiber, _U(radius), 0.1 * cm);
  double f_spacing_x = getAttrOrDefault(x_fiber, _Unicode(spacing_x), 0.122 * cm);
  double f_spacing_z = getAttrOrDefault(x_fiber, _Unicode(spacing_z), 0.134 * cm);
  std::string f_id_grid = getAttrOrDefault(x_fiber, _Unicode(identifier_grid), "grid");
  std::string f_id_fiber = getAttrOrDefault(x_fiber, _Unicode(identifier_fiber), "fiber");

  // Set up the readout grid for the fiber layers
  // Trapezoid is divided into segments with equal dz and equal number of divisions in x
  // Every segment is a polygon that can be attached later to the lightguide
  // The grid size is assumed to be ~2x2 cm (starting values). This is to be larger than
  // SiPM chip (for GlueX 13mmx13mm: 4x4 grid 3mmx3mm with 3600 50×50 μm pixels each)
  // See, e.g., https://arxiv.org/abs/1801.03088 Fig. 2d

  // Calculate number of divisions
  auto grid_div = getNdivisions(s_trd_x1, s_thick, 2.0*cm, 2.0*cm);
  // Calculate polygonal grid coordinates (vertices)
  auto grid_vtx = gridPoints(grid_div.first, grid_div.second, s_trd_x1, s_thick, hphi);
  Tube f_tube(0, f_radius, s_length);
  Volume f_vol("fiber_vol", f_tube, desc.material(x_fiber.materialStr()));

  vector<int> f_id_count(grid_div.first*grid_div.second, 0);
  auto f_pos = fiberPositions(f_radius, f_spacing_x, f_spacing_z, s_trd_x1, s_thick, hphi);
  // std::cout << f_pos.size() << " lines, ~" << f_pos.front().size() << " fibers each line" << std::endl;
  for (size_t il = 0; il < f_pos.size(); ++il) {
    auto &line = f_pos[il];
    if (line.empty()) {
      continue;
    }
    double l_pos_y = line.front().y();
    // use assembly as intermediate volume container to reduce number of daughter volumes
    Assembly lfibers(Form("fiber_array_line_%d", il));
    for (auto &p : line) {
      int f_grid_id = -1;
      int f_id = -1;
      // Check to which grid fiber belongs to
      for (auto &poly_vtx : grid_vtx) {
        if (p.y() != l_pos_y) {
          std::cerr << Form("Expected the same y position from a same line: %.2f, but got %.f", l_pos_y, p.y())
                    << std::endl;
          continue;
        }
        auto [grid_id, vtx_a, vtx_b, vtx_c, vtx_d] = poly_vtx;
        double poly_x[4] = {vtx_a.x(), vtx_b.x(), vtx_c.x(), vtx_d.x()};
        double poly_y[4] = {vtx_a.y(), vtx_b.y(), vtx_c.y(), vtx_d.y()};
        double f_xy[2] = {p.x(), p.y()};

        TGeoPolygon poly(4);
        poly.SetXY(poly_x, poly_y);
        poly.FinishPolygon();

        if(poly.Contains(f_xy)) {
          f_grid_id = grid_id;
          f_id = f_id_count[grid_id];
          f_id_count[grid_id]++;
        }
      }

      if ( x_fiber.isSensitive() ) {
        f_vol.setSensitiveDetector(sens);
      }
      f_vol.setAttributes(desc, x_fiber.regionStr(), x_fiber.limitsStr(), x_fiber.visStr());

      // Fiber placement
      // Transform3D f_tr(RotationZYX(0,0,M_PI*0.5),Position(p.x(), 0, p.y()));
      // PlacedVolume fiber_phv = s_vol.placeVolume(f_vol, Position(p.x(), 0., p.y()));
      PlacedVolume fiber_phv = lfibers.placeVolume(f_vol, Position(p.x(), 0., 0.));
      fiber_phv.addPhysVolID(f_id_grid, f_grid_id + 1).addPhysVolID(f_id_fiber, f_id + 1);
    }
    Transform3D l_tr(RotationZYX(0,0,M_PI*0.5),Position(0., 0, l_pos_y));
    s_vol.placeVolume(lfibers, l_tr);
  }
}

DECLARE_DETELEMENT(athena_EcalBarrelInterlayers, create_detector)
// DECLARE_DETELEMENT(athena_EcalBarrelInterlayers, create_detector)

