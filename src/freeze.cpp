#include <iostream>

#include <tfel/tfel.hpp>
#include <parameter/parameter.hpp>

class freeze {
public:
  freeze(const std::string& parameter_filename) {
    p.read_from_file(parameter_filename);

    particle_radius         = p.get_value<double>("particle-radius");
    alumina_diffusivity     = p.get_value<double>("alumina-diffusivity-coefficient");
    electrolyte_diffusivity = p.get_value<double>("liquid-electrolyte-diffusivity-coefficient");
    t_inj                   = p.get_value<double>("alumina-initial-temperature");
    t_init                  = p.get_value<double>("electrolyte-initial-temperature");
    t_end                   = p.get_value<double>("t-end");
    r_max                   = p.get_value<double>("r-max");
    n                       = p.get_value<int>   ("space-subdivisions");
    m                       = p.get_value<int>   ("time-subdivisions");
  }

  void run() {
    heat_cartesian_coordinates();
    heat_spherical_coordinates();    
  }
  
  
private:
  parameter::collection p;
  
  double particle_radius;
  double alumina_diffusivity;
  double electrolyte_diffusivity;
  double t_inj;
  double t_init;
  double t_end;
  double r_max;
  std::size_t n;
  std::size_t m;

  
private:
  double initial_condition(const double* x) {
    if (*x < particle_radius)
      return t_inj;
    else
      return t_init;
  }

  double diffusivity(const double* x) {
    if (*x < particle_radius)
      return alumina_diffusivity;
    else
      return electrolyte_diffusivity;
  }

  static double spherical_element(const double* x) {
    return (*x) * (*x);
  }

  void heat_cartesian_coordinates() {
    using namespace std::placeholders;

    const double delta_t(t_end / m);
  
    using cell_type = cell::edge;
    using mesh_type = fe_mesh<cell_type>;
    using fe_type = cell_type::fe::lagrange_p1;
    using fes_type = finite_element_space<fe_type>;
    using element_type = fes_type::element;

    mesh_type mesh(gen_segment_mesh(0.0, r_max, n));
    //submesh<cell_type> dmesh(mesh.get_boundary_submesh());
    submesh<cell_type> dmesh(mesh);

    fes_type fes(mesh, dmesh);

    element_type u_p(projector::lagrange<fe_type>(
      std::bind(&freeze::initial_condition, this, _1), fes));
  
    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto u(a.get_trial_function());
      auto v(a.get_test_function());

      std::function<double(const double*)> D(std::bind(&freeze::diffusivity, this, _1));
      
      a += integrate<quad::edge::gauss3>((1 / delta_t) * u * v
                                         + make_expr(D) * d<1>(u) * d<1>(v), mesh);
    }

    linear_form<fes_type> f(fes);
    for (std::size_t k(0); k < m; ++k) {
      exporter::ascii<mesh_type>("solution_cc_" + std::to_string(k) + ".dat",
                                 to_mesh_vertex_data<fe_type>(u_p));

    
      auto v(f.get_test_function());
    
      f.clear();
      f += integrate<quad::edge::gauss3>((1 / delta_t)
                                         * make_expr<fe_type>(u_p) * v,
                                         mesh);

      u_p = a.solve(f);
    }
    exporter::ascii<mesh_type>("solution_cc_" + std::to_string(m) + ".dat",
                               to_mesh_vertex_data<fe_type>(u_p));
  }

  void heat_spherical_coordinates() {
    using namespace std::placeholders;
  
    const double delta_t(t_end / m);
  
    using cell_type = cell::edge;
    using mesh_type = fe_mesh<cell_type>;
    using fe_type = cell_type::fe::lagrange_p1;
    using fes_type = finite_element_space<fe_type>;
    using element_type = fes_type::element;

    mesh_type mesh(gen_segment_mesh(0.0, r_max, n));
    //submesh<cell_type> dmesh(mesh.get_boundary_submesh());
    submesh<cell_type> dmesh(mesh);

    fes_type fes(mesh, dmesh);

    element_type u_p(
      projector::lagrange<fe_type>(
        std::bind(&freeze::initial_condition, this, _1),
        fes));
  
    bilinear_form<fes_type, fes_type> a(fes, fes); {
      auto u(a.get_trial_function());
      auto v(a.get_test_function());

      std::function<double(const double*)> D(std::bind(&freeze::diffusivity, this, _1));
      a += integrate<quad::edge::gauss3>((1 / delta_t) * make_expr(spherical_element)
                                         * u * v
                                         + make_expr(spherical_element) * make_expr(D)
                                         * d<1>(u) * d<1>(v),
                                         mesh);
    }

    linear_form<fes_type> f(fes);
    for (std::size_t k(0); k < m; ++k) {
      exporter::ascii<mesh_type>("solution_sc_" + std::to_string(k) + ".dat",
                                 to_mesh_vertex_data<fe_type>(u_p));

    
      auto v(f.get_test_function());
    
      f.clear();
      f += integrate<quad::edge::gauss3>((1 / delta_t) * make_expr(spherical_element)
                                         * make_expr<fe_type>(u_p) * v,
                                         mesh);

      u_p = a.solve(f);
    }
    exporter::ascii<mesh_type>("solution_sc_" + std::to_string(m) + ".dat",
                               to_mesh_vertex_data<fe_type>(u_p));
  }
  
};


int main(int argc, char *argv[]) {
  try {
    if (argc != 2)
      std::cout << "Usage: " << argv[0] << "  <configuration file>" << std::endl;
    else {
      freeze f(argv[1]);
      f.run();
    }
  }
  catch (const std::string& e) {
    std::cout << e << std::endl;
  }
  return 0;
}
